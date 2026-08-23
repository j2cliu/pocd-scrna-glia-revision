#!/usr/bin/env Rscript

# GLIA major revision — Figure 2 and Table S3 canonical data preparation
#
# Figure 2 displays direct animal-level expression of an author-selected,
# post hoc seven-transcript interferon-responsive panel. This script separates:
#   1. direct pseudobulk log2-CPM values;
#   2. DESeq2 log2-fold-change estimates; and
#   3. direct-expression leave-one-animal-out influence diagnostics.
#
# It recomputes all direct-expression effects from the six animal-level values.
# Existing independent-audit effect files are used only as computational
# cross-checks, not as additional biological replication.
#
# Example:
# Rscript scripts/82_prepare_glia_R1_figure2_data.R \
#   --audit-root /path/to/independent_audit \
#   --figure-root /path/to/figure2

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(jsonlite)
  library(readr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Unable to resolve the current script path.")
  }
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

get_arg <- function(flag, default = NULL, required = FALSE) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx)) {
    if (idx == length(args)) {
      stop("Missing value after ", flag)
    }
    return(args[[idx + 1L]])
  }
  if (required && is.null(default)) {
    stop("Required argument not supplied: ", flag)
  }
  default
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

cohens_d <- function(exposed, control) {
  nx <- length(exposed)
  ny <- length(control)
  pooled_variance <- (
    (nx - 1) * var(exposed) + (ny - 1) * var(control)
  ) / (nx + ny - 2)
  if (!is.finite(pooled_variance) || pooled_variance <= 0) {
    return(NA_real_)
  }
  (mean(exposed) - mean(control)) / sqrt(pooled_variance)
}

welch_ci <- function(exposed, control, confidence = 0.95) {
  nx <- length(exposed)
  ny <- length(control)
  exposed_component <- var(exposed) / nx
  control_component <- var(control) / ny
  variance <- exposed_component + control_component
  difference <- mean(exposed) - mean(control)
  if (!is.finite(variance) || variance <= 0) {
    return(c(
      mean_difference = difference,
      welch_se = 0,
      welch_df = NA_real_,
      ci_low = difference,
      ci_high = difference
    ))
  }
  standard_error <- sqrt(variance)
  degrees_freedom <- variance^2 / (
    exposed_component^2 / (nx - 1) +
      control_component^2 / (ny - 1)
  )
  critical <- qt((1 + confidence) / 2, df = degrees_freedom)
  c(
    mean_difference = difference,
    welch_se = standard_error,
    welch_df = degrees_freedom,
    ci_low = difference - critical * standard_error,
    ci_high = difference + critical * standard_error
  )
}

exact_permutation_p <- function(exposed, control) {
  pooled <- c(exposed, control)
  observed <- abs(mean(exposed) - mean(control))
  allocations <- combn(seq_along(pooled), length(exposed))
  permuted <- apply(allocations, 2L, function(exposed_idx) {
    abs(mean(pooled[exposed_idx]) - mean(pooled[-exposed_idx]))
  })
  c(
    hits = sum(permuted >= observed - 1e-12),
    allocations = length(permuted),
    p = mean(permuted >= observed - 1e-12)
  )
}

summarise_outcome <- function(
  values,
  sample_meta,
  outcome_type,
  outcome,
  dropped_animal = "(none)"
) {
  local <- tibble(
    sample = names(values),
    value = as.numeric(values)
  ) |>
    left_join(sample_meta, by = "sample")
  if (dropped_animal != "(none)") {
    local <- local |> filter(sample != dropped_animal)
  }
  exposed <- local$value[local$source_group == "Surgery"]
  control <- local$value[local$source_group == "Control"]
  ci <- welch_ci(exposed, control)
  permutation <- if (dropped_animal == "(none)") {
    exact_permutation_p(exposed, control)
  } else {
    c(hits = NA_real_, allocations = NA_real_, p = NA_real_)
  }
  tibble(
    outcome_type,
    outcome,
    dropped_animal,
    dropped_source_group = if (dropped_animal == "(none)") {
      NA_character_
    } else {
      sample_meta$source_group[sample_meta$sample == dropped_animal]
    },
    n_combined_exposure = length(exposed),
    n_oxygen_control = length(control),
    combined_exposure_mean = mean(exposed),
    oxygen_control_mean = mean(control),
    mean_difference = unname(ci[["mean_difference"]]),
    mean_difference_ci95_low = unname(ci[["ci_low"]]),
    mean_difference_ci95_high = unname(ci[["ci_high"]]),
    welch_se = unname(ci[["welch_se"]]),
    welch_df = unname(ci[["welch_df"]]),
    cohens_d = cohens_d(exposed, control),
    exact_permutation_hits = as.integer(permutation[["hits"]]),
    exact_permutation_allocations = as.integer(permutation[["allocations"]]),
    exact_permutation_p = unname(permutation[["p"]]),
    perfect_group_separation =
      min(exposed) > max(control) || max(exposed) < min(control)
  )
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
default_figure_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

audit_root <- normalizePath(
  get_arg("--audit-root", required = TRUE),
  mustWork = TRUE
)
figure_root <- normalizePath(
  get_arg("--figure-root", default_figure_root),
  mustWork = TRUE
)

panel_dir <- file.path(figure_root, "data", "panel_ready")
table_dir <- file.path(figure_root, "outputs", "tableS3")
manifest_dir <- file.path(figure_root, "manifests")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

audit_source <- function(filename) {
  path <- file.path(audit_root, filename)
  assert_true(file.exists(path), paste("Missing audit source:", filename))
  path
}

source_paths <- c(
  animal_gene_values = audit_source(
    "scrublet_negative_selected_panel_animal_values.csv"
  ),
  deseq2_full = audit_source("scrublet_negative_deseq2_full.csv"),
  count_audit = audit_source("scrublet_negative_count_audit.csv"),
  analysis_manifest = audit_source("scrublet_negative_analysis_manifest.csv"),
  full_effect_crosscheck = audit_source(
    "doublet_negative_gene_and_panel_effects.csv"
  ),
  loo_crosscheck = audit_source(
    "doublet_negative_gene_and_panel_loo.csv"
  ),
  panel_loo_crosscheck = audit_source(
    "scrublet_negative_selected_panel_effect_loo.csv"
  )
)

genes <- c("Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3")
n_panel_genes <- length(genes)
samples <- c("C1", "C2", "C3", "S1", "S2", "S3")
source_group_by_sample <- c(
  C1 = "Control", C2 = "Control", C3 = "Control",
  S1 = "Surgery", S2 = "Surgery", S3 = "Surgery"
)
display_group_by_source <- c(
  Control = "Oxygen control",
  Surgery = "Combined exposure"
)

panel_metadata <- tibble(
  panel_id = "author_selected_isg7",
  display_name = "Seven selected interferon-responsive transcripts",
  genes = paste(genes, collapse = ";"),
  n_genes = n_panel_genes,
  selection_status = "Author-selected, post hoc, exploratory",
  direct_expression_estimator = paste(
    "Equal-weight mean of gene-level pseudobulk log2-CPM values;",
    "log2(((gene UMI + 0.5)/(total retained UMI + 1)) x 1e6)"
  ),
  analysis_set = "7,371 Scrublet-negative frozen submitted microglia",
  biological_unit = "Animal/library",
  contrast = "Combined exposure minus oxygen control"
)

# ---- Animal-level direct expression ---------------------------------------

gene_values_source <- read_csv(
  source_paths[["animal_gene_values"]],
  show_col_types = FALSE
)
required_gene_value_columns <- c(
  "symbol", "sample", "gene_umi", "group", "total_retained_umi", "log2_cpm"
)
assert_true(
  identical(names(gene_values_source), required_gene_value_columns),
  "Unexpected selected-panel animal-value columns."
)
assert_true(nrow(gene_values_source) == 42L, "Expected seven genes x six animals.")
assert_true(
  setequal(gene_values_source$symbol, genes),
  "Selected-panel genes differ from the frozen seven-transcript list."
)
assert_true(
  setequal(gene_values_source$sample, samples),
  "Selected-panel samples differ from C1-C3/S1-S3."
)
assert_true(
  !anyDuplicated(gene_values_source[c("symbol", "sample")]),
  "Gene-by-animal rows are not unique."
)
assert_true(
  all(gene_values_source$group ==
    unname(source_group_by_sample[gene_values_source$sample])),
  "Gene-value source group labels do not match the frozen design."
)
recomputed_log2_cpm <- log2(
  (gene_values_source$gene_umi + 0.5) /
    (gene_values_source$total_retained_umi + 1) * 1e6
)
assert_true(
  max(abs(gene_values_source$log2_cpm - recomputed_log2_cpm)) < 1e-12,
  "Gene-level log2-CPM values do not match the stated pseudocount formula."
)

count_audit <- read_csv(source_paths[["count_audit"]], show_col_types = FALSE)
assert_true(sum(count_audit$cells_retained) == 7371L, "Expected 7,371 retained cells.")
assert_true(
  identical(count_audit$sample, samples),
  "Count-audit samples are not in the expected order."
)
umi_map <- setNames(count_audit$umi_retained, count_audit$sample)
assert_true(
  all(gene_values_source$total_retained_umi ==
    unname(umi_map[gene_values_source$sample])),
  "Selected-panel denominators differ from the count audit."
)

sample_meta <- tibble(
  sample = samples,
  source_group = unname(source_group_by_sample[samples]),
  display_group = unname(
    display_group_by_source[unname(source_group_by_sample[samples])]
  ),
  cells_retained = count_audit$cells_retained,
  total_retained_umi = count_audit$umi_retained
)

gene_values <- gene_values_source |>
  transmute(
    sample,
    source_group = group,
    display_group = unname(display_group_by_source[group]),
    gene = symbol,
    gene_order = match(symbol, genes),
    gene_umi = as.integer(gene_umi),
    total_retained_umi = as.integer(total_retained_umi),
    log2_cpm
  ) |>
  arrange(gene_order, match(sample, samples))

panel_values <- gene_values |>
  summarise(
    panel_mean_log2_cpm = mean(log2_cpm),
    .by = c(sample, source_group, display_group)
  ) |>
  arrange(match(sample, samples)) |>
  mutate(
    panel_id = "author_selected_isg7",
    n_genes = n_panel_genes,
    genes = paste(genes, collapse = ";")
  )

# ---- Direct effects and leave-one-animal-out diagnostics ------------------

outcome_values <- c(
  setNames(
    lapply(genes, function(gene_name) {
      values <- gene_values |>
        filter(gene == gene_name) |>
        arrange(match(sample, samples))
      setNames(values$log2_cpm, values$sample)
    }),
    genes
  ),
  list(
    author_selected_isg7 = setNames(
      panel_values$panel_mean_log2_cpm,
      panel_values$sample
    )
  )
)

direct_effects_all <- bind_rows(lapply(names(outcome_values), function(outcome) {
  outcome_type <- if (outcome == "author_selected_isg7") "panel" else "gene"
  bind_rows(lapply(c("(none)", samples), function(dropped) {
    summarise_outcome(
      values = outcome_values[[outcome]],
      sample_meta = sample_meta,
      outcome_type = outcome_type,
      outcome = outcome,
      dropped_animal = dropped
    )
  }))
})) |>
  mutate(
    scenario_order = match(dropped_animal, c("(none)", samples))
  ) |>
  arrange(outcome_type, match(outcome, c(genes, "author_selected_isg7")), scenario_order)

# Join each outcome to its full-cohort effect so every sign comparison is
# explicit and auditable.
full_signs <- direct_effects_all |>
  filter(dropped_animal == "(none)") |>
  select(outcome_type, outcome, full_difference = mean_difference)
direct_effects_all <- direct_effects_all |>
  left_join(full_signs, by = c("outcome_type", "outcome")) |>
  mutate(sign_matches_full = sign(mean_difference) == sign(full_difference)) |>
  select(-full_difference)

direct_effects_full <- direct_effects_all |>
  filter(dropped_animal == "(none)") |>
  select(-scenario_order)

panel_effect <- direct_effects_full |>
  filter(outcome_type == "panel", outcome == "author_selected_isg7")
panel_loo <- direct_effects_all |>
  filter(outcome_type == "panel", outcome == "author_selected_isg7") |>
  mutate(
    scenario = if_else(
      dropped_animal == "(none)",
      "Full cohort",
      paste("Drop", dropped_animal)
    ),
    attenuation_vs_full_percent = 100 * (
      panel_effect$mean_difference[[1L]] - mean_difference
    ) / panel_effect$mean_difference[[1L]]
  ) |>
  select(-scenario_order)

assert_true(
  abs(panel_effect$mean_difference - 1.0685736753232482) < 1e-12,
  "Full panel mean difference differs from the verified anchor."
)
assert_true(
  abs(panel_effect$cohens_d - 1.3247486233261767) < 1e-12,
  "Full panel Cohen's d differs from the verified anchor."
)
assert_true(
  abs(panel_effect$exact_permutation_p - 0.20) < 1e-12,
  "Full panel exact permutation p differs from 0.20."
)
drop_s3 <- panel_loo |> filter(dropped_animal == "S3")
assert_true(
  abs(drop_s3$mean_difference - 0.449272896448337) < 1e-12,
  "Drop-S3 panel difference differs from the verified anchor."
)
assert_true(
  abs(drop_s3$attenuation_vs_full_percent - 57.95583338580468) < 1e-9,
  "Drop-S3 attenuation differs from approximately 58%."
)
assert_true(
  all(panel_loo$mean_difference > 0),
  "At least one selected-panel leave-one-animal-out difference is not positive."
)

# ---- DESeq2 estimates ------------------------------------------------------

deseq2_full <- read_csv(source_paths[["deseq2_full"]], show_col_types = FALSE)
selected_deseq2 <- deseq2_full |>
  filter(symbol %in% genes)
assert_true(nrow(selected_deseq2) == 7L, "Expected one DESeq2 row per selected gene.")
assert_true(
  !anyDuplicated(selected_deseq2$symbol),
  "Selected gene symbols are duplicated in the DESeq2 display source."
)
selected_deseq2 <- selected_deseq2 |>
  mutate(gene_order = match(symbol, genes)) |>
  arrange(gene_order) |>
  transmute(
    ensembl_id = .data$gene,
    gene = symbol,
    gene_order,
    base_mean = baseMean,
    log2_fold_change = log2FoldChange,
    lfc_se = lfcSE,
    ci95_low,
    ci95_high,
    wald_statistic = stat,
    p_value = pvalue,
    adjusted_p_value = padj,
    contrast = "Combined exposure minus oxygen control",
    analysis_set = "7,371 Scrublet-negative cells; animal-level pseudobulk"
  )

expected_lfc <- c(
  Irf7 = 1.6342935867892807,
  Ifitm3 = 1.5666930166121293,
  Isg15 = 0.845979139688111,
  Mx1 = 1.3732402685906129,
  Ifit1 = 1.6811334650296175,
  Ifit2 = 1.222021945783385,
  Ifit3 = 1.1099937621902867
)
assert_true(
  max(abs(
    selected_deseq2$log2_fold_change -
      unname(expected_lfc[selected_deseq2$gene])
  )) < 1e-12,
  "Selected-gene DESeq2 estimates differ from the verified anchors."
)
assert_true(
  all(selected_deseq2$log2_fold_change > 0),
  "At least one selected-gene full-cohort DESeq2 estimate is not positive."
)

# ---- Computational cross-checks ------------------------------------------

full_crosscheck <- read_csv(
  source_paths[["full_effect_crosscheck"]],
  show_col_types = FALSE
) |>
  filter(
    (outcome_type == "gene" & outcome %in% genes) |
      (outcome_type == "panel" & outcome == "author_selected_isg7")
  ) |>
  select(
    outcome_type,
    outcome,
    cross_mean_difference = mean_difference,
    cross_ci_low = mean_difference_ci_low,
    cross_ci_high = mean_difference_ci_high,
    cross_cohens_d = cohens_d_point_estimate,
    cross_exact_p = exact_perm_p
  )
full_check <- direct_effects_full |>
  inner_join(full_crosscheck, by = c("outcome_type", "outcome"))
assert_true(nrow(full_check) == 8L, "Full-effect cross-check did not match eight outcomes.")
assert_true(
  max(abs(full_check$mean_difference - full_check$cross_mean_difference)) < 1e-12 &&
    max(abs(
      full_check$mean_difference_ci95_low - full_check$cross_ci_low
    )) < 1e-12 &&
    max(abs(
      full_check$mean_difference_ci95_high - full_check$cross_ci_high
    )) < 1e-12 &&
    max(abs(full_check$cohens_d - full_check$cross_cohens_d)) < 1e-12 &&
    max(abs(
      full_check$exact_permutation_p - full_check$cross_exact_p
    )) < 1e-12,
  "Direct full-effect recomputation differs from the independent cross-check."
)

loo_crosscheck <- read_csv(
  source_paths[["loo_crosscheck"]],
  show_col_types = FALSE
) |>
  filter(
    (outcome_type == "gene" & outcome %in% genes) |
      (outcome_type == "panel" & outcome == "author_selected_isg7")
  ) |>
  select(
    outcome_type,
    outcome,
    dropped_animal,
    cross_mean_difference = mean_difference,
    cross_ci_low = mean_difference_ci_low,
    cross_ci_high = mean_difference_ci_high,
    cross_cohens_d = cohens_d_point_estimate
  )
loo_check <- direct_effects_all |>
  filter(dropped_animal != "(none)") |>
  inner_join(
    loo_crosscheck,
    by = c("outcome_type", "outcome", "dropped_animal")
  )
assert_true(nrow(loo_check) == 48L, "LOO cross-check did not match 8 x 6 rows.")
assert_true(
  max(abs(loo_check$mean_difference - loo_check$cross_mean_difference)) < 1e-12 &&
    max(abs(
      loo_check$mean_difference_ci95_low - loo_check$cross_ci_low
    )) < 1e-11 &&
    max(abs(
      loo_check$mean_difference_ci95_high - loo_check$cross_ci_high
    )) < 1e-11 &&
    max(abs(loo_check$cohens_d - loo_check$cross_cohens_d)) < 1e-11,
  "Direct LOO recomputation differs from the independent cross-check."
)

panel_loo_source <- read_csv(
  source_paths[["panel_loo_crosscheck"]],
  show_col_types = FALSE
)
panel_loo_check <- panel_loo |>
  inner_join(
    panel_loo_source,
    by = "dropped_animal",
    suffix = c("_recomputed", "_source")
  )
assert_true(nrow(panel_loo_check) == 7L, "Panel LOO source cross-check is incomplete.")
assert_true(
  max(abs(
    panel_loo_check$mean_difference_recomputed -
      panel_loo_check$mean_difference_source
  )) < 1e-12,
  "Panel LOO mean differences disagree with the canonical audit output."
)

# ---- Table S3 -------------------------------------------------------------

table_s3_values <- bind_rows(
  gene_values |>
    transmute(
      outcome_type = "gene",
      outcome = gene,
      sample,
      display_group,
      gene_umi,
      total_retained_umi,
      value = log2_cpm,
      metric = "Gene-level pseudobulk log2-CPM"
    ),
  panel_values |>
    transmute(
      outcome_type = "panel",
      outcome = "author_selected_isg7",
      sample,
      display_group,
      gene_umi = NA_integer_,
      total_retained_umi = unname(umi_map[sample]),
      value = panel_mean_log2_cpm,
      metric = "Equal-gene-weight mean log2-CPM"
    )
) |>
  arrange(
    match(outcome_type, c("gene", "panel")),
    match(outcome, c(genes, "author_selected_isg7")),
    match(sample, samples)
  )

table_s3_effects <- direct_effects_full |>
  select(
    outcome_type,
    outcome,
    n_combined_exposure,
    n_oxygen_control,
    combined_exposure_mean,
    oxygen_control_mean,
    mean_difference,
    mean_difference_ci95_low,
    mean_difference_ci95_high,
    cohens_d,
    exact_permutation_p,
    perfect_group_separation
  ) |>
  left_join(
    selected_deseq2 |>
      transmute(
        outcome = gene,
        deseq2_log2_fold_change = log2_fold_change,
        deseq2_ci95_low = ci95_low,
        deseq2_ci95_high = ci95_high,
        deseq2_wald_statistic = wald_statistic,
        deseq2_p_value = p_value,
        deseq2_adjusted_p_value = adjusted_p_value
      ),
    by = "outcome"
  )

table_s3_loo <- direct_effects_all |>
  mutate(
    dropped_group = if_else(
      is.na(dropped_source_group),
      NA_character_,
      unname(display_group_by_source[dropped_source_group])
    )
  ) |>
  select(
    outcome_type,
    outcome,
    dropped_animal,
    dropped_group,
    n_combined_exposure,
    n_oxygen_control,
    combined_exposure_mean,
    oxygen_control_mean,
    mean_difference,
    mean_difference_ci95_low,
    mean_difference_ci95_high,
    cohens_d,
    exact_permutation_p,
    perfect_group_separation,
    sign_matches_full
  )

# ---- Write panel-ready data and provenance --------------------------------

output_paths <- c(
  panel_metadata = file.path(panel_dir, "fig02_panel_metadata.csv"),
  sample_metadata = file.path(panel_dir, "fig02_sample_metadata.csv"),
  gene_values = file.path(panel_dir, "fig02_gene_values.csv"),
  deseq2 = file.path(panel_dir, "fig02_deseq2.csv"),
  panel_values = file.path(panel_dir, "fig02_panel_values.csv"),
  panel_effect = file.path(panel_dir, "fig02_panel_effect.csv"),
  panel_loo = file.path(panel_dir, "fig02_panel_loo.csv"),
  table_s3_values_csv = file.path(
    table_dir, "TableS3_animal_values.csv"
  ),
  table_s3_values_tsv = file.path(
    table_dir, "TableS3_animal_values.tsv"
  ),
  table_s3_effects_csv = file.path(
    table_dir, "TableS3_direct_and_deseq2_effects.csv"
  ),
  table_s3_effects_tsv = file.path(
    table_dir, "TableS3_direct_and_deseq2_effects.tsv"
  ),
  table_s3_loo_csv = file.path(
    table_dir, "TableS3_leave_one_animal_out.csv"
  ),
  table_s3_loo_tsv = file.path(
    table_dir, "TableS3_leave_one_animal_out.tsv"
  )
)

write_csv(panel_metadata, output_paths[["panel_metadata"]])
write_csv(sample_meta, output_paths[["sample_metadata"]])
write_csv(gene_values, output_paths[["gene_values"]])
write_csv(selected_deseq2, output_paths[["deseq2"]])
write_csv(panel_values, output_paths[["panel_values"]])
write_csv(panel_effect, output_paths[["panel_effect"]])
write_csv(panel_loo, output_paths[["panel_loo"]])
write_csv(table_s3_values, output_paths[["table_s3_values_csv"]])
write_tsv(table_s3_values, output_paths[["table_s3_values_tsv"]])
write_csv(table_s3_effects, output_paths[["table_s3_effects_csv"]])
write_tsv(table_s3_effects, output_paths[["table_s3_effects_tsv"]])
write_csv(table_s3_loo, output_paths[["table_s3_loo_csv"]])
write_tsv(table_s3_loo, output_paths[["table_s3_loo_tsv"]])

source_manifest <- tibble(
  source_key = names(source_paths),
  source_file = basename(source_paths),
  role = c(
    "Canonical seven-gene direct animal-level values",
    "Canonical full-cohort DESeq2 estimates",
    "Retained-cell and UMI denominator validation",
    "Transcriptome-wide software/database audit",
    "Independent full direct-effect computational cross-check",
    "Independent gene/panel LOO computational cross-check",
    "Canonical selected-panel LOO cross-check"
  ),
  biological_replication = FALSE,
  sha256 = vapply(source_paths, sha256_file, character(1))
)
write_csv(
  source_manifest,
  file.path(manifest_dir, "fig02_source_manifest.csv")
)

output_manifest <- tibble(
  output_key = names(output_paths),
  relative_output = vapply(
    output_paths,
    function(path) sub(
      paste0("^", normalizePath(figure_root), "/?"),
      "",
      normalizePath(path)
    ),
    character(1)
  ),
  sha256 = vapply(output_paths, sha256_file, character(1))
)
write_csv(
  output_manifest,
  file.path(manifest_dir, "fig02_panel_ready_manifest.csv")
)

execution_manifest <- list(
  figure = "Figure 2",
  supplementary_table = "Table S3",
  purpose = paste(
    "Animal-level expression of author-selected, post hoc",
    "interferon-responsive transcripts"
  ),
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  analysis_set_cells = 7371L,
  biological_unit = "animal/library",
  groups = list(
    oxygen_control = c("C1", "C2", "C3"),
    combined_exposure = c("S1", "S2", "S3")
  ),
  selected_genes = as.list(genes),
  selection_status = "author-selected, post hoc, exploratory",
  direct_expression = paste(
    "log2(((gene UMI + 0.5)/(total retained UMI + 1)) x 1e6);",
    "panel is equal-weight mean across seven genes"
  ),
  direct_uncertainty = "Welch-Satterthwaite 95% mean-difference interval",
  effect_size = "Cohen's d using pooled within-group SD",
  exact_test = paste(
    "two-sided absolute mean-difference permutation over all",
    "20 three-versus-three allocations; no +1 correction"
  ),
  loo_policy = paste(
    "drop each animal once; descriptive influence diagnostics;",
    "no exact p value assigned to 3-vs-2 or 2-vs-3 rows"
  ),
  deseq2 = list(
    design = "~ group",
    reference = "Control",
    contrast = "Combined exposure minus oxygen control",
    filtered_features = 13926L,
    interval = "95% Wald interval"
  ),
  random_seed = NULL,
  r_version = R.version.string,
  packages = as.list(c(
    digest = as.character(packageVersion("digest")),
    dplyr = as.character(packageVersion("dplyr")),
    jsonlite = as.character(packageVersion("jsonlite")),
    readr = as.character(packageVersion("readr")),
    tidyr = as.character(packageVersion("tidyr"))
  ))
)
write_json(
  execution_manifest,
  file.path(manifest_dir, "fig02_data_preparation_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "fig02_data_preparation_sessionInfo.txt")
)

message("Figure 2 panel-ready data written to: ", panel_dir)
message("Table S3 machine-readable files written to: ", table_dir)
message(
  "Validated panel difference: ",
  sprintf("%.3f", panel_effect$mean_difference),
  "; drop-S3: ",
  sprintf("%.3f", drop_s3$mean_difference),
  " (",
  sprintf("%.1f", drop_s3$attenuation_vs_full_percent),
  "% attenuation)."
)
