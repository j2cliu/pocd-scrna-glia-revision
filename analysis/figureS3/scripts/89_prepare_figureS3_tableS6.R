#!/usr/bin/env Rscript

# GLIA major revision — Figure S3 and Table S6 canonical data preparation
#
# GSE289098 is a Cell Ranger v3.0.2-aggregated processed-count matrix derived
# from the same six GSE267933 libraries and cells.  This script therefore
# performs a within-cohort processing sensitivity analysis, not replication.
# It holds fixed the 7,371 Scrublet-negative cells, seven selected transcripts,
# animal/library unit, normalization formula, contrast, and leave-one-animal-
# out definitions used for Figure 2.

suppressPackageStartupMessages({
  library(DESeq2)
  library(digest)
  library(dplyr)
  library(jsonlite)
  library(readr)
  library(tibble)
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

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx)) {
    if (idx == length(args)) stop("Missing value after ", flag)
    return(args[[idx + 1L]])
  }
  default
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
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
  if (!is.finite(pooled_variance) || pooled_variance <= 0) return(NA_real_)
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
  payload,
  outcome_type,
  outcome,
  dropped_animal = "(none)"
) {
  local <- tibble(sample = names(values), value = as.numeric(values)) |>
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
    payload,
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

prefix_metrics <- function(data, prefix, keys) {
  metric_columns <- setdiff(names(data), keys)
  rename_with(data, ~paste0(prefix, .x), all_of(metric_columns))
}

read_count_matrix <- function(path, samples) {
  data <- read_csv(path, show_col_types = FALSE)
  assert_true(names(data)[[1L]] == "gene", paste("Unexpected count schema:", path))
  assert_true(
    identical(names(data)[-1L], samples),
    paste("Count columns are not ordered C1-C3/S1-S3:", path)
  )
  values <- as.matrix(data[, samples])
  assert_true(all(is.finite(values)), paste("Non-finite counts:", path))
  assert_true(all(values >= 0), paste("Negative counts:", path))
  assert_true(max(abs(values - round(values))) == 0, paste("Non-integer counts:", path))
  storage.mode(values) <- "integer"
  rownames(values) <- data$gene
  values
}

fit_deseq2 <- function(counts, sample_meta, payload, annotation) {
  col_data <- data.frame(
    group = factor(
      sample_meta$source_group,
      levels = c("Control", "Surgery")
    ),
    row.names = sample_meta$sample
  )
  dds <- DESeqDataSetFromMatrix(
    countData = counts[, sample_meta$sample, drop = FALSE],
    colData = col_data,
    design = ~group
  )
  dds <- DESeq(dds, quiet = TRUE)
  fit <- results(
    dds,
    contrast = c("group", "Surgery", "Control")
  ) |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    as_tibble() |>
    left_join(annotation, by = c("ensembl_id" = "gene")) |>
    transmute(
      payload,
      ensembl_id,
      gene = symbol,
      base_mean = baseMean,
      log2_fold_change = log2FoldChange,
      lfc_se = lfcSE,
      ci95_low = log2FoldChange - 1.96 * lfcSE,
      ci95_high = log2FoldChange + 1.96 * lfcSE,
      wald_statistic = stat,
      p_value = pvalue,
      adjusted_p_value = padj
    )
  list(
    fit = fit,
    size_factors = tibble(
      payload,
      sample = names(sizeFactors(dds)),
      size_factor = unname(sizeFactors(dds))
    )
  )
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
default_figure_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
figure_root <- normalizePath(
  get_arg("--figure-root", default_figure_root),
  mustWork = TRUE
)

preflight_dir <- file.path(figure_root, "data", "preflight")
panel_dir <- file.path(figure_root, "data", "panel_ready")
table_dir <- file.path(figure_root, "outputs", "tableS6")
manifest_dir <- file.path(figure_root, "manifests")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

figure2_root <- file.path(dirname(figure_root), "figure2")
figure2_panel_dir <- file.path(figure2_root, "data", "panel_ready")
figure2_table_dir <- file.path(figure2_root, "outputs", "tableS3")

input_paths <- c(
  primary_counts = file.path(
    preflight_dir,
    "gse267933_primary_commoncell_pseudobulk_counts.csv.gz"
  ),
  alternative_counts = file.path(
    preflight_dir,
    "gse289098_integrated_commoncell_pseudobulk_counts.csv.gz"
  ),
  annotation = file.path(preflight_dir, "gse289098_gene_annotation.csv"),
  selected_values = file.path(
    preflight_dir,
    "gse289098_selected_gene_payload_values.csv"
  ),
  panel_values = file.path(preflight_dir, "gse289098_panel_animal_values.csv"),
  common_count_audit = file.path(
    preflight_dir,
    "gse289098_common_cell_count_audit.csv"
  ),
  library_mapping = file.path(preflight_dir, "gse289098_library_mapping.csv"),
  barcode_audit = file.path(preflight_dir, "gse289098_barcode_audit.csv"),
  library_payload_audit = file.path(
    preflight_dir,
    "gse289098_payload_library_audit.csv"
  ),
  global_payload_audit = file.path(
    preflight_dir,
    "gse289098_payload_global_audit.csv"
  ),
  preflight_gate = file.path(preflight_dir, "gse289098_gate_table.csv"),
  fig02_gene_values = file.path(figure2_panel_dir, "fig02_gene_values.csv"),
  fig02_panel_values = file.path(figure2_panel_dir, "fig02_panel_values.csv"),
  fig02_panel_effect = file.path(figure2_panel_dir, "fig02_panel_effect.csv"),
  fig02_panel_loo = file.path(figure2_panel_dir, "fig02_panel_loo.csv"),
  fig02_deseq2 = file.path(figure2_panel_dir, "fig02_deseq2.csv"),
  table_s3_effects = file.path(
    figure2_table_dir,
    "TableS3_direct_and_deseq2_effects.csv"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing Figure S3 input(s):", paste(missing_inputs, collapse = ", "))
)

genes <- c("Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3")
samples <- c("C1", "C2", "C3", "S1", "S2", "S3")
payloads <- c("primary", "gse289098")
source_group_by_sample <- c(
  C1 = "Control", C2 = "Control", C3 = "Control",
  S1 = "Surgery", S2 = "Surgery", S3 = "Surgery"
)
display_group_by_source <- c(
  Control = "Oxygen control",
  Surgery = "Combined exposure"
)
payload_display <- c(
  primary = "GSE267933 primary counts",
  gse289098 = "GSE289098 aggregated counts"
)

selected_values <- read_csv(
  input_paths[["selected_values"]],
  show_col_types = FALSE
)
panel_source <- read_csv(input_paths[["panel_values"]], show_col_types = FALSE)
common_counts <- read_csv(
  input_paths[["common_count_audit"]],
  show_col_types = FALSE
)
library_mapping <- read_csv(
  input_paths[["library_mapping"]],
  show_col_types = FALSE
)
barcode_audit <- read_csv(input_paths[["barcode_audit"]], show_col_types = FALSE)
library_payload <- read_csv(
  input_paths[["library_payload_audit"]],
  show_col_types = FALSE
)
global_payload <- read_csv(
  input_paths[["global_payload_audit"]],
  show_col_types = FALSE
)
preflight_gate <- read_csv(input_paths[["preflight_gate"]], show_col_types = FALSE)
annotation <- read_csv(input_paths[["annotation"]], show_col_types = FALSE)

assert_true(nrow(selected_values) == 42L, "Expected seven genes x six animals.")
assert_true(
  nrow(distinct(selected_values, sample, gene)) == 42L &&
    setequal(selected_values$sample, samples) &&
    setequal(selected_values$gene, genes),
  "Selected-gene animal rows are incomplete or duplicated."
)
assert_true(
  all(selected_values$n_features_summed == 1L),
  "A selected symbol maps to multiple feature rows."
)
assert_true(
  nrow(panel_source) == 6L && identical(panel_source$sample, samples),
  "Panel animal rows are incomplete or out of order."
)
assert_true(
  nrow(common_counts) == 6L &&
    identical(common_counts$sample, samples) &&
    sum(common_counts$n_common_scrublet_negative_cells) == 7371L,
  "Common-cell audit does not preserve the exact 7,371-cell set."
)
assert_true(
  nrow(global_payload) == 1L &&
    global_payload$n_common_cells_all == 20684L &&
    global_payload$n_features_primary == 27998L &&
    global_payload$n_features_integrated == 27998L &&
    global_payload$n_common_scrublet_negative_microglia == 7371L &&
    global_payload$n_entries_integrated_greater == 0L &&
    global_payload$umi_removed == 3100300L &&
    global_payload$target_umi_removed == 1205287L &&
    global_payload$target_different_matrix_entries == 920379L,
  "Global payload audit differs from the frozen identity/count anchors."
)

sample_meta <- common_counts |>
  transmute(
    sample,
    source_group = group,
    display_group = unname(display_group_by_source[group]),
    n_common_cells = n_common_scrublet_negative_cells,
    primary_total_umi,
    alternative_total_umi = integrated_total_umi
  )
assert_true(
  all(sample_meta$source_group == unname(source_group_by_sample[samples])),
  "Source group labels do not match C1-C3/S1-S3."
)

recomputed_primary <- log2(
  (selected_values$primary_gene_umi + 0.5) /
    (selected_values$primary_total_umi + 1) * 1e6
)
recomputed_alternative <- log2(
  (selected_values$integrated_gene_umi + 0.5) /
    (selected_values$integrated_total_umi + 1) * 1e6
)
assert_true(
  max(abs(recomputed_primary - selected_values$primary_log2_cpm)) < 1e-12 &&
    max(abs(recomputed_alternative - selected_values$integrated_log2_cpm)) < 1e-12,
  "Gene log2-CPM values do not match the frozen pseudocount formula."
)

gene_values_long <- bind_rows(
  selected_values |>
    transmute(
      payload = "primary",
      sample,
      source_group = group,
      gene,
      gene_umi = primary_gene_umi,
      total_umi = primary_total_umi,
      log2_cpm = primary_log2_cpm
    ),
  selected_values |>
    transmute(
      payload = "gse289098",
      sample,
      source_group = group,
      gene,
      gene_umi = integrated_gene_umi,
      total_umi = integrated_total_umi,
      log2_cpm = integrated_log2_cpm
    )
) |>
  mutate(
    display_group = unname(display_group_by_source[source_group]),
    payload_display = unname(payload_display[payload]),
    gene_order = match(gene, genes)
  ) |>
  arrange(match(payload, payloads), gene_order, match(sample, samples))

panel_values_long <- gene_values_long |>
  summarise(
    panel_mean_log2_cpm = mean(log2_cpm),
    total_umi = first(total_umi),
    .by = c(payload, payload_display, sample, source_group, display_group)
  ) |>
  arrange(match(payload, payloads), match(sample, samples))

panel_check <- panel_values_long |>
  select(payload, sample, panel_mean_log2_cpm) |>
  pivot_wider(names_from = payload, values_from = panel_mean_log2_cpm) |>
  arrange(match(sample, samples))
assert_true(
  max(abs(panel_check$primary - panel_source$primary_panel_mean_log2_cpm)) < 1e-12 &&
    max(abs(panel_check$gse289098 - panel_source$integrated_panel_mean_log2_cpm)) < 1e-12,
  "Recomputed panel values differ from the Python preflight."
)

# Primary computational positive control against frozen Figure 2/Table S3.
fig02_gene_values <- read_csv(
  input_paths[["fig02_gene_values"]],
  show_col_types = FALSE
) |>
  select(sample, gene, figure2_log2_cpm = log2_cpm)
primary_gene_check <- gene_values_long |>
  filter(payload == "primary") |>
  left_join(fig02_gene_values, by = c("sample", "gene"))
assert_true(
  nrow(primary_gene_check) == 42L &&
    max(abs(primary_gene_check$log2_cpm - primary_gene_check$figure2_log2_cpm)) < 1e-12,
  "Primary gene values fail the Figure 2 positive control."
)
fig02_panel_values <- read_csv(
  input_paths[["fig02_panel_values"]],
  show_col_types = FALSE
)
primary_panel_check <- panel_values_long |>
  filter(payload == "primary") |>
  left_join(
    fig02_panel_values |>
      select(sample, figure2_panel_value = panel_mean_log2_cpm),
    by = "sample"
  )
assert_true(
  nrow(primary_panel_check) == 6L &&
    max(abs(
      primary_panel_check$panel_mean_log2_cpm -
        primary_panel_check$figure2_panel_value
    )) < 1e-12,
  "Primary panel values fail the Figure 2 positive control."
)

# Direct full and leave-one-animal-out effects under both count payloads.
outcome_effects <- bind_rows(lapply(payloads, function(payload_name) {
  local_gene <- gene_values_long |> filter(payload == payload_name)
  local_panel <- panel_values_long |> filter(payload == payload_name)
  outcome_values <- c(
    setNames(
      lapply(genes, function(gene_name) {
        local <- local_gene |>
          filter(gene == gene_name) |>
          arrange(match(sample, samples))
        setNames(local$log2_cpm, local$sample)
      }),
      genes
    ),
    list(
      author_selected_isg7 = setNames(
        local_panel$panel_mean_log2_cpm,
        local_panel$sample
      )
    )
  )
  bind_rows(lapply(names(outcome_values), function(outcome_name) {
    bind_rows(lapply(c("(none)", samples), function(dropped) {
      summarise_outcome(
        values = outcome_values[[outcome_name]],
        sample_meta = sample_meta,
        payload = payload_name,
        outcome_type = if (
          outcome_name == "author_selected_isg7"
        ) "panel" else "gene",
        outcome = outcome_name,
        dropped_animal = dropped
      )
    }))
  }))
})) |>
  mutate(
    payload_display = unname(payload_display[payload]),
    scenario_order = match(dropped_animal, c("(none)", samples)),
    scenario = if_else(
      dropped_animal == "(none)",
      "Full cohort",
      paste("Drop", dropped_animal)
    )
  ) |>
  arrange(
    match(payload, payloads),
    match(outcome_type, c("gene", "panel")),
    match(outcome, c(genes, "author_selected_isg7")),
    scenario_order
  )

full_effects <- outcome_effects |>
  filter(dropped_animal == "(none)")
panel_effects <- outcome_effects |>
  filter(outcome_type == "panel", outcome == "author_selected_isg7") |>
  group_by(payload) |>
  mutate(
    full_difference = mean_difference[dropped_animal == "(none)"],
    attenuation_vs_own_full_percent = 100 *
      (full_difference - mean_difference) / full_difference
  ) |>
  ungroup()

fig02_effects <- read_csv(
  input_paths[["table_s3_effects"]],
  show_col_types = FALSE
) |>
  select(
    outcome_type,
    outcome,
    figure2_difference = mean_difference,
    figure2_ci_low = mean_difference_ci95_low,
    figure2_ci_high = mean_difference_ci95_high,
    figure2_d = cohens_d,
    figure2_exact_p = exact_permutation_p
  )
primary_full_check <- full_effects |>
  filter(payload == "primary") |>
  left_join(fig02_effects, by = c("outcome_type", "outcome"))
assert_true(
  nrow(primary_full_check) == 8L &&
    max(abs(primary_full_check$mean_difference - primary_full_check$figure2_difference)) < 1e-12 &&
    max(abs(primary_full_check$mean_difference_ci95_low - primary_full_check$figure2_ci_low)) < 1e-12 &&
    max(abs(primary_full_check$mean_difference_ci95_high - primary_full_check$figure2_ci_high)) < 1e-12 &&
    max(abs(primary_full_check$cohens_d - primary_full_check$figure2_d)) < 1e-12 &&
    max(abs(primary_full_check$exact_permutation_p - primary_full_check$figure2_exact_p)) < 1e-12,
  "Primary full effects fail the Figure 2/Table S3 positive control."
)

fig02_panel_effect_source <- read_csv(
  input_paths[["fig02_panel_effect"]],
  show_col_types = FALSE
)
primary_panel_full_check <- full_effects |>
  filter(
    payload == "primary",
    outcome_type == "panel",
    outcome == "author_selected_isg7"
  )
assert_true(
  nrow(fig02_panel_effect_source) == 1L &&
    nrow(primary_panel_full_check) == 1L &&
    abs(
      primary_panel_full_check$mean_difference -
        fig02_panel_effect_source$mean_difference
    ) < 1e-12 &&
    abs(
      primary_panel_full_check$cohens_d - fig02_panel_effect_source$cohens_d
    ) < 1e-12 &&
    abs(
      primary_panel_full_check$exact_permutation_p -
        fig02_panel_effect_source$exact_permutation_p
    ) < 1e-12,
  "Primary full panel effect fails the Figure 2 positive control."
)

fig02_panel_loo <- read_csv(
  input_paths[["fig02_panel_loo"]],
  show_col_types = FALSE
) |>
  select(
    dropped_animal,
    figure2_difference = mean_difference,
    figure2_ci_low = mean_difference_ci95_low,
    figure2_ci_high = mean_difference_ci95_high,
    figure2_d = cohens_d
  )
primary_loo_check <- panel_effects |>
  filter(payload == "primary") |>
  left_join(fig02_panel_loo, by = "dropped_animal")
assert_true(
  nrow(primary_loo_check) == 7L &&
    max(abs(primary_loo_check$mean_difference - primary_loo_check$figure2_difference)) < 1e-12 &&
    max(abs(primary_loo_check$mean_difference_ci95_low - primary_loo_check$figure2_ci_low)) < 1e-12 &&
    max(abs(primary_loo_check$mean_difference_ci95_high - primary_loo_check$figure2_ci_high)) < 1e-12 &&
    max(abs(primary_loo_check$cohens_d - primary_loo_check$figure2_d)) < 1e-12,
  "Primary panel LOO effects fail the Figure 2 positive control."
)

# DESeq2 uses one frozen primary-derived 13,926-feature universe for both
# payloads. Alternative-payload filtering is not allowed to change the model.
primary_counts <- read_count_matrix(input_paths[["primary_counts"]], samples)
alternative_counts <- read_count_matrix(
  input_paths[["alternative_counts"]],
  samples
)
assert_true(
  identical(rownames(primary_counts), rownames(alternative_counts)) &&
    nrow(primary_counts) == 13926L,
  "DESeq2 payloads do not share the frozen 13,926-feature universe."
)
assert_true(
  all(alternative_counts <= primary_counts),
  "Alternative common-cell pseudobulk has a count increase."
)
annotation_model <- annotation |>
  filter(gene %in% rownames(primary_counts))
assert_true(
  nrow(annotation_model) == 13926L && !anyDuplicated(annotation_model$gene),
  "Model annotation does not map one-to-one to the frozen universe."
)

deseq_primary <- fit_deseq2(
  primary_counts,
  sample_meta,
  "primary",
  annotation_model
)
deseq_alternative <- fit_deseq2(
  alternative_counts,
  sample_meta,
  "gse289098",
  annotation_model
)
deseq_selected_long <- bind_rows(
  deseq_primary$fit,
  deseq_alternative$fit
) |>
  filter(gene %in% genes) |>
  mutate(
    gene_order = match(gene, genes),
    payload_display = unname(payload_display[payload])
  ) |>
  arrange(match(payload, payloads), gene_order)
assert_true(
  nrow(deseq_selected_long) == 14L &&
    nrow(distinct(deseq_selected_long, payload, gene)) == 14L,
  "Selected-gene DESeq2 results are incomplete or duplicated."
)

fig02_deseq2 <- read_csv(input_paths[["fig02_deseq2"]], show_col_types = FALSE) |>
  select(
    gene,
    figure2_lfc = log2_fold_change,
    figure2_lfc_se = lfc_se,
    figure2_p = p_value,
    figure2_padj = adjusted_p_value
  )
primary_deseq_check <- deseq_selected_long |>
  filter(payload == "primary") |>
  left_join(fig02_deseq2, by = "gene")
assert_true(
  nrow(primary_deseq_check) == 7L &&
    max(abs(primary_deseq_check$log2_fold_change - primary_deseq_check$figure2_lfc)) < 1e-10 &&
    max(abs(primary_deseq_check$lfc_se - primary_deseq_check$figure2_lfc_se)) < 1e-10 &&
    max(abs(primary_deseq_check$p_value - primary_deseq_check$figure2_p)) < 1e-10 &&
    max(abs(primary_deseq_check$adjusted_p_value - primary_deseq_check$figure2_padj), na.rm = TRUE) < 1e-10,
  "Primary DESeq2 results fail the Figure 2 positive control."
)

# Paired payload tables.
paired_animal_values <- bind_rows(
  selected_values |>
    transmute(
      outcome_type = "gene",
      outcome = gene,
      sample,
      source_group = group,
      display_group = unname(display_group_by_source[group]),
      primary_gene_umi,
      alternative_gene_umi = integrated_gene_umi,
      primary_total_umi,
      alternative_total_umi = integrated_total_umi,
      primary_value = primary_log2_cpm,
      alternative_value = integrated_log2_cpm,
      alternative_minus_primary = integrated_minus_primary_log2_cpm,
      metric = "Gene-level pseudobulk log2-CPM"
    ),
  panel_check |>
    left_join(sample_meta, by = "sample") |>
    transmute(
      outcome_type = "panel",
      outcome = "author_selected_isg7",
      sample,
      source_group,
      display_group,
      primary_gene_umi = NA_real_,
      alternative_gene_umi = NA_real_,
      primary_total_umi,
      alternative_total_umi,
      primary_value = primary,
      alternative_value = gse289098,
      alternative_minus_primary = gse289098 - primary,
      metric = "Equal-gene-weight mean log2-CPM"
    )
) |>
  arrange(
    match(outcome_type, c("gene", "panel")),
    match(outcome, c(genes, "author_selected_isg7")),
    match(sample, samples)
  )

full_metrics <- c(
  "combined_exposure_mean", "oxygen_control_mean", "mean_difference",
  "mean_difference_ci95_low", "mean_difference_ci95_high", "welch_se",
  "welch_df", "cohens_d", "exact_permutation_hits",
  "exact_permutation_allocations", "exact_permutation_p",
  "perfect_group_separation"
)
paired_full_effects <- full_effects |>
  select(payload, outcome_type, outcome, all_of(full_metrics)) |>
  pivot_wider(
    names_from = payload,
    values_from = all_of(full_metrics),
    names_glue = "{payload}_{.value}"
  ) |>
  mutate(
    alternative_minus_primary_difference =
      gse289098_mean_difference - primary_mean_difference,
    sign_concordant =
      sign(gse289098_mean_difference) == sign(primary_mean_difference)
  ) |>
  arrange(match(outcome, c(genes, "author_selected_isg7")))

panel_effects_paired <- panel_effects |>
  select(
    payload, dropped_animal, dropped_source_group, scenario, scenario_order,
    n_combined_exposure, n_oxygen_control, mean_difference,
    mean_difference_ci95_low, mean_difference_ci95_high, cohens_d,
    exact_permutation_p, perfect_group_separation,
    attenuation_vs_own_full_percent
  ) |>
  pivot_wider(
    names_from = payload,
    values_from = c(
      n_combined_exposure, n_oxygen_control, mean_difference,
      mean_difference_ci95_low, mean_difference_ci95_high, cohens_d,
      exact_permutation_p, perfect_group_separation,
      attenuation_vs_own_full_percent
    ),
    names_glue = "{payload}_{.value}"
  ) |>
  mutate(
    alternative_minus_primary_difference =
      gse289098_mean_difference - primary_mean_difference,
    primary_sign = sign(primary_mean_difference),
    alternative_sign = sign(gse289098_mean_difference),
    sign_concordant = primary_sign == alternative_sign
  ) |>
  arrange(scenario_order)

deseq_metrics <- c(
  "base_mean", "log2_fold_change", "lfc_se", "ci95_low", "ci95_high",
  "wald_statistic", "p_value", "adjusted_p_value"
)
paired_deseq2 <- deseq_selected_long |>
  select(payload, ensembl_id, gene, gene_order, all_of(deseq_metrics)) |>
  pivot_wider(
    names_from = payload,
    values_from = all_of(deseq_metrics),
    names_glue = "{payload}_{.value}"
  ) |>
  mutate(
    alternative_minus_primary_log2_fold_change =
      gse289098_log2_fold_change - primary_log2_fold_change,
    sign_concordant =
      sign(gse289098_log2_fold_change) == sign(primary_log2_fold_change)
  ) |>
  arrange(gene_order)

assert_true(
  all(panel_effects_paired$sign_concordant) &&
    all(panel_effects_paired$primary_mean_difference > 0) &&
    all(panel_effects_paired$gse289098_mean_difference > 0),
  "Panel full/LOO direction is not concordant across payloads."
)
assert_true(
  all(paired_full_effects$sign_concordant),
  "At least one full gene/panel direction differs across payloads."
)

gene_effects_paired <- outcome_effects |>
  filter(outcome_type == "gene") |>
  select(
    payload, outcome, dropped_animal, dropped_source_group, scenario,
    scenario_order, n_combined_exposure, n_oxygen_control, mean_difference,
    mean_difference_ci95_low, mean_difference_ci95_high, cohens_d,
    exact_permutation_p, perfect_group_separation
  ) |>
  pivot_wider(
    names_from = payload,
    values_from = c(
      n_combined_exposure, n_oxygen_control, mean_difference,
      mean_difference_ci95_low, mean_difference_ci95_high, cohens_d,
      exact_permutation_p, perfect_group_separation
    ),
    names_glue = "{payload}_{.value}"
  ) |>
  mutate(
    alternative_minus_primary_difference =
      gse289098_mean_difference - primary_mean_difference,
    sign_concordant =
      sign(primary_mean_difference) == sign(gse289098_mean_difference)
  ) |>
  arrange(match(outcome, genes), scenario_order)
gene_sign_discordance <- gene_effects_paired |> filter(!sign_concordant)
assert_true(
  nrow(gene_sign_discordance) == 1L &&
    gene_sign_discordance$outcome == "Isg15" &&
    gene_sign_discordance$dropped_animal == "S3",
  "Gene-level sign-discordance diagnostic differs from the verified anchor."
)

full_panel_primary <- panel_effects_paired |>
  filter(dropped_animal == "(none)") |>
  pull(primary_mean_difference)
full_panel_alternative <- panel_effects_paired |>
  filter(dropped_animal == "(none)") |>
  pull(gse289098_mean_difference)
assert_true(
  abs(full_panel_primary - 1.0685736753232482) < 1e-12 &&
    abs(full_panel_alternative - 1.0982960515733415) < 1e-12,
  "Full panel sensitivity anchors are not reproduced."
)

final_gate <- preflight_gate |>
  mutate(
    status = if_else(gate == "Estimand concordance", "PASS", status),
    evidence = if_else(
      gate == "Estimand concordance",
      paste(
        "Full and all six leave-one-animal-out panel differences remain positive;",
        "alternative-minus-primary full contrast =",
        sprintf("%+.6f", full_panel_alternative - full_panel_primary)
      ),
      evidence
    )
  )

metadata <- tibble(
  figure = "Figure S3",
  table = "Table S6",
  accession = "GSE289098",
  role = "Same-cohort processed-count sensitivity; no biological replication",
  common_estimand_build_gate = "PASS",
  provenance_input_gate = "PARTIAL: two malformed nonanalytic local download stubs",
  panel_concordance_gate = "PASS",
  n_animals = 6L,
  n_animals_per_group = 3L,
  n_common_all_cells = 20684L,
  n_common_scrublet_negative_microglia = 7371L,
  n_common_features = 27998L,
  frozen_deseq2_feature_universe = 13926L,
  alternative_features_ge10 = 13831L,
  selected_panel_status = "Author-selected, post hoc, exploratory",
  contrast = "Combined exposure minus oxygen control",
  full_primary_panel_difference = full_panel_primary,
  full_alternative_panel_difference = full_panel_alternative,
  alternative_minus_primary_panel_difference =
    full_panel_alternative - full_panel_primary,
  gene_influence_note = paste(
    "One near-zero LOO gene estimate crosses sign:",
    "Isg15 after omitting S3; 48/49 gene full-plus-LOO scenarios agree"
  ),
  claim_ceiling = paste(
    "Technical processing sensitivity within the same six libraries;",
    "not replication, independent validation, or robustness across cohorts"
  )
)

table_s6a <- bind_cols(
  global_payload |>
    select(
      -malformed_local_series_matrix_stub,
      -malformed_local_raw_tar_stub
    ),
  tibble(
    common_estimand_build_gate = "PASS",
    panel_full_and_loo_direction_gate = "PASS",
    gene_full_direction_gate = "PASS",
    gene_full_plus_loo_sign_agreement = "48 of 49",
    analytic_inputs = paste(
      "Matrix, barcodes, features, and protocol files verified;",
      "checksums retained in the analysis manifest"
    ),
    interpretation = paste(
      "Same six libraries/cells; Cell Ranger v3.0.2-aggregated",
      "processed counts; no added biological replication"
    )
  )
)

table_s6b <- library_mapping |>
  left_join(
    barcode_audit |>
      select(
        sample, primary_all_cells, integrated_mapped_cells,
        submitted_microglia, scrublet_negative_microglia,
        missing_primary_to_integrated, duplicate_integrated_matches
      ),
    by = "sample"
  ) |>
  left_join(
    library_payload |>
      select(
        sample, primary_total_umi_all_cells = primary_total_umi,
        alternative_total_umi_all_cells = integrated_total_umi,
        umi_removed_all_cells = umi_removed,
        umi_removed_percent_all_cells = umi_removed_percent,
        primary_nnz_all_cells = primary_nnz,
        alternative_nnz_all_cells = integrated_nnz,
        nnz_removed_all_cells = nnz_removed,
        cells_with_any_payload_difference,
        bit_identical_library,
        median_cell_umi_ratio
      ),
    by = "sample"
  ) |>
  left_join(
    common_counts |>
      select(
        sample,
        n_common_scrublet_negative_cells,
        primary_total_umi_common_cells = primary_total_umi,
        alternative_total_umi_common_cells = integrated_total_umi,
        alternative_to_primary_umi_ratio_common_cells =
          integrated_to_primary_umi_ratio
      ),
    by = "sample"
  ) |>
  arrange(match(sample, samples))

deseq_size_factors <- bind_rows(
  deseq_primary$size_factors,
  deseq_alternative$size_factors
)

data_dictionary <- tribble(
  ~part, ~filename, ~unit, ~description,
  "A", "TableS6A_global_payload_and_gate_audit", "Count payload", paste(
    "Matrix dimensions, exact cell/feature overlap, count losses, frozen",
    "feature universes, malformed-stub note, and final sensitivity gates"
  ),
  "B", "TableS6B_library_mapping_and_counts", "Animal/library", paste(
    "GSE267933 identifiers, GSE289098 suffix mapping, all-cell and frozen",
    "microglial counts, and per-library count-payload changes"
  ),
  "C", "TableS6C_paired_animal_values", "Animal/library x outcome", paste(
    "Paired gene and seven-transcript-panel values under both count payloads;",
    "alternative-minus-primary differences are descriptive"
  ),
  "D", "TableS6D_paired_direct_full_effects", "Outcome", paste(
    "Full animal-level group means, mean difference, Welch interval,",
    "Cohen's d, and exact 3-v-3 permutation P under each payload"
  ),
  "E", "TableS6E_panel_full_and_leave_one_out", "Omission scenario", paste(
    "Matched full and six systematic leave-one-animal-out panel estimates;",
    "omission rows are influence diagnostics, not independent tests"
  ),
  "F", "TableS6F_selected_gene_deseq2", "Transcript", paste(
    "Selected-transcript DESeq2 estimates under both payloads using the same",
    "primary-derived 13,926-feature universe; padj is payload-specific"
  ),
  "G", "TableS6G_deseq2_size_factors", "Animal/library x payload",
  "DESeq2 median-ratio size factors for computational audit",
  "H", "TableS6H_gene_full_and_leave_one_out", "Transcript x omission scenario",
  paste(
    "Matched full and leave-one-animal-out direct-expression effects for all",
    "seven genes; supports the reported 48-of-49 sign agreement"
  )
)

panel_output_paths <- c(
  metadata = file.path(panel_dir, "figS03_metadata.csv"),
  library_audit = file.path(panel_dir, "figS03_library_audit.csv"),
  animal_panel_values = file.path(panel_dir, "figS03_animal_panel_values.csv"),
  gene_delta_heatmap = file.path(panel_dir, "figS03_gene_delta_heatmap.csv"),
  panel_effects_full_loo = file.path(
    panel_dir,
    "figS03_panel_effects_full_loo.csv"
  ),
  final_gate = file.path(panel_dir, "figS03_gate_table.csv")
)

table_bases <- c(
  A = "TableS6A_global_payload_and_gate_audit",
  B = "TableS6B_library_mapping_and_counts",
  C = "TableS6C_paired_animal_values",
  D = "TableS6D_paired_direct_full_effects",
  E = "TableS6E_panel_full_and_leave_one_out",
  F = "TableS6F_selected_gene_deseq2",
  G = "TableS6G_deseq2_size_factors",
  H = "TableS6H_gene_full_and_leave_one_out",
  dictionary = "TableS6_data_dictionary"
)
table_data <- list(
  A = table_s6a,
  B = table_s6b,
  C = paired_animal_values,
  D = paired_full_effects,
  E = panel_effects_paired,
  F = paired_deseq2,
  G = deseq_size_factors,
  H = gene_effects_paired,
  dictionary = data_dictionary
)

write_csv(metadata, panel_output_paths[["metadata"]])
write_csv(table_s6b, panel_output_paths[["library_audit"]])
write_csv(
  paired_animal_values |> filter(outcome_type == "panel"),
  panel_output_paths[["animal_panel_values"]]
)
write_csv(
  paired_animal_values |>
    filter(outcome_type == "gene") |>
    transmute(
      gene = outcome,
      gene_order = match(outcome, genes),
      sample,
      source_group,
      display_group,
      alternative_minus_primary_log2_cpm = alternative_minus_primary
    ),
  panel_output_paths[["gene_delta_heatmap"]]
)
write_csv(panel_effects_paired, panel_output_paths[["panel_effects_full_loo"]])
write_csv(final_gate, panel_output_paths[["final_gate"]])

table_output_paths <- character()
for (key in names(table_bases)) {
  csv_path <- file.path(table_dir, paste0(table_bases[[key]], ".csv"))
  tsv_path <- file.path(table_dir, paste0(table_bases[[key]], ".tsv"))
  write_csv(table_data[[key]], csv_path)
  write_tsv(table_data[[key]], tsv_path)
  table_output_paths[paste0(key, "_csv")] <- csv_path
  table_output_paths[paste0(key, "_tsv")] <- tsv_path
}

all_output_paths <- c(panel_output_paths, table_output_paths)
source_manifest <- tibble(
  input = names(input_paths),
  path = unname(input_paths),
  role = c(
    "Frozen primary common-cell pseudobulk counts",
    "GSE289098 common-cell pseudobulk counts",
    "Exact feature ID and symbol mapping",
    "Paired selected-gene counts and direct values",
    "Paired animal-level panel preflight values",
    "Per-animal common-cell counts and denominators",
    "GSM/BioSample/SRX/suffix provenance mapping",
    "Full barcode and retained-cell overlap audit",
    "Per-library all-cell count-payload audit",
    "Global count-payload audit",
    "Preflight gate table",
    "Frozen Figure 2 gene-value positive control",
    "Frozen Figure 2 panel-value positive control",
    "Frozen Figure 2 full-panel positive control",
    "Frozen Figure 2 panel-LOO positive control",
    "Frozen Figure 2 DESeq2 positive control",
    "Frozen Table S3 direct-effect positive control"
  ),
  biological_replication_added = FALSE,
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv(source_manifest, file.path(manifest_dir, "figS03_source_manifest.csv"))

output_manifest <- tibble(
  output = names(all_output_paths),
  relative_path = vapply(
    all_output_paths,
    function(path) sub(
      paste0("^", normalizePath(figure_root), "/?"),
      "",
      normalizePath(path)
    ),
    character(1)
  ),
  sha256 = vapply(all_output_paths, sha256_file, character(1))
)
write_csv(
  output_manifest,
  file.path(manifest_dir, "figS03_panel_table_output_manifest.csv")
)

execution <- list(
  figure = "Figure S3",
  supplementary_table = "Table S6",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  role = "same-cohort processed-count sensitivity only",
  biological_unit = "animal/library",
  groups = list(
    oxygen_control = c("C1", "C2", "C3"),
    combined_exposure = c("S1", "S2", "S3")
  ),
  common_cells = 7371L,
  common_features = 27998L,
  selected_genes = as.list(genes),
  selected_panel_status = "author-selected, post hoc, exploratory",
  direct_expression = paste(
    "log2(((gene UMI + 0.5)/(all-feature retained UMI + 1)) x 1e6);",
    "panel is the equal-weight mean of seven gene-level values"
  ),
  direct_uncertainty = "Welch-Satterthwaite 95% mean-difference interval",
  effect_size = "Cohen's d using pooled within-group SD",
  exact_test = paste(
    "two-sided absolute mean-difference permutation over all 20",
    "three-versus-three allocations; full cohort only; no +1 correction"
  ),
  loo = "six systematic one-animal omissions; descriptive influence analysis",
  deseq2 = list(
    design = "~group; Surgery versus Control",
    feature_universe = 13926L,
    universe_rule = "primary payload total count >=10, frozen for both fits",
    adjusted_p_values = paste(
      "BH values are payload-specific and are not an equivalence or",
      "concordance gate"
    )
  ),
  positive_control = "Primary rerun reproduces Figure 2/Table S3",
  gates = final_gate,
  gene_influence_note = metadata$gene_influence_note[[1L]],
  claim_ceiling = metadata$claim_ceiling[[1L]]
)
write_json(
  execution,
  file.path(manifest_dir, "figS03_data_preparation_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "figS03_data_preparation_sessionInfo.txt")
)

message("Figure S3 / Table S6 data preparation complete.")
message(
  sprintf(
    "Panel full difference: primary %.6f; GSE289098 %.6f; delta %+.6f",
    full_panel_primary,
    full_panel_alternative,
    full_panel_alternative - full_panel_primary
  )
)
message("Common-estimand build gate: PASS; panel full/LOO concordance: PASS.")
