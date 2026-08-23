#!/usr/bin/env Rscript

# GLIA major revision — Figure 4 and Table S5 canonical data preparation
#
# GSE283401 is used as external between-experiment context, not as replication
# or validation of GSE267933. The 6 h and 48 h experiments contain different
# animals and were deposited as separate count matrices. Consequently, the
# time-by-treatment coefficient is compatible with a temporal difference but
# cannot identify one independently of matrix-aligned, unrecorded processing or
# sequencing-batch effects.
#
# This script:
#   1. parses the deposited SOFT metadata and audits all hippocampal libraries;
#   2. fits old-animal treatment contrasts separately at 6 h and 48 h using
#      their within-experiment >=10-count filters;
#   4. fits the old-animal ~ time + treatment + time:treatment model;
#   5. runs the same mouse MSigDB Hallmark GSEA specification for all three
#      estimands, with no additional human-to-mouse ortholog mapping;
#   6. writes panel-ready data, Table S5 Parts A–E, and provenance manifests.
#
# Example:
# Rscript scripts/86_prepare_glia_R1_figure4_data.R \
#   --project-root /path/to/pocd_scrna \
#   --figure-root /path/to/figure4

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(DESeq2)
  library(digest)
  library(dplyr)
  library(jsonlite)
  library(msigdbr)
  library(org.Mm.eg.db)
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

assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  assert_true(
    length(missing) == 0L,
    paste(label, "is missing columns:", paste(missing, collapse = ", "))
  )
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

write_both <- function(data, csv_path, tsv_path) {
  write_csv(data, csv_path)
  write_tsv(data, tsv_path)
}

extract_soft_value <- function(block, prefix, required = TRUE) {
  hit <- grep(paste0("^", prefix), block, value = TRUE)
  if (length(hit) == 0L) {
    if (required) {
      stop("SOFT block is missing: ", prefix)
    }
    return(NA_character_)
  }
  sub(paste0("^", prefix), "", hit[[1L]])
}

extract_soft_relation <- function(block, relation_label) {
  prefix <- paste0(
    "!Sample_relation = ",
    relation_label,
    ": .*"
  )
  hit <- grep(paste0("^", prefix), block, value = TRUE)
  if (length(hit) != 1L) {
    stop("Expected one ", relation_label, " relation in a SOFT block.")
  }
  sub("^.*[/=]", "", hit[[1L]])
}

parse_soft_samples <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE = ", lines)
  assert_true(length(starts) > 0L, "No sample blocks found in SOFT file.")
  ends <- c(starts[-1L] - 1L, length(lines))

  rows <- lapply(seq_along(starts), function(i) {
    block <- lines[starts[[i]]:ends[[i]]]
    characteristics <- sub(
      "^!Sample_characteristics_ch1 = ",
      "",
      grep("^!Sample_characteristics_ch1 = ", block, value = TRUE)
    )
    characteristic_keys <- sub(":.*$", "", characteristics)
    characteristic_values <- sub("^[^:]+: ?", "", characteristics)
    characteristic_map <- setNames(
      characteristic_values,
      characteristic_keys
    )

    facs_hit <- grep(
      "^!Sample_description = FACS[0-9]+_(HP|HM)$",
      block,
      value = TRUE
    )
    assert_true(
      length(facs_hit) == 1L,
      "Expected exactly one canonical FACS identifier in each SOFT block."
    )

    tibble(
      gsm = extract_soft_value(block, "!Sample_geo_accession = "),
      title = extract_soft_value(block, "!Sample_title = "),
      facs_id = sub("^!Sample_description = ", "", facs_hit[[1L]]),
      biosample = extract_soft_relation(block, "BioSample"),
      sra_experiment = extract_soft_relation(block, "SRA"),
      tissue = unname(characteristic_map[["tissue"]]),
      region = unname(characteristic_map[["region"]]),
      cell_type = unname(characteristic_map[["cell type"]]),
      genotype = unname(characteristic_map[["genotype"]]),
      age_source = unname(characteristic_map[["age"]]),
      treatment_source = unname(characteristic_map[["treatment"]]),
      time_source = unname(characteristic_map[["timepoint"]]),
      platform = extract_soft_value(block, "!Sample_platform_id = "),
      instrument = extract_soft_value(block, "!Sample_instrument_model = "),
      library_strategy = extract_soft_value(
        block,
        "!Sample_library_strategy = "
      )
    )
  })

  bind_rows(rows)
}

map_symbols <- function(ensembl_ids) {
  clean_ids <- sub("\\.\\d+$", "", ensembl_ids)
  unname(
    AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = clean_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
}

format_gene_result <- function(result_object, estimand, contrast_label) {
  result_object |>
    as.data.frame() |>
    rownames_to_column("ensembl") |>
    as_tibble() |>
    mutate(
      symbol = map_symbols(ensembl),
      estimand = estimand,
      contrast = contrast_label,
      ci95_low = log2FoldChange - qnorm(0.975) * lfcSE,
      ci95_high = log2FoldChange + qnorm(0.975) * lfcSE,
      direction = case_when(
        is.na(log2FoldChange) ~ NA_character_,
        log2FoldChange > 0 ~ "More positive in combined exposure",
        log2FoldChange < 0 ~ "More positive in source-labelled control",
        TRUE ~ "Zero"
      )
    ) |>
    relocate(
      estimand,
      contrast,
      ensembl,
      symbol,
      baseMean,
      log2FoldChange,
      lfcSE,
      ci95_low,
      ci95_high,
      stat,
      pvalue,
      padj,
      direction
    )
}

run_hallmark_gsea <- function(
  gene_result,
  hallmark,
  estimand,
  contrast_label
) {
  ranks <- gene_result |>
    filter(!is.na(symbol), !is.na(stat)) |>
    group_by(symbol) |>
    slice_max(order_by = abs(stat), n = 1L, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(stat))
  gene_list <- setNames(ranks$stat, ranks$symbol)
  assert_true(
    !anyDuplicated(names(gene_list)),
    paste("Duplicated symbols remain in", estimand, "GSEA ranks.")
  )

  set.seed(42)
  fit <- clusterProfiler::GSEA(
    gene_list,
    TERM2GENE = hallmark,
    pvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    eps = 1e-30,
    seed = TRUE,
    verbose = FALSE
  ) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(
      estimand = estimand,
      contrast = contrast_label,
      fdr_lt_0_05 = p.adjust < 0.05,
      nes_rank_descending = min_rank(desc(NES)),
      absolute_nes_rank = min_rank(desc(abs(NES)))
    ) |>
    relocate(estimand, contrast)

  assert_true(
    nrow(fit) == 50L && n_distinct(fit$ID) == 50L,
    paste("Expected all 50 Hallmark sets for", estimand, "GSEA.")
  )
  fit
}

# ---- Paths and immutable inputs -------------------------------------------

script_path <- get_script_path()
script_dir <- dirname(script_path)
default_figure_root <- normalizePath(
  file.path(script_dir, ".."),
  mustWork = TRUE
)
project_root <- normalizePath(
  get_arg("--project-root", required = TRUE),
  mustWork = TRUE
)
figure_root <- normalizePath(
  get_arg("--figure-root", default_figure_root),
  mustWork = TRUE
)

source_dir <- file.path(project_root, "data", "raw", "GSE283401")
panel_dir <- file.path(figure_root, "data", "panel_ready")
table_dir <- file.path(figure_root, "outputs", "tableS5")
manifest_dir <- file.path(figure_root, "manifests")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

source_paths <- c(
  soft = file.path(source_dir, "GSE283401_family.soft"),
  counts_6h = file.path(source_dir, "GSE283401_hp_6h_counts.csv.gz"),
  counts_48h = file.path(source_dir, "GSE283401_hp_48h_counts.csv.gz")
)
assert_true(
  all(file.exists(source_paths)),
  "One or more GSE283401 source files are missing."
)

expected_source_hashes <- c(
  soft = "9a48201c185162cd717b0463a0d39981108a19c883646e3ebf98c73b311228b0",
  counts_6h =
    "6c69ce150ce70b0b4758d002c13e736cc82558c8ac312332d9e8c935b019e363",
  counts_48h =
    "0a9cdc0554c895689217192abb1eb50a8dabfe69edb198dfdd455bb0f394470f"
)
observed_source_hashes <- vapply(
  source_paths,
  sha256_file,
  FUN.VALUE = character(1)
)
assert_true(
  identical(observed_source_hashes, expected_source_hashes),
  "A GSE283401 source checksum differs from the frozen audit."
)

# ---- SOFT metadata and count-matrix design audit --------------------------

soft_meta <- parse_soft_samples(source_paths[["soft"]])
assert_true(
  nrow(soft_meta) == 96L &&
    n_distinct(soft_meta$gsm) == 96L &&
    n_distinct(soft_meta$facs_id) == 96L,
  "Expected 96 one-to-one GSE283401 SOFT sample records."
)

hippocampus_meta <- soft_meta |>
  filter(region == "hippocampus") |>
  mutate(
    age = case_when(
      age_source == "3-5 months" ~ "Young (3–5 months)",
      age_source == "20-22 months" ~ "Old (20–22 months)",
      TRUE ~ NA_character_
    ),
    age_order = if_else(age_source == "3-5 months", 1L, 2L),
    treatment = case_when(
      treatment_source == "control" ~ "Source-labelled control",
      treatment_source == "surgery" ~ "Isoflurane + laparotomy",
      TRUE ~ NA_character_
    ),
    treatment_order = if_else(treatment_source == "control", 1L, 2L),
    time = case_when(
      time_source == "6 hours" ~ "6 h",
      time_source == "48 hours" ~ "48 h",
      TRUE ~ NA_character_
    ),
    time_order = if_else(time_source == "6 hours", 1L, 2L),
    facs_core = sub("_(HP|HM)$", "", facs_id),
    deposited_matrix = if_else(
      time == "6 h",
      "GSE283401_hp_6h_counts.csv.gz",
      "GSE283401_hp_48h_counts.csv.gz"
    ),
    included_in_old_model = age_source == "20-22 months",
    biological_unit = "Animal/library"
  ) |>
  arrange(time_order, age_order, treatment_order, facs_id)

assert_true(
  nrow(hippocampus_meta) == 48L &&
    all(!is.na(hippocampus_meta$age)) &&
    all(!is.na(hippocampus_meta$treatment)) &&
    all(!is.na(hippocampus_meta$time)),
  "The 48 hippocampal SOFT records were not completely classified."
)
assert_true(
  n_distinct(hippocampus_meta$facs_core) == 48L,
  "A hippocampal FACS identifier appears at more than one time point."
)

counts_6h_source <- read_csv(
  source_paths[["counts_6h"]],
  show_col_types = FALSE
)
counts_48h_source <- read_csv(
  source_paths[["counts_48h"]],
  show_col_types = FALSE
)
assert_columns(counts_6h_source, "gene", "6 h count matrix")
assert_columns(counts_48h_source, "gene", "48 h count matrix")
assert_true(
  nrow(counts_6h_source) == 54532L &&
    nrow(counts_48h_source) == 54532L &&
    identical(counts_6h_source$gene, counts_48h_source$gene) &&
    !anyDuplicated(counts_6h_source$gene),
  "The two hippocampal count matrices do not share the frozen gene rows."
)

count_columns_6h <- setdiff(names(counts_6h_source), "gene")
count_columns_48h <- setdiff(names(counts_48h_source), "gene")
assert_true(
  length(count_columns_6h) == 23L &&
    length(count_columns_48h) == 25L &&
    !length(intersect(count_columns_6h, count_columns_48h)) &&
    setequal(
      c(count_columns_6h, count_columns_48h),
      hippocampus_meta$facs_id
    ),
  "Count-matrix columns do not map one-to-one to the 48 hippocampal samples."
)

counts_6h <- as.matrix(counts_6h_source[, count_columns_6h])
counts_48h <- as.matrix(counts_48h_source[, count_columns_48h])
rownames(counts_6h) <- counts_6h_source$gene
rownames(counts_48h) <- counts_48h_source$gene
storage.mode(counts_6h) <- "integer"
storage.mode(counts_48h) <- "integer"
assert_true(
  all(is.finite(counts_6h)) &&
    all(is.finite(counts_48h)) &&
    all(counts_6h >= 0L) &&
    all(counts_48h >= 0L),
  "Count matrices contain non-finite or negative values."
)

library_totals <- c(colSums(counts_6h), colSums(counts_48h))
hippocampus_meta <- hippocampus_meta |>
  mutate(
    matrix_order = match(
      facs_id,
      c(count_columns_6h, count_columns_48h)
    ),
    library_total = unname(library_totals[facs_id])
  ) |>
  arrange(matrix_order)
assert_true(
  all(!is.na(hippocampus_meta$matrix_order)) &&
    all(!is.na(hippocampus_meta$library_total)),
  "Library totals were not recovered for every hippocampal sample."
)

design_counts <- hippocampus_meta |>
  count(
    time,
    time_order,
    age,
    age_order,
    treatment,
    treatment_order,
    name = "n_animals"
  ) |>
  arrange(time_order, age_order, treatment_order)
expected_design <- tribble(
  ~time, ~age, ~treatment, ~n_animals,
  "6 h", "Young (3–5 months)", "Source-labelled control", 4L,
  "6 h", "Young (3–5 months)", "Isoflurane + laparotomy", 4L,
  "6 h", "Old (20–22 months)", "Source-labelled control", 8L,
  "6 h", "Old (20–22 months)", "Isoflurane + laparotomy", 7L,
  "48 h", "Young (3–5 months)", "Source-labelled control", 4L,
  "48 h", "Young (3–5 months)", "Isoflurane + laparotomy", 4L,
  "48 h", "Old (20–22 months)", "Source-labelled control", 8L,
  "48 h", "Old (20–22 months)", "Isoflurane + laparotomy", 9L
)
assert_true(
  identical(
    design_counts |>
      select(time, age, treatment, n_animals),
    expected_design
  ),
  "The recovered age-by-treatment-by-time design differs from the audit."
)

old_meta <- hippocampus_meta |>
  filter(included_in_old_model) |>
  mutate(
    time_factor = factor(time, levels = c("6 h", "48 h")),
    treatment_factor = factor(
      treatment_source,
      levels = c("control", "surgery")
    )
  ) |>
  as.data.frame()
rownames(old_meta) <- old_meta$facs_id
assert_true(
  nrow(old_meta) == 32L &&
    sum(old_meta$time == "6 h" & old_meta$treatment_source == "control") ==
      8L &&
    sum(old_meta$time == "6 h" & old_meta$treatment_source == "surgery") ==
      7L &&
    sum(old_meta$time == "48 h" & old_meta$treatment_source == "control") ==
      8L &&
    sum(old_meta$time == "48 h" & old_meta$treatment_source == "surgery") ==
      9L,
  "The old-animal analysis design is not 8/7 at 6 h and 8/9 at 48 h."
)

# ---- Feature universes ----------------------------------------------------

old_counts <- cbind(
  counts_6h[, old_meta$facs_id[old_meta$time == "6 h"], drop = FALSE],
  counts_48h[, old_meta$facs_id[old_meta$time == "48 h"], drop = FALSE]
)
old_counts <- old_counts[, old_meta$facs_id, drop = FALSE]
interaction_keep <- rowSums(old_counts) >= 10L
interaction_features <- rownames(old_counts)[interaction_keep]
old_counts_filtered <- old_counts[interaction_features, , drop = FALSE]
assert_true(
  length(interaction_features) == 30280L,
  "The combined old-animal >=10-count feature universe is not 30,280."
)

# ---- DESeq2: separate timepoint contrasts and formal interaction ----------

fit_timepoint <- function(time_label) {
  local_meta <- old_meta[old_meta$time == time_label, , drop = FALSE]
  local_counts_source <- if (time_label == "6 h") counts_6h else counts_48h
  local_counts <- local_counts_source[, rownames(local_meta), drop = FALSE]
  local_keep <- rowSums(local_counts) >= 10L
  local_counts <- local_counts[local_keep, , drop = FALSE]
  col_data <- data.frame(
    treatment = factor(
      local_meta$treatment_source,
      levels = c("control", "surgery")
    ),
    row.names = rownames(local_meta)
  )
  dds <- DESeqDataSetFromMatrix(
    countData = local_counts,
    colData = col_data,
    design = ~ treatment
  )
  dds <- DESeq(dds, quiet = TRUE)
  result_object <- results(
    dds,
    contrast = c("treatment", "surgery", "control"),
    alpha = 0.05
  )
  result <- format_gene_result(
    result_object,
    estimand = paste0("treatment_", gsub(" ", "", time_label)),
    contrast_label = paste0(
      time_label,
      ": isoflurane plus laparotomy minus source-labelled control"
    )
  )
  attr(result, "n_features") <- nrow(local_counts)
  result
}

gene_6h <- fit_timepoint("6 h")
gene_48h <- fit_timepoint("48 h")
assert_true(
  nrow(gene_6h) == 26317L && nrow(gene_48h) == 23732L,
  "The separate 6 h/48 h >=10-count filters differ from 26,317/23,732."
)

interaction_col_data <- data.frame(
  time = factor(
    ifelse(old_meta$time == "6 h", "6h", "48h"),
    levels = c("6h", "48h")
  ),
  treatment = factor(
    old_meta$treatment_source,
    levels = c("control", "surgery")
  ),
  row.names = old_meta$facs_id
)
interaction_dds <- DESeqDataSetFromMatrix(
  countData = old_counts_filtered,
  colData = interaction_col_data,
  design = ~ time + treatment + time:treatment
)
interaction_dds <- DESeq(interaction_dds, quiet = TRUE)
interaction_name <- grep(
  "time48h.*treatmentsurgery",
  resultsNames(interaction_dds),
  value = TRUE
)
assert_true(
  length(interaction_name) == 1L,
  paste(
    "Could not identify the time-by-treatment coefficient:",
    paste(resultsNames(interaction_dds), collapse = ", ")
  )
)
interaction_gene <- format_gene_result(
  results(interaction_dds, name = interaction_name, alpha = 0.05),
  estimand = "time_by_treatment",
  contrast_label = paste(
    "(48 h exposure − 48 h control) −",
    "(6 h exposure − 6 h control)"
  )
) |>
  mutate(
    direction = case_when(
      is.na(log2FoldChange) ~ NA_character_,
      log2FoldChange > 0 ~
        "Treatment contrast more positive in deposited 48 h experiment",
      log2FoldChange < 0 ~
        "Treatment contrast more positive in deposited 6 h experiment",
      TRUE ~ "Zero interaction"
    )
  )

assert_true(
  nrow(interaction_gene) == length(interaction_features),
  "The interaction result does not use the 30,280-feature universe."
)

# ---- Mouse MSigDB Hallmark GSEA, hard-fail specification ------------------

hallmark_raw <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "MH"
)
assert_true(
  identical(sort(unique(hallmark_raw$db_version)), "2026.1.Mm"),
  paste(
    "Expected mouse MSigDB 2026.1.Mm; received:",
    paste(sort(unique(hallmark_raw$db_version)), collapse = ";")
  )
)
hallmark <- hallmark_raw |>
  select(gs_name, gene_symbol) |>
  distinct()
assert_true(
  n_distinct(hallmark$gs_name) == 50L &&
    nrow(hallmark) == 7191L,
  "Mouse MSigDB Hallmark membership differs from the frozen database audit."
)

hallmark_6h <- run_hallmark_gsea(
  gene_6h,
  hallmark,
  estimand = "treatment_6h",
  contrast_label =
    "6 h: isoflurane plus laparotomy minus source-labelled control"
)
hallmark_48h <- run_hallmark_gsea(
  gene_48h,
  hallmark,
  estimand = "treatment_48h",
  contrast_label =
    "48 h: isoflurane plus laparotomy minus source-labelled control"
)
hallmark_interaction <- run_hallmark_gsea(
  interaction_gene,
  hallmark,
  estimand = "time_by_treatment",
  contrast_label = paste(
    "(48 h exposure − 48 h control) −",
    "(6 h exposure − 6 h control)"
  )
)
hallmark_all_estimands <- bind_rows(
  hallmark_6h,
  hallmark_48h,
  hallmark_interaction
)

# ---- Panel-ready selected results and verified anchors --------------------

focus_genes <- c(
  "Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3"
)
focus_gene_order <- setNames(seq_along(focus_genes), focus_genes)
selected_gene_timepoint <- bind_rows(gene_6h, gene_48h) |>
  filter(symbol %in% focus_genes) |>
  group_by(estimand, symbol) |>
  slice_max(baseMean, n = 1L, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    transcript_order = unname(focus_gene_order[symbol]),
    time = if_else(estimand == "treatment_6h", "6 h", "48 h"),
    time_order = if_else(time == "6 h", 1L, 2L),
    n_control = if_else(time == "6 h", 8L, 8L),
    n_combined_exposure = if_else(time == "6 h", 7L, 9L),
    inference_status = paste(
      "Separately fitted gene-level descriptive context;",
      "not a within-animal time course"
    ),
    panel_status = "Author-selected, post hoc, exploratory"
  ) |>
  arrange(transcript_order, time_order)
assert_true(
  nrow(selected_gene_timepoint) == 14L &&
    n_distinct(selected_gene_timepoint$symbol) == 7L,
  "Expected seven selected transcripts at each of two time points."
)

selected_gene_interaction <- interaction_gene |>
  filter(symbol %in% focus_genes) |>
  group_by(symbol) |>
  slice_max(baseMean, n = 1L, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    transcript_order = unname(focus_gene_order[symbol]),
    time = "48 h versus 6 h interaction",
    time_order = 3L,
    n_control = 16L,
    n_combined_exposure = 16L,
    inference_status = paste(
      "Between-experiment interaction; not a within-animal time course"
    ),
    panel_status = "Author-selected, post hoc, exploratory"
  ) |>
  arrange(transcript_order)
selected_gene_all_estimands <- bind_rows(
  selected_gene_timepoint,
  selected_gene_interaction
)

expected_gene_anchors <- tribble(
  ~symbol, ~lfc_6h, ~lfc_48h, ~lfc_interaction,
  "Irf7", -0.5030606386528876, 0.4220254975908904, 0.9267446232509424,
  "Ifitm3", -1.029771242336015, 0.328500796416897, 1.3467311236102575,
  "Isg15", -0.5665919929859579, -0.32714856142447146, 0.22564485001934348,
  "Mx1", -0.9061491429994939, -0.2647087835342375, 0.6273423593733181,
  "Ifit1", -1.1388598144782807, 0.09734235347969465, 1.2944684599337841,
  "Ifit2", -0.3093376181163432, 1.3515671136934, 1.668016499066212,
  "Ifit3", -1.246091438011821, 0.7602468852222085, 2.003279077313509
)
observed_gene_anchors <- selected_gene_all_estimands |>
  select(symbol, estimand, log2FoldChange) |>
  pivot_wider(names_from = estimand, values_from = log2FoldChange) |>
  rename(
    observed_lfc_6h = treatment_6h,
    observed_lfc_48h = treatment_48h,
    observed_lfc_interaction = time_by_treatment
  ) |>
  left_join(expected_gene_anchors, by = "symbol")
gene_anchor_max_delta <- max(abs(c(
    observed_gene_anchors$observed_lfc_6h -
      observed_gene_anchors$lfc_6h,
    observed_gene_anchors$observed_lfc_48h -
      observed_gene_anchors$lfc_48h,
    observed_gene_anchors$observed_lfc_interaction -
      observed_gene_anchors$lfc_interaction
  )))
if (!is.finite(gene_anchor_max_delta) || gene_anchor_max_delta >= 1e-10) {
  print(as.data.frame(observed_gene_anchors), row.names = FALSE, digits = 16)
}
assert_true(
  is.finite(gene_anchor_max_delta) && gene_anchor_max_delta < 1e-10,
  "Selected-transcript DESeq2 anchors differ from the independent audit."
)

focus_ids <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
focus_display <- c(
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon alpha response",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon gamma response",
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA signaling via NF-kB",
  HALLMARK_INFLAMMATORY_RESPONSE =
    "Hallmark inflammatory-response gene set"
)
hallmark_interaction_focus <- hallmark_interaction |>
  filter(ID %in% focus_ids) |>
  mutate(
    hallmark_display = unname(focus_display[ID]),
    hallmark_order = match(ID, focus_ids),
    interpretation = paste(
      "Positive NES means the treatment-associated contrast is more positive",
      "in the deposited 48 h experiment than in the deposited 6 h experiment"
    ),
    limitation = paste(
      "Between-experiment interaction; different animals;",
      "time is inseparable from deposited matrix and any unrecorded",
      "time-aligned processing or sequencing-batch effects"
    )
  ) |>
  arrange(hallmark_order)
assert_true(
  nrow(hallmark_interaction_focus) == 4L,
  "The four focused interaction Hallmarks were not recovered."
)

anchor <- hallmark_interaction_focus |>
  select(ID, NES, p.adjust)
expected_anchor <- tibble(
  ID = focus_ids,
  expected_nes = c(
    2.624158343520711,
    2.365226167247039,
    1.3432884692591072,
    1.708914471375484
  ),
  expected_fdr = c(
    7.313829969354819e-11,
    2.357760270412666e-11,
    0.11355003162542474,
    0.001289550290487773
  )
)
anchor_check <- anchor |>
  left_join(expected_anchor, by = "ID")
assert_true(
  max(abs(anchor_check$NES - anchor_check$expected_nes)) < 1e-10 &&
    max(abs(anchor_check$p.adjust - anchor_check$expected_fdr)) < 1e-10,
  "Canonical interaction Hallmark anchors differ from the independent audit."
)

expected_hallmark_timepoint <- tribble(
  ~estimand, ~ID, ~expected_nes, ~expected_fdr,
  "treatment_6h", focus_ids[[1]], -2.394195, 1.38918e-7,
  "treatment_6h", focus_ids[[2]], -2.100494, 3.89843e-7,
  "treatment_6h", focus_ids[[3]], -1.741438, 0.00109587,
  "treatment_6h", focus_ids[[4]], -1.594133, 0.00642635,
  "treatment_48h", focus_ids[[1]], 1.688735, 0.0191733,
  "treatment_48h", focus_ids[[2]], 1.588051, 0.0191733,
  "treatment_48h", focus_ids[[3]], -1.293214, 0.282486,
  "treatment_48h", focus_ids[[4]], 1.012537, 0.715398
)
observed_hallmark_timepoint <- hallmark_all_estimands |>
  filter(
    estimand %in% c("treatment_6h", "treatment_48h"),
    ID %in% focus_ids
  ) |>
  select(estimand, ID, NES, p.adjust) |>
  left_join(expected_hallmark_timepoint, by = c("estimand", "ID"))
assert_true(
  max(abs(
    observed_hallmark_timepoint$NES -
      observed_hallmark_timepoint$expected_nes
  )) < 1e-5 &&
    max(abs(
      observed_hallmark_timepoint$p.adjust -
        observed_hallmark_timepoint$expected_fdr
    )) < 1e-5,
  "Time-specific Hallmark QA anchors differ from the independent audit."
)

study_comparison <- tribble(
  ~dataset, ~role, ~age, ~assay, ~exposure, ~sampling, ~biological_n,
  "GSE267933", "Primary cohort", "18-month male mice",
  "10x scRNA-seq",
  "2.5% sevoflurane for 30 min + laparotomy",
  "24 h", "3 control; 3 exposed",
  "GSE283401", "External context", "20–22-month male mice",
  "FACS microglial bulk RNA-seq",
  "1.2% isoflurane in 30% O2 for 2 h + laparotomy",
  "6 h from exposure start", "8 control; 7 exposed",
  "GSE283401", "External context", "20–22-month male mice",
  "FACS microglial bulk RNA-seq",
  "1.2% isoflurane in 30% O2 for 2 h + laparotomy",
  "48 h from exposure start", "8 control; 9 exposed"
) |>
  mutate(
    animal_relation = if_else(
      dataset == "GSE267933",
      "Primary six-animal cohort",
      "Different animals at 6 h and 48 h"
    ),
    claim_role = if_else(
      dataset == "GSE267933",
      "Primary exploratory association",
      "Independent context; not replication or validation"
    )
  )

figure_metadata <- tibble(
  figure_id = "Figure 4",
  title = paste(
    "External bulk-microglial experiments provide bounded context for",
    "time-aligned differences in treatment contrast"
  ),
  external_accession = "GSE283401",
  biological_unit = "Animal/library",
  old_model_n = 32L,
  old_6h_control_n = 8L,
  old_6h_exposure_n = 7L,
  old_48h_control_n = 8L,
  old_48h_exposure_n = 9L,
  design = "~ time + treatment + time:treatment",
  interaction = paste(
    "(48 h exposure − 48 h control) −",
    "(6 h exposure − 6 h control)"
  ),
  database = "Mouse MSigDB 2026.1.Mm Hallmarks",
  ranking = "Signed DESeq2 Wald statistic",
  gene_feature_universe_6h = nrow(gene_6h),
  gene_feature_universe_48h = nrow(gene_48h),
  gene_feature_universe_interaction = length(interaction_features),
  claim_ceiling = paste(
    "Between-experiment difference in treatment-associated contrast",
    "compatible with, but not identifying, a temporal difference"
  ),
  limitation = paste(
    "Different animals at 6 h and 48 h; time is inseparable from deposited",
    "matrix and any unrecorded time-aligned processing or sequencing-batch",
    "effects; not a trajectory, replication, or validation"
  )
)

# ---- Table S5 and panel-ready files ---------------------------------------

sample_manifest_table <- hippocampus_meta |>
  select(
    gsm,
    facs_id,
    facs_core,
    biosample,
    sra_experiment,
    matrix_order,
    age_source,
    age,
    treatment_source,
    treatment,
    time_source,
    time,
    region,
    cell_type,
    genotype,
    platform,
    instrument,
    library_strategy,
    deposited_matrix,
    library_total,
    included_in_old_model,
    biological_unit
  )

add_gene_table_notes <- function(
  data,
  feature_description,
  estimand_description
) {
  data |>
    mutate(
      feature_universe = feature_description,
      symbol_mapping = paste0(
        "org.Mm.eg.db ",
        as.character(packageVersion("org.Mm.eg.db")),
        "; multiVals=first"
      ),
      estimand_note = estimand_description
    )
}

table_s5_b <- add_gene_table_notes(
  gene_6h,
  "26,317 Ensembl features with >=10 total counts across 15 old 6 h animals",
  paste(
    "Old animals, 6 h: isoflurane-plus-laparotomy versus",
    "source-labelled control"
  )
)
table_s5_c <- add_gene_table_notes(
  gene_48h,
  "23,732 Ensembl features with >=10 total counts across 17 old 48 h animals",
  paste(
    "Old animals, 48 h: isoflurane-plus-laparotomy versus",
    "source-labelled control"
  )
)
table_s5_d <- add_gene_table_notes(
  interaction_gene,
  paste(
    "30,280 shared Ensembl features with >=10 total counts across all",
    "32 old animals"
  ),
  paste(
    "Positive interaction = treatment contrast more positive in the",
    "deposited 48 h than 6 h experiment; not within-animal change"
  )
)
table_s5_e <- hallmark_interaction |>
  mutate(
    database = "Mouse MSigDB 2026.1.Mm Hallmarks",
    ranking = "Signed DESeq2 Wald statistic",
    bh_family = "All 50 Hallmark sets within the interaction estimand",
    interpretation_limit = paste(
      "Between-experiment interaction; not within-animal change, trajectory,",
      "replication, validation, or causal attribution to time"
    )
  )

table_s5_dictionary <- tibble(
  part = c("A", "B", "C", "D", "E"),
  file_stem = c(
    "TableS5A_gse283401_sample_manifest",
    "TableS5B_gene_treatment_contrast_6h",
    "TableS5C_gene_treatment_contrast_48h",
    "TableS5D_gene_time_by_treatment_interaction",
    "TableS5E_hallmark_time_by_treatment_interaction"
  ),
  contents = c(
    "All 48 hippocampal source libraries; 32 old-animal model rows flagged",
    "Old-animal 6 h treatment-versus-control gene-level DESeq2 results",
    "Old-animal 48 h treatment-versus-control gene-level DESeq2 results",
    "Old-animal time-by-treatment gene-level DESeq2 interaction",
    "All 50 mouse MSigDB Hallmark interaction results"
  ),
  unit_or_scale = c(
    "Animal/library",
    "Unshrunk DESeq2 log2 fold change",
    "Unshrunk DESeq2 log2 fold change",
    "DESeq2 interaction log2 fold change",
    "GSEA normalized enrichment score"
  ),
  interpretation_limit = c(
    "Young animals document the source design but are not in the old model",
    "Separately fitted descriptive gene-level context",
    "Separately fitted descriptive gene-level context",
    paste(
      "Different animals and separate matrices; compatible with but does not",
      "identify a temporal difference"
    ),
    paste(
      "BH across all 50 sets; NES does not measure protein or pathway activity"
    )
  )
)

panel_paths <- c(
  metadata = file.path(panel_dir, "fig04_metadata.csv"),
  sample_manifest = file.path(panel_dir, "fig04_sample_manifest.csv"),
  design_counts = file.path(panel_dir, "fig04_design_counts.csv"),
  study_comparison = file.path(panel_dir, "fig04_study_comparison.csv"),
  selected_gene_timepoint = file.path(
    panel_dir,
    "fig04_selected_gene_timepoint.csv"
  ),
  selected_gene_all_estimands = file.path(
    panel_dir,
    "fig04_selected_gene_all_estimands.csv"
  ),
  hallmark_interaction_focus = file.path(
    panel_dir,
    "fig04_hallmark_interaction_focus.csv"
  ),
  hallmark_all_estimands = file.path(
    panel_dir,
    "fig04_hallmark_all_estimands.csv"
  )
)
write_csv(figure_metadata, panel_paths[["metadata"]])
write_csv(sample_manifest_table, panel_paths[["sample_manifest"]])
write_csv(design_counts, panel_paths[["design_counts"]])
write_csv(study_comparison, panel_paths[["study_comparison"]])
write_csv(
  selected_gene_timepoint,
  panel_paths[["selected_gene_timepoint"]]
)
write_csv(
  selected_gene_all_estimands,
  panel_paths[["selected_gene_all_estimands"]]
)
write_csv(
  hallmark_interaction_focus,
  panel_paths[["hallmark_interaction_focus"]]
)
write_csv(
  hallmark_all_estimands,
  panel_paths[["hallmark_all_estimands"]]
)

table_data <- list(
  A_gse283401_sample_manifest = sample_manifest_table,
  B_gene_treatment_contrast_6h = table_s5_b,
  C_gene_treatment_contrast_48h = table_s5_c,
  D_gene_time_by_treatment_interaction = table_s5_d,
  E_hallmark_time_by_treatment_interaction = table_s5_e,
  data_dictionary = table_s5_dictionary
)
table_stems <- c(
  A_gse283401_sample_manifest = "TableS5A_gse283401_sample_manifest",
  B_gene_treatment_contrast_6h = "TableS5B_gene_treatment_contrast_6h",
  C_gene_treatment_contrast_48h = "TableS5C_gene_treatment_contrast_48h",
  D_gene_time_by_treatment_interaction =
    "TableS5D_gene_time_by_treatment_interaction",
  E_hallmark_time_by_treatment_interaction =
    "TableS5E_hallmark_time_by_treatment_interaction",
  data_dictionary = "TableS5_data_dictionary"
)
table_paths <- character()
for (key in names(table_data)) {
  csv_path <- file.path(table_dir, paste0(table_stems[[key]], ".csv"))
  tsv_path <- file.path(table_dir, paste0(table_stems[[key]], ".tsv"))
  write_both(table_data[[key]], csv_path, tsv_path)
  table_paths <- c(
    table_paths,
    setNames(csv_path, paste0(key, "_csv")),
    setNames(tsv_path, paste0(key, "_tsv"))
  )
}

# ---- Manifests, execution record, and session information -----------------

source_manifest <- tibble(
  source_key = names(source_paths),
  source_file = basename(source_paths),
  absolute_path = unname(source_paths),
  sha256 = unname(observed_source_hashes),
  source_type = c(
    "GEO SOFT metadata",
    "Raw gene-count matrix; hippocampus 6 h",
    "Raw gene-count matrix; hippocampus 48 h"
  )
)
write_csv(
  source_manifest,
  file.path(manifest_dir, "fig04_source_manifest.csv")
)

all_data_paths <- c(panel_paths, table_paths)
data_output_manifest <- tibble(
  output_key = names(all_data_paths),
  relative_output = sub(
    paste0("^", figure_root, "/"),
    "",
    unname(all_data_paths)
  ),
  sha256 = vapply(
    unname(all_data_paths),
    sha256_file,
    FUN.VALUE = character(1)
  )
)
write_csv(
  data_output_manifest,
  file.path(manifest_dir, "fig04_data_output_manifest.csv")
)

execution_record <- list(
  figure = "Figure 4",
  supplementary_table = "Table S5",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  external_accession = "GSE283401",
  biological_unit = "animal/library",
  n_hippocampal_source_libraries = nrow(hippocampus_meta),
  n_old_model_libraries = nrow(old_meta),
  old_design = list(
    six_h = list(control = 8L, combined_exposure = 7L),
    forty_eight_h = list(control = 8L, combined_exposure = 9L)
  ),
  source_treatment_label = "surgery",
  manuscript_treatment_label = "1.2% isoflurane plus laparotomy",
  feature_filters = list(
    six_h = paste(
      ">=10 total counts across the 15 old-animal 6 h libraries;",
      "26,317 features"
    ),
    forty_eight_h = paste(
      ">=10 total counts across the 17 old-animal 48 h libraries;",
      "23,732 features"
    ),
    interaction = paste(
      "Shared Ensembl rows with >=10 total counts across all 32 old-animal",
      "libraries; 30,280 features"
    )
  ),
  separate_models = list(
    six_h = "~ treatment",
    forty_eight_h = "~ treatment"
  ),
  interaction_model = "~ time + treatment + time:treatment",
  interaction_coefficient = interaction_name,
  interaction_orientation = paste(
    "(48 h exposure - 48 h control) -",
    "(6 h exposure - 6 h control)"
  ),
  hallmark = list(
    database = "mouse MSigDB 2026.1.Mm",
    msigdbr_call = paste(
      'msigdbr(db_species="MM", species="Mus musculus",',
      'collection="MH")'
    ),
    fallback = FALSE,
    gene_sets = 50L,
    ranking = "signed DESeq2 Wald statistic",
    duplicate_symbol_rule = "largest absolute Wald statistic within estimand",
    minGSSize = 10L,
    maxGSSize = 500L,
    pvalueCutoff = 1,
    eps = 1e-30,
    seed = 42L,
    bh_family = "all 50 Hallmark sets within each estimand"
  ),
  claim_ceiling = paste(
    "Between-experiment difference in treatment-associated contrast",
    "compatible with, but not identifying, a temporal difference"
  ),
  immutable_limitation = paste(
    "Different animals at 6 h and 48 h; time is inseparable from deposited",
    "matrix and any unrecorded time-aligned processing or sequencing-batch",
    "effects; not a within-animal trajectory, replication, or validation"
  ),
  verified_interaction_hallmark_anchors = lapply(
    seq_len(nrow(hallmark_interaction_focus)),
    function(i) {
      list(
        id = hallmark_interaction_focus$ID[[i]],
        nes = hallmark_interaction_focus$NES[[i]],
        fdr = hallmark_interaction_focus$p.adjust[[i]]
      )
    }
  ),
  r_version = R.version.string,
  packages = as.list(setNames(
    vapply(
      c(
        "AnnotationDbi", "clusterProfiler", "DESeq2", "digest", "dplyr",
        "jsonlite", "msigdbr", "org.Mm.eg.db", "readr", "tibble", "tidyr"
      ),
      function(pkg) as.character(packageVersion(pkg)),
      FUN.VALUE = character(1)
    ),
    c(
      "AnnotationDbi", "clusterProfiler", "DESeq2", "digest", "dplyr",
      "jsonlite", "msigdbr", "org.Mm.eg.db", "readr", "tibble", "tidyr"
    )
  ))
)
write_json(
  execution_record,
  file.path(manifest_dir, "fig04_data_preparation_execution.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 16
)

session_path <- file.path(
  manifest_dir,
  "fig04_data_preparation_sessionInfo.txt"
)
session_connection <- file(session_path, open = "wt")
sink(session_connection)
print(sessionInfo())
sink()
close(session_connection)

cat("Figure 4 panel-ready data written to:", panel_dir, "\n")
cat("Table S5 machine-readable files written to:", table_dir, "\n")
cat(
  sprintf(
    "Old-animal design: 6 h 8 control/7 exposed; 48 h 8 control/9 exposed.\n"
  )
)
cat(
  sprintf(
    paste0(
      "Feature universes: 6 h %d; 48 h %d; interaction %d ",
      "Ensembl features.\n"
    ),
    nrow(gene_6h),
    nrow(gene_48h),
    length(interaction_features)
  )
)
cat("Verified interaction Hallmark anchors:\n")
print(
  as.data.frame(
    hallmark_interaction_focus |>
      select(ID, setSize, NES, pvalue, p.adjust)
  ),
  row.names = FALSE,
  digits = 5
)
