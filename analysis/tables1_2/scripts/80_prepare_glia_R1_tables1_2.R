#!/usr/bin/env Rscript

# GLIA major revision — canonical preparation for revised Tables 1 and 2
#
# Table 1: per-animal Scrublet exclusion and UMI audit of the frozen
#          GSE267933 microglial cell set.
# Table 2: per-animal composition of the six non-Rare submitted partitions
#          and animal-level exposure contrasts.
#
# This script uses the analysis-complete Figure 1 panel-ready files as its
# frozen inputs. It independently recomputes every percentage and every
# Table 2 animal-level statistic before creating numeric and display files.
#
# Example:
# Rscript scripts/80_prepare_glia_R1_tables1_2.R \
#   --figure-root /path/to/results_rebuild/figure1 \
#   --table-root /path/to/results_rebuild/tables1_2

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

fmt_int <- function(x) {
  ifelse(is.na(x), "\u2014", formatC(x, format = "f", digits = 0, big.mark = ","))
}

fmt_num <- function(x, digits = 1L) {
  ifelse(
    is.na(x),
    "\u2014",
    formatC(x, format = "f", digits = digits, drop0trailing = FALSE)
  )
}

fmt_num_comma <- function(x, digits = 1L) {
  ifelse(
    is.na(x),
    "\u2014",
    formatC(
      x,
      format = "f",
      digits = digits,
      big.mark = ",",
      drop0trailing = FALSE
    )
  )
}

cohens_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  vx <- var(x)
  vy <- var(y)
  pooled_sd <- sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
  if (pooled_sd == 0) {
    return(NA_real_)
  }
  (mean(x) - mean(y)) / pooled_sd
}

welch_ci <- function(x, y, alpha = 0.05) {
  nx <- length(x)
  ny <- length(y)
  vx <- var(x)
  vy <- var(y)
  difference <- mean(x) - mean(y)
  se <- sqrt(vx / nx + vy / ny)
  if (se == 0) {
    return(c(difference = difference, low = NA_real_, high = NA_real_))
  }
  df <- (vx / nx + vy / ny)^2 /
    ((vx / nx)^2 / (nx - 1) + (vy / ny)^2 / (ny - 1))
  critical <- qt(1 - alpha / 2, df = df)
  c(
    difference = difference,
    low = difference - critical * se,
    high = difference + critical * se
  )
}

exact_permutation_p <- function(x, y) {
  all_values <- c(x, y)
  allocations <- combn(seq_along(all_values), length(x))
  observed <- abs(mean(x) - mean(y))
  permuted <- apply(allocations, 2L, function(exposed_idx) {
    abs(mean(all_values[exposed_idx]) - mean(all_values[-exposed_idx]))
  })
  c(
    p = mean(permuted >= observed - 1e-12),
    n_permutations = length(permuted)
  )
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
default_table_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

figure_root <- normalizePath(
  get_arg("--figure-root", required = TRUE),
  mustWork = TRUE
)
table_root <- normalizePath(
  get_arg("--table-root", default_table_root),
  mustWork = TRUE
)

data_dir <- file.path(table_root, "data", "canonical")
output_dir <- file.path(table_root, "outputs")
manifest_dir <- file.path(table_root, "manifests")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

source_paths <- c(
  sample_qc = file.path(
    figure_root, "data", "panel_ready", "fig01_sample_qc.csv"
  ),
  composition = file.path(
    figure_root, "data", "panel_ready", "fig01_composition.csv"
  ),
  composition_effects = file.path(
    figure_root, "data", "panel_ready", "fig01_composition_effects.csv"
  ),
  figure_manifest = file.path(
    figure_root, "manifests", "fig01_panel_ready_manifest.csv"
  )
)
assert_true(all(file.exists(source_paths)), "At least one frozen Figure 1 input is missing.")

# Confirm that the three table inputs still match the Figure 1 output manifest.
figure_manifest <- read_csv(
  source_paths[["figure_manifest"]],
  show_col_types = FALSE
)
for (key in c("sample_qc", "composition", "composition_effects")) {
  expected_hash <- figure_manifest$sha256[figure_manifest$output_key == key]
  assert_true(length(expected_hash) == 1L, paste("Missing manifest hash for", key))
  assert_true(
    identical(sha256_file(source_paths[[key]]), expected_hash),
    paste("Frozen Figure 1 input hash changed:", key)
  )
}

expected_samples <- c("C1", "C2", "C3", "S1", "S2", "S3")
display_groups <- c(
  Control = "Oxygen control",
  Surgery = "Combined exposure"
)

# ---- Table 1 ---------------------------------------------------------------

sample_qc <- read_csv(source_paths[["sample_qc"]], show_col_types = FALSE)
assert_true(
  identical(sample_qc$sample, expected_samples),
  "Table 1 samples are not in the expected C1-C3/S1-S3 order."
)
assert_true(
  identical(sample_qc$group, c(rep("Control", 3L), rep("Surgery", 3L))),
  "Table 1 source-group labels are inconsistent with the frozen design."
)

sample_qc <- sample_qc |>
  mutate(
    exposure_group = unname(display_groups[group]),
    percent_cells_removed_recomputed =
      100 * predicted_doublets_removed / cells_before,
    percent_umi_removed_recomputed = 100 * umi_removed / umi_before
  )

assert_true(
  all(sample_qc$cells_before - sample_qc$predicted_doublets_removed ==
    sample_qc$cells_retained),
  "Table 1 cell arithmetic failed."
)
assert_true(
  all(sample_qc$umi_before - sample_qc$umi_removed == sample_qc$umi_retained),
  "Table 1 UMI arithmetic failed."
)
assert_true(sum(sample_qc$cells_before) == 7461L, "Expected 7,461 cells before exclusion.")
assert_true(
  sum(sample_qc$predicted_doublets_removed) == 90L,
  "Expected 90 Scrublet-predicted doublets."
)
assert_true(sum(sample_qc$cells_retained) == 7371L, "Expected 7,371 retained cells.")
assert_true(
  sum(sample_qc$umi_before) == 58430610,
  "Unexpected integer-safe UMI total before exclusion."
)
assert_true(
  sum(sample_qc$umi_removed) == 1168668,
  "Unexpected integer-safe removed-UMI total."
)
assert_true(
  sum(sample_qc$umi_retained) == 57261942,
  "Unexpected integer-safe retained-UMI total."
)

table01_numeric <- sample_qc |>
  transmute(
    animal_library = sample,
    exposure_group,
    cells_before,
    predicted_doublets_removed,
    percent_cells_removed = percent_cells_removed_recomputed,
    cells_retained,
    microglial_umi_before = umi_before,
    microglial_umi_removed = umi_removed,
    percent_microglial_umi_removed = percent_umi_removed_recomputed,
    microglial_umi_retained = umi_retained,
    median_umi_per_retained_cell = median_umi_retained_cell
  )

table01_total <- table01_numeric |>
  summarise(
    animal_library = "Total",
    exposure_group = "All animals",
    cells_before = sum(cells_before),
    predicted_doublets_removed = sum(predicted_doublets_removed),
    percent_cells_removed =
      100 * sum(predicted_doublets_removed) / sum(cells_before),
    cells_retained = sum(cells_retained),
    microglial_umi_before = sum(microglial_umi_before),
    microglial_umi_removed = sum(microglial_umi_removed),
    percent_microglial_umi_removed =
      100 * sum(microglial_umi_removed) / sum(microglial_umi_before),
    microglial_umi_retained = sum(microglial_umi_retained),
    median_umi_per_retained_cell = NA_real_
  )
table01_numeric_with_total <- bind_rows(table01_numeric, table01_total)

table01_display <- table01_numeric_with_total |>
  transmute(
    `Animal/library` = animal_library,
    `Exposure group` = exposure_group,
    `Cells before exclusion, n` = fmt_int(cells_before),
    `Predicted doublets removed, n (%)` = paste0(
      fmt_int(predicted_doublets_removed), " (",
      fmt_num(percent_cells_removed, 2L), ")"
    ),
    `Cells retained, n` = fmt_int(cells_retained),
    `Microglial UMIs before exclusion, n` = fmt_int(microglial_umi_before),
    `Microglial UMIs removed, n (%)` = paste0(
      fmt_int(microglial_umi_removed), " (",
      fmt_num(percent_microglial_umi_removed, 2L), ")"
    ),
    `Microglial UMIs retained, n` = fmt_int(microglial_umi_retained),
    `Median UMI per retained cell` = fmt_num_comma(
      median_umi_per_retained_cell, 1L
    )
  )

# ---- Table 2 ---------------------------------------------------------------

composition <- read_csv(source_paths[["composition"]], show_col_types = FALSE)
composition_effects_source <- read_csv(
  source_paths[["composition_effects"]],
  show_col_types = FALSE
)

assert_true(nrow(composition) == 36L, "Expected six partitions x six animals.")
assert_true(
  setequal(composition$partition, 0:5),
  "Table 2 must contain exactly submitted partitions 0-5."
)
assert_true(
  identical(sort(unique(composition$sample)), expected_samples),
  "Table 2 animal identifiers do not match the frozen design."
)
assert_true(sum(composition$n_cells) == 7367L, "Expected 7,367 non-Rare cells.")

denominators <- composition |>
  summarise(non_rare_total = sum(n_cells), .by = sample) |>
  arrange(match(sample, expected_samples))
expected_denominators <- c(1134L, 915L, 1271L, 1301L, 1599L, 1147L)
assert_true(
  identical(denominators$sample, expected_samples) &&
    identical(as.integer(denominators$non_rare_total), expected_denominators),
  "Table 2 non-Rare denominators differ from the independently audited values."
)

composition <- composition |>
  left_join(denominators, by = "sample") |>
  mutate(
    pct_recomputed = 100 * n_cells / non_rare_total,
    sample = factor(sample, levels = expected_samples)
  ) |>
  arrange(partition, sample)
assert_true(
  max(abs(composition$pct_recomputed - composition$pct_of_animal)) < 1e-10,
  "Table 2 source percentages do not equal n/non-Rare denominator."
)

table02_statistics <- bind_rows(lapply(0:5, function(partition_id) {
  block <- composition |>
    filter(partition == partition_id) |>
    arrange(sample)
  exposure <- block$pct_recomputed[block$source_group == "Surgery"]
  control <- block$pct_recomputed[block$source_group == "Control"]
  ci <- welch_ci(exposure, control)
  permutation <- exact_permutation_p(exposure, control)
  tibble(
    partition = partition_id,
    submitted_trace_label = unique(block$submitted_trace_label),
    oxygen_control_mean_percent = mean(control),
    combined_exposure_mean_percent = mean(exposure),
    difference_percentage_points = unname(ci[["difference"]]),
    ci95_low = unname(ci[["low"]]),
    ci95_high = unname(ci[["high"]]),
    cohens_d = cohens_d(exposure, control),
    exact_permutation_p = unname(permutation[["p"]]),
    n_permutations = as.integer(permutation[["n_permutations"]]),
    perfect_separation =
      min(exposure) > max(control) || max(exposure) < min(control)
  )
}))

# Cross-check all recomputed fields against the frozen Figure 1 effect file.
effect_check <- table02_statistics |>
  inner_join(
    composition_effects_source,
    by = c("partition", "submitted_trace_label")
  )
assert_true(nrow(effect_check) == 6L, "Table 2 effect cross-check did not match six rows.")
assert_true(
  all(round(effect_check$oxygen_control_mean_percent, 3) ==
    effect_check$control_mean_pct),
  "Control means disagree with the frozen Figure 1 effects."
)
assert_true(
  all(round(effect_check$combined_exposure_mean_percent, 3) ==
    effect_check$combined_exposure_mean_pct),
  "Combined-exposure means disagree with the frozen Figure 1 effects."
)
assert_true(
  all(round(effect_check$difference_percentage_points.x, 3) ==
    effect_check$difference_percentage_points.y),
  "Mean differences disagree with the frozen Figure 1 effects."
)
assert_true(
  all(round(effect_check$ci95_low.x, 3) == effect_check$ci95_low.y) &&
    all(round(effect_check$ci95_high.x, 3) == effect_check$ci95_high.y),
  "Welch confidence intervals disagree with the frozen Figure 1 effects."
)
assert_true(
  all(round(effect_check$cohens_d.x, 3) == effect_check$cohens_d.y),
  "Cohen's d values disagree with the frozen Figure 1 effects."
)
assert_true(
  all(effect_check$exact_permutation_p.x ==
    effect_check$exact_permutation_p.y),
  "Exact permutation p values disagree with the frozen Figure 1 effects."
)
assert_true(
  all(table02_statistics$n_permutations == 20L),
  "Expected all 20 three-versus-three label allocations."
)
assert_true(
  !any(table02_statistics$perfect_separation),
  "Unexpected perfect animal-level separation."
)

table02_composition_long <- composition |>
  transmute(
    partition,
    submitted_trace_label,
    animal_library = as.character(sample),
    exposure_group = display_group,
    n_cells,
    non_rare_denominator = non_rare_total,
    percent_of_animal = pct_recomputed
  )

table02_part_a <- table02_composition_long |>
  mutate(
    partition_label = paste0(partition, " (", submitted_trace_label, ")"),
    value = paste0(fmt_int(n_cells), " (", fmt_num(percent_of_animal, 1L), ")")
  ) |>
  select(partition, partition_label, animal_library, value) |>
  pivot_wider(names_from = animal_library, values_from = value) |>
  arrange(partition) |>
  select(-partition) |>
  rename(`Submitted partition (traceability label)` = partition_label)

table02_total_row <- tibble(
  `Submitted partition (traceability label)` = "Non-Rare total",
  C1 = paste0(fmt_int(expected_denominators[[1L]]), " (100.0)"),
  C2 = paste0(fmt_int(expected_denominators[[2L]]), " (100.0)"),
  C3 = paste0(fmt_int(expected_denominators[[3L]]), " (100.0)"),
  S1 = paste0(fmt_int(expected_denominators[[4L]]), " (100.0)"),
  S2 = paste0(fmt_int(expected_denominators[[5L]]), " (100.0)"),
  S3 = paste0(fmt_int(expected_denominators[[6L]]), " (100.0)")
)
table02_part_a <- bind_rows(table02_part_a, table02_total_row)

table02_part_b <- table02_statistics |>
  transmute(
    `Submitted partition (traceability label)` =
      paste0(partition, " (", submitted_trace_label, ")"),
    `Oxygen control, mean %` = fmt_num(oxygen_control_mean_percent, 1L),
    `Combined exposure, mean %` = fmt_num(
      combined_exposure_mean_percent, 1L
    ),
    `Difference (95% CI), percentage points` = paste0(
      fmt_num(difference_percentage_points, 1L), " (",
      fmt_num(ci95_low, 1L), " to ", fmt_num(ci95_high, 1L), ")"
    ),
    `Cohen's d` = fmt_num(cohens_d, 2L),
    `Exact P` = fmt_num(exact_permutation_p, 2L)
  )

# ---- Titles and notes ------------------------------------------------------

table_text <- list(
  table1 = list(
    title = paste(
      "Table 1. Per-animal Scrublet exclusion and UMI audit of the frozen",
      "GSE267933 microglial cell set"
    ),
    note = paste(
      "Values refer to the frozen submitted microglial cell set, not to all",
      "barcodes in each sequencing library. Scrublet was applied separately",
      "to each animal/library; predicted doublets are model-based calls and",
      "are not confirmed doublets. Percentages removed were calculated",
      "relative to the corresponding before-exclusion cell or UMI total.",
      "The retained cells constitute the 7,371-cell label-independent",
      "molecular-analysis set. C1-C3 were oxygen controls, whereas S1-S3",
      "received the combined exposure of 2.5% sevoflurane in 50% O2 for",
      "30 min plus laparotomy. No anesthesia-only group was available.",
      "UMI, unique molecular identifier."
    )
  ),
  table2 = list(
    title = paste(
      "Table 2. Per-animal composition of the six non-Rare submitted",
      "partitions and conditional exposure contrasts (GSE267933)"
    ),
    part_a = "Part A. Per-animal composition",
    part_b = "Part B. Animal-level group contrasts",
    note = paste(
      "Entries in Part A are cell count (percentage). Percentages were",
      "calculated within each animal after excluding the submitted",
      "partition 6/Rare label (94 cells) from the frozen 7,461-cell set,",
      "yielding 7,367 cells; Scrublet status was not used to define this",
      "composition set. This differs intentionally from the 7,371-cell",
      "molecular-analysis set in Table 1. Biological names shown in",
      "parentheses are retained only for traceability to the submitted",
      "manuscript and are not treated as validated biological states.",
      "Part B uses the animal/library as the biological unit. Differences",
      "are the unweighted combined-exposure-minus-oxygen-control mean",
      "differences. Confidence intervals are Welch-Satterthwaite 95%",
      "intervals; Cohen's d is a point estimate using the pooled",
      "within-group standard deviation. Exact two-sided P values were",
      "obtained from all 20 allocations of three of the six animals to the",
      "exposed group. These unadjusted tests are descriptive; the six",
      "partition percentages are compositional and not independent.",
      "CI, confidence interval."
    )
  )
)

# ---- Deterministic outputs and manifests ----------------------------------

output_paths <- c(
  table01_numeric = file.path(data_dir, "table01_qc_numeric.csv"),
  table02_composition_long = file.path(
    data_dir, "table02_composition_long.csv"
  ),
  table02_statistics = file.path(data_dir, "table02_statistics.csv"),
  table01_display_csv = file.path(
    output_dir, "Table1_microglial_analysis_set_audit.csv"
  ),
  table01_display_tsv = file.path(
    output_dir, "Table1_microglial_analysis_set_audit.tsv"
  ),
  table02_part_a_csv = file.path(
    output_dir, "Table2_partA_per_animal_composition.csv"
  ),
  table02_part_a_tsv = file.path(
    output_dir, "Table2_partA_per_animal_composition.tsv"
  ),
  table02_part_b_csv = file.path(
    output_dir, "Table2_partB_animal_level_contrasts.csv"
  ),
  table02_part_b_tsv = file.path(
    output_dir, "Table2_partB_animal_level_contrasts.tsv"
  ),
  table_text = file.path(data_dir, "tables1_2_text.json")
)

write_csv(table01_numeric_with_total, output_paths[["table01_numeric"]])
write_csv(table02_composition_long, output_paths[["table02_composition_long"]])
write_csv(table02_statistics, output_paths[["table02_statistics"]])
write_csv(table01_display, output_paths[["table01_display_csv"]])
write_tsv(table01_display, output_paths[["table01_display_tsv"]])
write_csv(table02_part_a, output_paths[["table02_part_a_csv"]])
write_tsv(table02_part_a, output_paths[["table02_part_a_tsv"]])
write_csv(table02_part_b, output_paths[["table02_part_b_csv"]])
write_tsv(table02_part_b, output_paths[["table02_part_b_tsv"]])
write_json(
  table_text,
  output_paths[["table_text"]],
  pretty = TRUE,
  auto_unbox = TRUE
)

source_manifest <- tibble(
  source_key = names(source_paths),
  relative_source = c(
    "figure1/data/panel_ready/fig01_sample_qc.csv",
    "figure1/data/panel_ready/fig01_composition.csv",
    "figure1/data/panel_ready/fig01_composition_effects.csv",
    "figure1/manifests/fig01_panel_ready_manifest.csv"
  ),
  role = c(
    "Table 1 frozen per-animal cell and UMI audit",
    "Table 2 frozen per-animal non-Rare counts",
    "Table 2 Figure 1 cross-check values",
    "Frozen-input hash verification"
  ),
  sha256 = vapply(source_paths, sha256_file, character(1))
)
write_csv(
  source_manifest,
  file.path(manifest_dir, "tables1_2_source_manifest.csv")
)

output_manifest <- tibble(
  output_key = names(output_paths),
  relative_output = vapply(
    output_paths,
    function(path) sub(
      paste0("^", normalizePath(table_root), "/?"),
      "",
      normalizePath(path)
    ),
    character(1)
  ),
  sha256 = vapply(output_paths, sha256_file, character(1))
)
write_csv(
  output_manifest,
  file.path(manifest_dir, "tables1_2_output_manifest.csv")
)

execution_manifest <- list(
  tables = c("Table 1", "Table 2"),
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  statistical_unit = "animal/library",
  cell_sets = list(
    table1_before_exclusion = 7461L,
    table1_scrublet_negative = 7371L,
    table2_non_rare_conditional = 7367L
  ),
  table2_statistics = list(
    contrast = "combined exposure minus oxygen control",
    confidence_interval = "Welch-Satterthwaite 95% mean-difference interval",
    effect_size = "Cohen's d using pooled within-group SD",
    exact_test = paste(
      "two-sided absolute mean-difference permutation over all",
      "20 three-versus-three allocations; no +1 correction"
    )
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
  file.path(manifest_dir, "tables1_2_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "tables1_2_sessionInfo.txt")
)

message("Canonical Table 1 and Table 2 files written to: ", output_dir)
message("Validated cell sets: Table 1, 7,461 -> 7,371; Table 2, 7,367.")
