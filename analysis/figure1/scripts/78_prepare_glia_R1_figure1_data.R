#!/usr/bin/env Rscript

# GLIA major revision — Figure 1 panel-ready data preparation
#
# This script consolidates and validates the frozen source files used by
# Figure 1. It does not recompute clustering or alter submitted labels.
#
# Example:
# Rscript scripts/78_prepare_glia_R1_figure1_data.R \
#   --project-root /path/to/pocd_scrna \
#   --audit-root /path/to/independent_audit \
#   --figure-root /path/to/figure1

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

script_path <- get_script_path()
script_dir <- dirname(script_path)
default_figure_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

project_root <- normalizePath(
  get_arg("--project-root", required = TRUE),
  mustWork = TRUE
)
audit_root <- normalizePath(
  get_arg("--audit-root", required = TRUE),
  mustWork = TRUE
)
figure_root <- normalizePath(
  get_arg("--figure-root", default_figure_root),
  mustWork = TRUE
)

panel_dir <- file.path(figure_root, "data", "panel_ready")
manifest_dir <- file.path(figure_root, "manifests")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

project_source <- function(relative_path) {
  path <- file.path(project_root, relative_path)
  assert_true(file.exists(path), paste("Missing project source:", relative_path))
  path
}

audit_source <- function(filename) {
  path <- file.path(audit_root, filename)
  assert_true(file.exists(path), paste("Missing audit source:", filename))
  path
}

source_paths <- c(
  umap = project_source("data/results/r_export_v3/fig1_umap.csv"),
  umap_provenance = project_source("provenance/export_figures_v3_provenance.json"),
  composition = project_source("data/results/revision_R1/B_per_animal_composition.csv"),
  composition_effects = project_source("data/results/revision_R1/B_donor_level_tests.csv"),
  composition_provenance = project_source("provenance/B_per_animal_composition_provenance.json"),
  count_audit = audit_source("scrublet_negative_count_audit.csv"),
  seed_stability = audit_source("subtyped_stored_graph_seed_stability.csv")
)

source_roles <- c(
  umap = "Panel B: frozen UMAP coordinates and submitted traceability labels",
  umap_provenance = "Panel B: export provenance",
  composition = "Panel D: per-animal non-Rare conditional composition",
  composition_effects = "Panel D: animal-level group differences and uncertainty",
  composition_provenance = "Panel D: analysis provenance",
  count_audit = "Panel A: submitted and Scrublet-negative analysis-set counts",
  seed_stability = "Panel C: corrected stored-graph seed sensitivity"
)

expected_samples <- c("C1", "C2", "C3", "S1", "S2", "S3")
expected_source_groups <- c(
  C1 = "Control", C2 = "Control", C3 = "Control",
  S1 = "Surgery", S2 = "Surgery", S3 = "Surgery"
)
display_groups <- c(
  Control = "Oxygen control",
  Surgery = "Combined exposure"
)
trace_to_partition <- c(
  "Inflammatory" = 0L,
  "Transitional-A" = 1L,
  "Homeostatic-A" = 2L,
  "Homeostatic-B" = 3L,
  "Homeostatic-C" = 4L,
  "Transitional-B" = 5L,
  "Rare" = 6L
)
expected_partition_counts <- c(
  "0" = 1650L,
  "1" = 1381L,
  "2" = 1285L,
  "3" = 1224L,
  "4" = 1209L,
  "5" = 618L,
  "6" = 94L
)

# ---- Panel A: study design and analysis-set counts -------------------------

count_audit <- read_csv(source_paths[["count_audit"]], show_col_types = FALSE)
assert_true(
  identical(count_audit$sample, expected_samples),
  "Count-audit samples are not C1-C3/S1-S3 in the expected order."
)
assert_true(
  identical(count_audit$group, unname(expected_source_groups[expected_samples])),
  "Count-audit source group labels do not match the frozen design."
)
assert_true(sum(count_audit$cells_before) == 7461L, "Expected 7,461 submitted microglia.")
assert_true(
  sum(count_audit$predicted_doublets_removed) == 90L,
  "Expected 90 Scrublet-predicted doublets."
)
assert_true(sum(count_audit$cells_retained) == 7371L, "Expected 7,371 retained cells.")

design <- tibble::tribble(
  ~source_group, ~display_group, ~samples, ~n_animals, ~exposure,
  "Control", "Oxygen control", "C1-C3", 3L,
  "50% oxygen for 30 min; no anesthesia or laparotomy",
  "Surgery", "Combined exposure", "S1-S3", 3L,
  "2.5% sevoflurane in 50% oxygen for 30 min plus exploratory laparotomy"
) |>
  mutate(
    dataset = "GSE267933",
    age = "18 months",
    sex = "male",
    strain = "C57BL/6",
    tissue = "hippocampus",
    sampling = "24 h after exposure",
    assay = "10x Genomics Chromium Single Cell 3' v2"
  )

analysis_sets <- tibble::tribble(
  ~analysis_set, ~n_cells, ~rule, ~figure_panels,
  "Submitted-partition audit", 7461L,
  "All frozen submitted microglia", "1B-1C",
  "Conditional composition", 7367L,
  "Partition 6 (Rare; 94 cells) excluded; Scrublet status not used", "1D",
  "Primary molecular analyses", 7371L,
  "90 Scrublet-predicted doublets excluded without conditioning on submitted partition", "Figures 2-4"
)

# ---- Panel B: UMAP and submitted numeric partitions -----------------------

umap <- read_csv(source_paths[["umap"]], show_col_types = FALSE, name_repair = "minimal")
assert_true(ncol(umap) == 6L, "UMAP source should contain six columns.")
names(umap)[1L] <- "barcode"
required_umap <- c("barcode", "UMAP1", "UMAP2", "subtype", "group", "sample")
assert_true(identical(names(umap), required_umap), "Unexpected UMAP source columns.")
assert_true(nrow(umap) == 7461L, "UMAP source must contain 7,461 cells.")
assert_true(!anyDuplicated(umap$barcode), "UMAP cell barcodes are not unique.")
assert_true(!anyNA(umap), "UMAP source contains missing values.")
assert_true(
  setequal(unique(umap$subtype), names(trace_to_partition)),
  "UMAP traceability labels do not match the frozen seven-label mapping."
)

umap <- umap |>
  mutate(
    partition = unname(trace_to_partition[subtype]),
    source_group = group,
    display_group = unname(display_groups[source_group])
  ) |>
  select(
    barcode, UMAP1, UMAP2, partition,
    submitted_trace_label = subtype,
    sample, source_group, display_group
  ) |>
  arrange(barcode)

observed_partition_counts <- table(umap$partition)
assert_true(
  identical(
    as.integer(observed_partition_counts[names(expected_partition_counts)]),
    as.integer(expected_partition_counts)
  ),
  "Submitted partition counts do not match the frozen object."
)
observed_sample_groups <- umap |>
  distinct(sample, source_group) |>
  arrange(match(sample, expected_samples))
assert_true(
  identical(observed_sample_groups$sample, expected_samples) &&
    identical(
      observed_sample_groups$source_group,
      unname(expected_source_groups[expected_samples])
    ),
  "UMAP sample-to-group mapping is inconsistent with the design."
)

# ---- Panel C: stored-graph seed sensitivity -------------------------------

seed_stability <- read_csv(source_paths[["seed_stability"]], show_col_types = FALSE)
expected_seed_stability <- tibble::tribble(
  ~seed, ~n_clusters, ~ari_vs_original,
  0L, 7L, 1.0000000000000000,
  1L, 6L, 0.49092414553881564,
  7L, 7L, 0.5858529047117915,
  42L, 7L, 0.5127939875974998,
  123L, 6L, 0.4033111527516049,
  2024L, 7L, 0.5165138553705313
)
assert_true(
  identical(as.integer(seed_stability$seed), expected_seed_stability$seed) &&
    identical(as.integer(seed_stability$n_clusters), expected_seed_stability$n_clusters) &&
    max(abs(seed_stability$ari_vs_original - expected_seed_stability$ari_vs_original)) < 1e-12,
  "Stored-graph seed sensitivity differs from the independently verified values."
)
seed_stability <- seed_stability |>
  mutate(
    diagnostic = "Fixed stored graph; Leiden resolution 0.4",
    submitted_seed = seed == 0L
  )

# ---- Panel D: conditional per-animal composition --------------------------

composition <- read_csv(source_paths[["composition"]], show_col_types = FALSE)
required_composition <- c("sample", "group", "mg_subtype2", "n_cells", "pct_of_animal")
assert_true(
  identical(names(composition), required_composition),
  "Unexpected per-animal composition columns."
)
assert_true(nrow(composition) == 36L, "Composition must have six labels x six animals.")
assert_true(sum(composition$n_cells) == 7367L, "Conditional composition must contain 7,367 cells.")
assert_true(
  setequal(unique(composition$mg_subtype2), names(trace_to_partition)[1:6]),
  "Composition must contain exactly the six non-Rare submitted labels."
)
composition_sums <- composition |>
  summarise(
    n = sum(n_cells),
    pct = sum(pct_of_animal),
    .by = sample
  ) |>
  arrange(match(sample, expected_samples))
assert_true(
  identical(composition_sums$sample, expected_samples) &&
    max(abs(composition_sums$pct - 100)) < 1e-8,
  "Conditional composition percentages do not sum to 100% within every animal."
)

composition <- composition |>
  mutate(
    partition = unname(trace_to_partition[mg_subtype2]),
    source_group = group,
    display_group = unname(display_groups[source_group]),
    sample = factor(sample, levels = expected_samples)
  ) |>
  arrange(sample, partition) |>
  mutate(sample = as.character(sample)) |>
  select(
    sample, source_group, display_group, partition,
    submitted_trace_label = mg_subtype2,
    n_cells, pct_of_animal
  )

composition_effects <- read_csv(
  source_paths[["composition_effects"]],
  show_col_types = FALSE
)
assert_true(nrow(composition_effects) == 6L, "Expected six composition effects.")
assert_true(
  setequal(composition_effects$subtype, names(trace_to_partition)[1:6]),
  "Composition effects do not cover the six non-Rare submitted labels."
)
assert_true(
  abs(
    composition_effects$diff_pct_surgery_minus_control[
      composition_effects$subtype == "Transitional-A"
    ] - 6.130
  ) < 1e-9,
  "Transitional-A composition anchor differs from +6.130 percentage points."
)
assert_true(
  abs(
    composition_effects$diff_pct_surgery_minus_control[
      composition_effects$subtype == "Inflammatory"
    ] - 2.316
  ) < 1e-9,
  "Inflammatory composition anchor differs from +2.316 percentage points."
)

composition_effects <- composition_effects |>
  mutate(partition = unname(trace_to_partition[subtype])) |>
  arrange(partition) |>
  select(
    partition,
    submitted_trace_label = subtype,
    control_mean_pct,
    combined_exposure_mean_pct = surgery_mean_pct,
    difference_percentage_points = diff_pct_surgery_minus_control,
    ci95_low = welch_ci95_low,
    ci95_high = welch_ci95_high,
    cohens_d,
    exact_permutation_p = exact_perm_p_two_sided,
    n_permutations,
    perfect_separation
  )
assert_true(
  !any(composition_effects$perfect_separation),
  "Unexpected perfect animal-level separation in a submitted label."
)

# ---- Write deterministic panel-ready files --------------------------------

output_paths <- c(
  design = file.path(panel_dir, "fig01_design.csv"),
  analysis_sets = file.path(panel_dir, "fig01_analysis_sets.csv"),
  sample_qc = file.path(panel_dir, "fig01_sample_qc.csv"),
  umap = file.path(panel_dir, "fig01_umap.csv"),
  seed_stability = file.path(panel_dir, "fig01_seed_stability.csv"),
  composition = file.path(panel_dir, "fig01_composition.csv"),
  composition_effects = file.path(panel_dir, "fig01_composition_effects.csv")
)

write_csv(design, output_paths[["design"]])
write_csv(analysis_sets, output_paths[["analysis_sets"]])
write_csv(count_audit, output_paths[["sample_qc"]])
write_csv(umap, output_paths[["umap"]])
write_csv(seed_stability, output_paths[["seed_stability"]])
write_csv(composition, output_paths[["composition"]])
write_csv(composition_effects, output_paths[["composition_effects"]])

source_manifest <- tibble(
  source_key = names(source_paths),
  relative_source = c(
    "data/results/r_export_v3/fig1_umap.csv",
    "provenance/export_figures_v3_provenance.json",
    "data/results/revision_R1/B_per_animal_composition.csv",
    "data/results/revision_R1/B_donor_level_tests.csv",
    "provenance/B_per_animal_composition_provenance.json",
    "independent_audit/scrublet_negative_count_audit.csv",
    "independent_audit/subtyped_stored_graph_seed_stability.csv"
  ),
  role = unname(source_roles[names(source_paths)]),
  sha256 = vapply(source_paths, sha256_file, character(1))
)
output_manifest <- tibble(
  output_key = names(output_paths),
  relative_output = file.path("data", "panel_ready", basename(output_paths)),
  sha256 = vapply(output_paths, sha256_file, character(1))
)

write_csv(source_manifest, file.path(manifest_dir, "fig01_source_manifest.csv"))
write_csv(output_manifest, file.path(manifest_dir, "fig01_panel_ready_manifest.csv"))

execution_manifest <- list(
  figure = "Figure 1",
  purpose = "Study design and audit of the submitted microglial partition",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  source_manifest_sha256 = sha256_file(
    file.path(manifest_dir, "fig01_source_manifest.csv")
  ),
  panel_ready_manifest_sha256 = sha256_file(
    file.path(manifest_dir, "fig01_panel_ready_manifest.csv")
  ),
  cell_sets = list(
    submitted_partition_audit = 7461L,
    conditional_non_rare_composition = 7367L,
    scrublet_negative_molecular_analysis = 7371L
  ),
  partition_mapping = as.list(trace_to_partition),
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
  file.path(manifest_dir, "fig01_data_preparation_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "fig01_data_preparation_sessionInfo.txt")
)

message("Figure 1 panel-ready data written to: ", panel_dir)
message("Validated cell sets: 7,461 submitted; 7,367 non-Rare composition; 7,371 Scrublet-negative.")
