#!/usr/bin/env Rscript

# GLIA major revision — Figure S2 and Table S2 canonical data preparation
#
# This script consolidates two deliberately distinct algorithmic diagnostics:
#   1. Leiden seed sensitivity on the exact neighborhood graph stored in the
#      submitted 7,461-cell microglial object; and
#   2. conditional 80% cell-subsampling and resolution diagnostics in which a
#      neighborhood graph is rebuilt from inherited whole-cell PC coordinates.
#
# The latter are not tests on the submitted stored graph, end-to-end pipeline
# stability, bootstrap biological replication, or validation of cell states.

suppressPackageStartupMessages({
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
  if (length(file_arg) != 1L) stop("Unable to resolve the current script path.")
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

assert_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  assert_true(
    length(missing_columns) == 0L,
    paste(label, "is missing:", paste(missing_columns, collapse = ", "))
  )
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
default_figure_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
figure_root <- normalizePath(
  get_arg("--figure-root", default_figure_root),
  mustWork = TRUE
)
default_workspace_root <- normalizePath(
  file.path(figure_root, "..", "..", ".."),
  mustWork = TRUE
)
workspace_root <- normalizePath(
  get_arg("--workspace-root", default_workspace_root),
  mustWork = TRUE
)
project_root_arg <- get_arg("--project-root")
assert_true(!is.null(project_root_arg), "Required argument not supplied: --project-root")
project_root <- normalizePath(
  project_root_arg,
  mustWork = TRUE
)

panel_dir <- file.path(figure_root, "data", "panel_ready")
table_dir <- file.path(figure_root, "outputs", "tableS2")
manifest_dir <- file.path(figure_root, "manifests")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

cluster_stability_dir <- file.path(
  project_root,
  "data", "results", "cluster_stability"
)
provenance_dir <- file.path(project_root, "provenance")

input_paths <- c(
  stored_graph_seed = file.path(
    workspace_root,
    "independent_audit", "subtyped_stored_graph_seed_stability.csv"
  ),
  stored_graph_pairwise = file.path(
    workspace_root,
    "independent_audit", "subtyped_stored_graph_pairwise_ari.csv"
  ),
  stored_graph_audit_script = file.path(
    workspace_root,
    "independent_audit", "audit_cluster_stability_inputs.py"
  ),
  subsampling_partition = file.path(
    cluster_stability_dir,
    "bootstrap_jaccard_by_subtype.csv"
  ),
  subsampling_cluster_counts = file.path(
    cluster_stability_dir,
    "bootstrap_cluster_counts.csv"
  ),
  resolution_summary = file.path(
    cluster_stability_dir,
    "resolution_sweep_stability.csv"
  ),
  resolution_cluster = file.path(
    cluster_stability_dir,
    "resolution_sweep_percluster.csv"
  ),
  subsampling_provenance = file.path(
    provenance_dir,
    "cluster_bootstrap_provenance.json"
  ),
  resolution_provenance = file.path(
    provenance_dir,
    "resolution_sweep_provenance.json"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing Figure S2/Table S2 input(s):", paste(missing_inputs, collapse = ", "))
)

# Files known to be unsuitable as evidence about the submitted stored graph.
# They are recorded for exclusion only and are never read below.
excluded_inputs <- c(
  file.path(cluster_stability_dir, "leiden_stability_isolated.csv"),
  file.path(cluster_stability_dir, "perturbation_ari.csv")
)

stored_seed <- read_csv(input_paths[["stored_graph_seed"]], show_col_types = FALSE)
stored_pair <- read_csv(input_paths[["stored_graph_pairwise"]], show_col_types = FALSE)
subsampling_partition_source <- read_csv(
  input_paths[["subsampling_partition"]],
  show_col_types = FALSE
)
subsampling_counts_source <- read_csv(
  input_paths[["subsampling_cluster_counts"]],
  show_col_types = FALSE
)
resolution_summary_source <- read_csv(
  input_paths[["resolution_summary"]],
  show_col_types = FALSE
)
resolution_cluster_source <- read_csv(
  input_paths[["resolution_cluster"]],
  show_col_types = FALSE
)
subsampling_provenance <- read_json(input_paths[["subsampling_provenance"]])
resolution_provenance <- read_json(input_paths[["resolution_provenance"]])

assert_columns(
  stored_seed,
  c("seed", "n_clusters", "ari_vs_original"),
  "stored-graph seed results"
)
assert_columns(
  stored_pair,
  c("seed_a", "seed_b", "ari"),
  "stored-graph pairwise ARI"
)
assert_columns(
  subsampling_partition_source,
  c(
    "subtype", "n_cells", "jaccard_mean_seed_varies",
    "jaccard_sd_seed_varies", "jaccard_median_seed_varies",
    "recovery_rate_jaccard_gt_0.5", "jaccard_mean_seed_fixed"
  ),
  "subsampling per-partition output"
)
assert_columns(
  subsampling_counts_source,
  c("arm", "n_clusters"),
  "subsampling cluster counts"
)
assert_columns(
  resolution_summary_source,
  c(
    "resolution", "n_clusters_full_cohort", "n_clusters_bootstrap_mode",
    "jaccard_mean_across_clusters", "jaccard_median_across_clusters",
    "jaccard_min", "jaccard_max", "n_clusters_ge_0.75",
    "n_clusters_ge_0.60"
  ),
  "resolution summary"
)
assert_columns(
  resolution_cluster_source,
  c("resolution", "cluster", "n_cells", "jaccard_mean"),
  "resolution per-cluster output"
)

seeds <- c(0L, 1L, 7L, 42L, 123L, 2024L)
expected_seed_ari <- c(
  1.0,
  0.49092414553881564,
  0.5858529047117915,
  0.5127939875974998,
  0.4033111527516049,
  0.5165138553705313
)
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)

assert_true(
  nrow(stored_seed) == 6L &&
    identical(as.integer(stored_seed$seed), seeds) &&
    identical(as.integer(stored_seed$n_clusters), c(7L, 6L, 7L, 7L, 6L, 7L)) &&
    max(abs(stored_seed$ari_vs_original - expected_seed_ari)) < 1e-12,
  "Correct stored-graph seed anchors do not match the independent audit."
)
assert_true(
  nrow(stored_pair) == choose(length(seeds), 2L) &&
    all(stored_pair$seed_a %in% seeds) &&
    all(stored_pair$seed_b %in% seeds) &&
    all(stored_pair$seed_a != stored_pair$seed_b) &&
    nrow(distinct(stored_pair, seed_a, seed_b)) == choose(length(seeds), 2L),
  "Stored-graph pairwise ARI is not the complete six-seed upper triangle."
)
assert_true(
  sha256_file(input_paths[["stored_graph_seed"]]) ==
    "531db343727ecb498134b805cc5b8e776d3917a10d51f176e8edc9444ebdab3f" &&
    sha256_file(input_paths[["stored_graph_pairwise"]]) ==
      "7ed501ef6808cf6f8a8bcc91a4a6b78ecaa5f931e64ea999e63c5cbe348ea8c9",
  "Stored-graph independent-audit files changed after the canonical audit."
)

expected_project_hashes <- c(
  subsampling_partition =
    "a47a196c8103e6303e7b221237033ef1952dfbbe66b5379131bbc3f2b9ccdf2e",
  subsampling_cluster_counts =
    "3df379ae3dd8c19bd50f12f77638c7bb70fe8225b811a9d64ec493c4fbfe6153",
  resolution_summary =
    "9bcf8d5e6de6f4449cfb70226a045ad41523ea52a27e818189f53c9d6a030be1",
  resolution_cluster =
    "c40ab094f0e3f7e37d4f2ecd7304020de0907c46cda19556177da7502dc6c5c3"
)
actual_project_hashes <- vapply(
  input_paths[names(expected_project_hashes)],
  sha256_file,
  character(1)
)
assert_true(
  identical(unname(actual_project_hashes), unname(expected_project_hashes)),
  "One or more conditional diagnostic outputs changed from their provenance records."
)
assert_true(
  subsampling_provenance$output_hashes$bootstrap_jaccard_by_subtype.csv ==
    expected_project_hashes[["subsampling_partition"]] &&
    subsampling_provenance$output_hashes$bootstrap_cluster_counts.csv ==
      expected_project_hashes[["subsampling_cluster_counts"]] &&
    resolution_provenance$output_hashes$resolution_sweep_stability.csv ==
      expected_project_hashes[["resolution_summary"]] &&
    resolution_provenance$output_hashes$resolution_sweep_percluster.csv ==
      expected_project_hashes[["resolution_cluster"]],
  "Project provenance JSON does not match the diagnostic CSV files."
)

assert_true(
  nrow(subsampling_partition_source) == 7L &&
    sum(subsampling_partition_source$n_cells) == 7461L,
  "Subsampling per-partition output does not cover all 7,461 submitted cells."
)
assert_true(
  nrow(subsampling_counts_source) == 150L &&
    sum(subsampling_counts_source$arm == "seed varies") == 100L &&
    sum(subsampling_counts_source$arm == "seed fixed") == 50L,
  "Subsampling run counts must be 100 variable-seed and 50 fixed-seed runs."
)
assert_true(
  nrow(resolution_summary_source) == 8L &&
    max(abs(resolution_summary_source$resolution - resolutions)) < 1e-12,
  "Resolution summary does not contain the frozen eight-resolution sweep."
)
assert_true(
  nrow(resolution_cluster_source) == 61L &&
    all(
      resolution_cluster_source |>
        group_by(resolution) |>
        summarise(total_cells = sum(n_cells), .groups = "drop") |>
        pull(total_cells) == 7461L
    ),
  "Resolution per-cluster rows do not partition all 7,461 cells at each resolution."
)

partition_map <- tribble(
  ~partition, ~submitted_trace_label,
  0L, "Inflammatory",
  1L, "Transitional-A",
  2L, "Homeostatic-A",
  3L, "Homeostatic-B",
  4L, "Homeostatic-C",
  5L, "Transitional-B",
  6L, "Rare"
)
assert_true(
  setequal(
    subsampling_partition_source$subtype,
    partition_map$submitted_trace_label
  ),
  "Submitted trace labels do not map one-to-one to partitions 0-6."
)

# Table S2A: correct submitted stored-graph seed sensitivity.
table_s2a <- stored_seed |>
  transmute(
    graph_scope = "Exact graph stored in adata_microglia_subtyped.h5ad",
    leiden_resolution = 0.4,
    seed = as.integer(seed),
    n_partitions = as.integer(n_clusters),
    adjusted_rand_index_vs_submitted = ari_vs_original,
    reproduces_submitted_partition = abs(ari_vs_original - 1) < 1e-12
  )

# Table S2B: complete unique pairwise comparisons among the same fixed-graph runs.
table_s2b <- stored_pair |>
  transmute(
    graph_scope = "Exact graph stored in adata_microglia_subtyped.h5ad",
    leiden_resolution = 0.4,
    seed_a = as.integer(seed_a),
    seed_b = as.integer(seed_b),
    adjusted_rand_index = ari
  ) |>
  arrange(match(seed_a, seeds), match(seed_b, seeds))

# Table S2C: conditional 80% cell-subsampling results, mapped back to numeric
# submitted partitions. Traceability names are retained only in parentheses.
table_s2c <- subsampling_partition_source |>
  select(-any_of("interpretation")) |>
  left_join(partition_map, by = c("subtype" = "submitted_trace_label")) |>
  transmute(
    graph_scope = paste(
      "Graph rebuilt per 80% cell subsample from inherited whole-cell",
      "PC1-PC20; 15 neighbors; Leiden resolution 0.4"
    ),
    partition = as.integer(partition),
    submitted_trace_label = subtype,
    n_cells_in_submitted_partition = as.integer(n_cells),
    n_variable_seed_runs = 100L,
    mean_best_match_jaccard_variable_seed = jaccard_mean_seed_varies,
    sd_best_match_jaccard_variable_seed = jaccard_sd_seed_varies,
    median_best_match_jaccard_variable_seed = jaccard_median_seed_varies,
    proportion_variable_seed_runs_jaccard_gt_0.5 =
      recovery_rate_jaccard_gt_0.5,
    n_fixed_seed_runs = 50L,
    mean_best_match_jaccard_fixed_seed = jaccard_mean_seed_fixed,
    variable_seed_mean_ge_0.75 = jaccard_mean_seed_varies >= 0.75
  ) |>
  arrange(partition)

# Table S2D: complete cluster-count distribution for both subsampling arms.
subsampling_arm_metadata <- tribble(
  ~source_arm, ~display_arm, ~n_runs, ~neighbor_graph_seed, ~leiden_seed,
  "seed varies", "Graph and Leiden seeds varied", 100L,
  "Iteration-specific (0-99)", "Iteration-specific (0-99)",
  "seed fixed", "Graph and Leiden seeds fixed", 50L,
  "0", "0"
)
cluster_count_range <- seq.int(
  min(subsampling_counts_source$n_clusters),
  max(subsampling_counts_source$n_clusters)
)
table_s2d <- subsampling_counts_source |>
  count(source_arm = arm, n_partitions_recovered = n_clusters, name = "n_runs_at_count") |>
  complete(
    source_arm = subsampling_arm_metadata$source_arm,
    n_partitions_recovered = cluster_count_range,
    fill = list(n_runs_at_count = 0L)
  ) |>
  left_join(subsampling_arm_metadata, by = "source_arm") |>
  transmute(
    graph_scope = paste(
      "Graph rebuilt per 80% cell subsample from inherited whole-cell",
      "PC1-PC20; 15 neighbors; Leiden resolution 0.4"
    ),
    arm = display_arm,
    n_runs,
    neighbor_graph_seed,
    leiden_seed,
    n_partitions_recovered = as.integer(n_partitions_recovered),
    n_runs_at_count = as.integer(n_runs_at_count),
    percentage_of_runs = 100 * n_runs_at_count / n_runs
  ) |>
  arrange(match(arm, subsampling_arm_metadata$display_arm), n_partitions_recovered)

# Table S2E/F: conditional reconstructed-graph resolution sweep.
table_s2e <- resolution_summary_source |>
  transmute(
    graph_scope = paste(
      "Full-cohort reference graph and each 80% subsample graph rebuilt",
      "from inherited whole-cell PC1-PC20; 15 neighbors"
    ),
    n_subsampling_runs = 30L,
    reference_and_iteration_seed_scheme = paste(
      "Reference graph/Leiden seed 0; subsample graph/Leiden seed varies 0-29"
    ),
    leiden_resolution = resolution,
    n_partitions_full_cohort = as.integer(n_clusters_full_cohort),
    modal_n_partitions_across_subsamples = as.integer(n_clusters_bootstrap_mode),
    mean_best_match_jaccard_across_full_cohort_partitions =
      jaccard_mean_across_clusters,
    median_best_match_jaccard_across_full_cohort_partitions =
      jaccard_median_across_clusters,
    minimum_best_match_jaccard = jaccard_min,
    maximum_best_match_jaccard = jaccard_max,
    n_partitions_mean_jaccard_ge_0.75 = as.integer(n_clusters_ge_0.75),
    n_partitions_mean_jaccard_ge_0.60 = as.integer(n_clusters_ge_0.60)
  )

table_s2f <- resolution_cluster_source |>
  transmute(
    graph_scope = paste(
      "Full-cohort reference graph and each 80% subsample graph rebuilt",
      "from inherited whole-cell PC1-PC20; 15 neighbors"
    ),
    n_subsampling_runs = 30L,
    leiden_resolution = resolution,
    reconstructed_partition_id = as.character(cluster),
    n_cells = as.integer(n_cells),
    mean_best_match_jaccard = jaccard_mean,
    mean_jaccard_ge_0.75 = jaccard_mean >= 0.75,
    small_partition_le_150_cells = n_cells <= 150L
  ) |>
  arrange(leiden_resolution, as.integer(reconstructed_partition_id))

# Figure-ready tables.
seed_pairs_symmetric <- bind_rows(
  table_s2b |>
    transmute(seed_row = seed_a, seed_column = seed_b, adjusted_rand_index),
  table_s2b |>
    transmute(seed_row = seed_b, seed_column = seed_a, adjusted_rand_index),
  tibble(
    seed_row = seeds,
    seed_column = seeds,
    adjusted_rand_index = 1
  )
) |>
  arrange(match(seed_row, seeds), match(seed_column, seeds)) |>
  left_join(
    table_s2a |>
      select(
        seed_row = seed,
        n_partitions,
        adjusted_rand_index_vs_submitted
      ),
    by = "seed_row"
  )

fig_s2_subsampling_partition <- table_s2c |>
  mutate(
    display_partition = paste0(
      partition,
      " (", submitted_trace_label,
      if_else(partition == 6L, "; doublet-enriched", ""),
      ")"
    )
  )

fig_s2_resolution <- table_s2e |>
  transmute(
    resolution = leiden_resolution,
    n_partitions_full_cohort,
    modal_n_partitions_across_subsamples,
    mean_jaccard = mean_best_match_jaccard_across_full_cohort_partitions,
    median_jaccard = median_best_match_jaccard_across_full_cohort_partitions,
    minimum_jaccard = minimum_best_match_jaccard,
    maximum_jaccard = maximum_best_match_jaccard,
    n_partitions_ge_0.75 = n_partitions_mean_jaccard_ge_0.75,
    n_partitions_ge_0.60 = n_partitions_mean_jaccard_ge_0.60
  )

metadata <- tibble(
  submitted_cells = 7461L,
  stored_graph_resolution = 0.4,
  stored_graph_seeds = paste(seeds, collapse = ","),
  stored_graph_nondefault_ari_min = min(table_s2a$adjusted_rand_index_vs_submitted[-1]),
  stored_graph_nondefault_ari_max = max(table_s2a$adjusted_rand_index_vs_submitted[-1]),
  stored_graph_pairwise_ari_mean = mean(table_s2b$adjusted_rand_index),
  stored_graph_pairwise_ari_min = min(table_s2b$adjusted_rand_index),
  stored_graph_pairwise_ari_max = max(table_s2b$adjusted_rand_index),
  subsampling_fraction = 0.8,
  subsampling_cells = round(7461 * 0.8),
  subsampling_variable_seed_runs = 100L,
  subsampling_fixed_seed_runs = 50L,
  resolution_runs_per_setting = 30L,
  claim_ceiling = paste(
    "Algorithmic sensitivity and conditional reconstructed-graph diagnostics only;",
    "not biological replication, end-to-end stability, or state validation"
  )
)

data_dictionary <- tribble(
  ~part, ~filename, ~unit, ~description,
  "A", "TableS2A_fixed_stored_graph_seed_results", "Leiden seed", paste(
    "Partition count and adjusted Rand index against the submitted partition",
    "for six Leiden seeds on the exact stored neighborhood graph"
  ),
  "B", "TableS2B_fixed_stored_graph_pairwise_ari", "Unique seed pair", paste(
    "All 15 pairwise adjusted Rand indices among the same six fixed-graph runs"
  ),
  "C", "TableS2C_subsampling_submitted_partition_jaccard", "Submitted partition", paste(
    "Best-match Jaccard summaries after 80% cell subsampling and graph",
    "reconstruction from inherited PC coordinates at resolution 0.4"
  ),
  "D", "TableS2D_subsampling_cluster_count_distribution", "Arm x partition count", paste(
    "Recovered-partition count distribution for 100 variable-seed and 50",
    "fixed-seed reconstructed-graph subsampling runs"
  ),
  "E", "TableS2E_resolution_sweep_summary", "Leiden resolution", paste(
    "Full-cohort and modal subsample partition counts and across-partition",
    "best-match Jaccard summaries for eight reconstructed-graph resolutions"
  ),
  "F", "TableS2F_resolution_sweep_per_partition", "Reconstructed partition", paste(
    "Per-partition size and mean best-match Jaccard at each resolution; IDs",
    "are resolution-specific and do not map to submitted partition numbers"
  )
)

panel_output_paths <- c(
  metadata = file.path(panel_dir, "figS02_metadata.csv"),
  fixed_seed_matrix = file.path(panel_dir, "figS02_fixed_seed_matrix.csv"),
  subsampling_partition = file.path(
    panel_dir,
    "figS02_subsampling_partition.csv"
  ),
  cluster_count_distribution = file.path(
    panel_dir,
    "figS02_cluster_count_distribution.csv"
  ),
  resolution_summary = file.path(panel_dir, "figS02_resolution_summary.csv")
)

write_csv(metadata, panel_output_paths[["metadata"]])
write_csv(seed_pairs_symmetric, panel_output_paths[["fixed_seed_matrix"]])
write_csv(
  fig_s2_subsampling_partition,
  panel_output_paths[["subsampling_partition"]]
)
write_csv(table_s2d, panel_output_paths[["cluster_count_distribution"]])
write_csv(fig_s2_resolution, panel_output_paths[["resolution_summary"]])

table_bases <- c(
  A = "TableS2A_fixed_stored_graph_seed_results",
  B = "TableS2B_fixed_stored_graph_pairwise_ari",
  C = "TableS2C_subsampling_submitted_partition_jaccard",
  D = "TableS2D_subsampling_cluster_count_distribution",
  E = "TableS2E_resolution_sweep_summary",
  F = "TableS2F_resolution_sweep_per_partition",
  dictionary = "TableS2_data_dictionary"
)
table_data <- list(
  A = table_s2a,
  B = table_s2b,
  C = table_s2c,
  D = table_s2d,
  E = table_s2e,
  F = table_s2f,
  dictionary = data_dictionary
)

table_output_paths <- character()
for (key in names(table_bases)) {
  csv_path <- file.path(table_dir, paste0(table_bases[[key]], ".csv"))
  tsv_path <- file.path(table_dir, paste0(table_bases[[key]], ".tsv"))
  write_csv(table_data[[key]], csv_path)
  write_tsv(table_data[[key]], tsv_path)
  table_output_paths[paste0(key, "_csv")] <- csv_path
  table_output_paths[paste0(key, "_tsv")] <- tsv_path
}

source_roles <- c(
  "Correct seed sensitivity on the exact graph stored in the submitted object",
  "Correct pairwise seed ARI on the exact stored graph",
  "Independent read-only audit that generated the correct stored-graph outputs",
  paste(
    "Conditional per-submitted-partition Jaccard summaries after graph",
    "reconstruction from inherited PCs"
  ),
  "Conditional reconstructed-graph subsampling cluster counts",
  "Conditional reconstructed-graph resolution summary",
  "Conditional reconstructed-graph per-partition resolution values",
  "Provenance for the conditional subsampling outputs",
  "Provenance for the conditional resolution outputs"
)
source_manifest <- tibble(
  input = names(input_paths),
  path = unname(input_paths),
  role = source_roles,
  tests_submitted_stored_graph = c(TRUE, TRUE, TRUE, rep(FALSE, 6)),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write_csv(source_manifest, file.path(manifest_dir, "figS02_source_manifest.csv"))

all_output_paths <- c(panel_output_paths, table_output_paths)
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
  file.path(manifest_dir, "figS02_panel_table_output_manifest.csv")
)

execution <- list(
  figure = "Figure S2",
  supplementary_table = "Table S2",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  submitted_cells = 7461L,
  fixed_stored_graph = list(
    resolution = 0.4,
    seeds = as.list(seeds),
    output = "Parts A-B and Figure S2A",
    interpretation = paste(
      "Algorithmic Leiden-seed sensitivity on the exact graph stored in",
      "adata_microglia_subtyped.h5ad"
    )
  ),
  conditional_reconstructed_graph = list(
    pc_basis = "Inherited whole-cell PC1-PC20",
    n_neighbors = 15L,
    cell_subsampling_fraction = 0.8,
    resolution_0.4_runs = list(variable_seed = 100L, fixed_seed = 50L),
    resolution_sweep = list(
      settings = as.list(resolutions),
      runs_per_setting = 30L
    ),
    output = "Parts C-F and Figure S2B-C",
    interpretation = paste(
      "Conditional graph-reconstruction diagnostics; not the submitted",
      "stored graph, end-to-end stability, or biological replication"
    )
  ),
  excluded_as_submitted_graph_evidence = as.list(excluded_inputs),
  submitted_partition_labels = paste(
    "Numeric partitions are primary; retrospective biological names appear",
    "only as traceability labels"
  ),
  claim_ceiling = metadata$claim_ceiling[[1L]]
)
write_json(
  execution,
  file.path(manifest_dir, "figS02_data_preparation_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "figS02_data_preparation_sessionInfo.txt")
)

message("Figure S2 / Table S2 data preparation complete.")
message(sprintf(
  "Fixed stored graph: non-default seed ARI %.3f-%.3f; pairwise mean %.3f.",
  metadata$stored_graph_nondefault_ari_min,
  metadata$stored_graph_nondefault_ari_max,
  metadata$stored_graph_pairwise_ari_mean
))
message(
  "Conditional graph reconstruction: 100 variable-seed + 50 fixed-seed ",
  "subsamples; 30 subsamples at each of eight resolutions."
)
