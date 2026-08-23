#!/usr/bin/env Rscript

# GLIA major revision — canonical Figure S2 plotting and assembly
#
# Panel A uses the exact submitted stored graph. Panels B-C use graphs rebuilt
# from inherited whole-cell PC coordinates after cell subsampling. The panels
# are intentionally separated because they do not test the same computational
# object. All plotting, assembly, and publication export are performed in R.

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(ragg)
  library(readr)
  library(systemfonts)
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
panel_dir <- file.path(figure_root, "data", "panel_ready")
output_dir <- file.path(figure_root, "outputs")
manifest_dir <- file.path(figure_root, "manifests")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

width_mm <- as.numeric(get_arg("--width-mm", "180"))
height_mm <- as.numeric(get_arg("--height-mm", "160"))
dpi <- as.integer(get_arg("--dpi", "300"))
preview_suffix <- get_arg("--preview-suffix", "")
assert_true(
  width_mm == 180 && height_mm == 160,
  "Canonical Figure S2 dimensions are frozen at 180 x 160 mm."
)
assert_true(dpi == 300L, "Submission TIFF is frozen at 300 dpi.")
assert_true(
  grepl("^(_[A-Za-z0-9][A-Za-z0-9._-]*)?$", preview_suffix),
  "Invalid preview suffix."
)

input_paths <- c(
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
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing Figure S2 input(s):", paste(missing_inputs, collapse = ", "))
)

metadata <- read_csv(input_paths[["metadata"]], show_col_types = FALSE)
seed_matrix <- read_csv(
  input_paths[["fixed_seed_matrix"]],
  show_col_types = FALSE
)
subsampling_partition <- read_csv(
  input_paths[["subsampling_partition"]],
  show_col_types = FALSE
)
cluster_count_distribution <- read_csv(
  input_paths[["cluster_count_distribution"]],
  show_col_types = FALSE
)
resolution_summary <- read_csv(
  input_paths[["resolution_summary"]],
  show_col_types = FALSE
)

assert_columns(
  metadata,
  c(
    "submitted_cells", "stored_graph_resolution", "stored_graph_seeds",
    "stored_graph_nondefault_ari_min", "stored_graph_nondefault_ari_max",
    "stored_graph_pairwise_ari_mean", "subsampling_fraction",
    "subsampling_variable_seed_runs", "subsampling_fixed_seed_runs",
    "resolution_runs_per_setting", "claim_ceiling"
  ),
  "figS02_metadata.csv"
)
assert_columns(
  seed_matrix,
  c(
    "seed_row", "seed_column", "adjusted_rand_index", "n_partitions",
    "adjusted_rand_index_vs_submitted"
  ),
  "figS02_fixed_seed_matrix.csv"
)
assert_columns(
  subsampling_partition,
  c(
    "partition", "submitted_trace_label", "n_cells_in_submitted_partition",
    "mean_best_match_jaccard_variable_seed",
    "sd_best_match_jaccard_variable_seed",
    "mean_best_match_jaccard_fixed_seed", "display_partition"
  ),
  "figS02_subsampling_partition.csv"
)
assert_columns(
  cluster_count_distribution,
  c(
    "arm", "n_runs", "n_partitions_recovered", "n_runs_at_count",
    "percentage_of_runs"
  ),
  "figS02_cluster_count_distribution.csv"
)
assert_columns(
  resolution_summary,
  c(
    "resolution", "n_partitions_full_cohort",
    "modal_n_partitions_across_subsamples", "mean_jaccard",
    "median_jaccard", "minimum_jaccard", "maximum_jaccard",
    "n_partitions_ge_0.75"
  ),
  "figS02_resolution_summary.csv"
)

seeds <- c(0L, 1L, 7L, 42L, 123L, 2024L)
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
assert_true(
  nrow(metadata) == 1L &&
    metadata$submitted_cells == 7461L &&
    metadata$stored_graph_resolution == 0.4 &&
    metadata$subsampling_fraction == 0.8 &&
    metadata$subsampling_variable_seed_runs == 100L &&
    metadata$subsampling_fixed_seed_runs == 50L &&
    metadata$resolution_runs_per_setting == 30L,
  "Figure S2 metadata fail the frozen scope/run anchors."
)
assert_true(
  nrow(seed_matrix) == 36L &&
    setequal(seed_matrix$seed_row, seeds) &&
    setequal(seed_matrix$seed_column, seeds) &&
    nrow(distinct(seed_matrix, seed_row, seed_column)) == 36L &&
    all(seed_matrix$adjusted_rand_index >= 0) &&
    all(seed_matrix$adjusted_rand_index <= 1),
  "Figure S2A matrix is not the complete symmetric six-seed matrix."
)
assert_true(
  nrow(subsampling_partition) == 7L &&
    identical(as.integer(subsampling_partition$partition), 0:6) &&
    sum(subsampling_partition$n_cells_in_submitted_partition) == 7461L &&
    sum(
      subsampling_partition$mean_best_match_jaccard_variable_seed >= 0.75
    ) == 1L,
  "Figure S2B per-partition summaries fail the frozen anchors."
)
assert_true(
  nrow(cluster_count_distribution) == 10L &&
    all(
      cluster_count_distribution |>
        group_by(arm) |>
        summarise(total = sum(percentage_of_runs), .groups = "drop") |>
        pull(total) == 100
    ),
  "Figure S2B cluster-count distributions do not sum to 100%."
)
assert_true(
  nrow(resolution_summary) == 8L &&
    max(abs(resolution_summary$resolution - resolutions)) < 1e-12 &&
    identical(
      as.integer(resolution_summary$n_partitions_full_cohort),
      c(3L, 4L, 6L, 7L, 9L, 10L, 10L, 12L)
    ),
  "Figure S2C resolution summary fails the frozen anchors."
)

# Shared visual specification.
font_match <- systemfonts::match_fonts("Arial")
font_family <- if (nrow(font_match) > 0L && nzchar(font_match$path[[1L]])) {
  "Arial"
} else {
  "Helvetica"
}
if (font_family != "Arial") warning("Arial was not found; using Helvetica.")

purple <- "#6A51A3"
purple_dark <- "#4A3780"
purple_pale <- "#EFEAF7"
neutral_dark <- "#333333"
neutral_mid <- "#737373"
neutral_light <- "#D4D4D4"
neutral_faint <- "#ECECEC"
neutral_pale <- "#F7F7F7"

base_theme <- theme_classic(base_family = font_family, base_size = 8.2) +
  theme(
    text = element_text(color = neutral_dark),
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 7.8),
    axis.text = element_text(size = 7.0, color = neutral_dark),
    legend.title = element_text(size = 7.0),
    legend.text = element_text(size = 6.9),
    plot.margin = margin(6, 7, 5, 7)
  )

# Panel A: all pairwise adjusted Rand indices on the exact graph stored in the
# submitted object. Marginal columns report each run versus submitted labels.
seed_annotations <- seed_matrix |>
  distinct(seed_row, n_partitions, adjusted_rand_index_vs_submitted) |>
  mutate(
    y = length(seeds) + 1L - match(seed_row, seeds),
    ari_label = sprintf("%.2f", adjusted_rand_index_vs_submitted),
    k_label = as.character(n_partitions)
  )

panel_a_data <- seed_matrix |>
  mutate(
    x = match(seed_column, seeds),
    y = length(seeds) + 1L - match(seed_row, seeds),
    value_label = sprintf("%.2f", adjusted_rand_index),
    text_color = if_else(adjusted_rand_index >= 0.72, "white", neutral_dark)
  )

p_a <- ggplot(panel_a_data, aes(x = x, y = y)) +
  geom_tile(
    aes(fill = adjusted_rand_index),
    width = 0.94,
    height = 0.94,
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = value_label, color = text_color),
    family = font_family,
    size = 2.1,
    fontface = "bold"
  ) +
  geom_vline(
    xintercept = 6.62,
    color = neutral_light,
    linewidth = 0.45
  ) +
  geom_text(
    data = seed_annotations,
    aes(x = 7.25, y = y, label = ari_label),
    inherit.aes = FALSE,
    family = font_family,
    size = 2.15,
    color = neutral_dark
  ) +
  geom_text(
    data = seed_annotations,
    aes(x = 8.15, y = y, label = k_label),
    inherit.aes = FALSE,
    family = font_family,
    size = 2.15,
    color = neutral_dark
  ) +
  annotate(
    "text", x = 3.5, y = 6.72, label = "Pairwise ARI",
    family = font_family, size = 2.25, fontface = "bold",
    color = neutral_dark
  ) +
  annotate(
    "text", x = 7.25, y = 6.72, label = "ARI vs\nsubmitted",
    family = font_family, size = 2.05, fontface = "bold",
    lineheight = 0.85, color = neutral_dark
  ) +
  annotate(
    "text", x = 8.15, y = 6.72, label = "k",
    family = font_family, size = 2.1, fontface = "bold",
    color = neutral_dark
  ) +
  scale_fill_gradientn(
    colours = c("#F2F2F2", "#CFC3E2", "#8C6DB6", purple_dark),
    values = scales::rescale(c(0.35, 0.55, 0.75, 1.0)),
    limits = c(0.35, 1.0),
    guide = "none"
  ) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = seq_along(seeds),
    labels = seeds,
    limits = c(0.45, 8.45),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq_along(seeds),
    labels = rev(seeds),
    limits = c(0.45, 6.82),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(
    tag = "A",
    x = "Leiden seed (column)",
    y = "Leiden seed (row)"
  ) +
  base_theme +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(10, 5, 5, 8)
  )

# Panel B, upper: submitted-partition best-match Jaccard after graph
# reconstruction for each 80% cell subsample. Variable-seed error bars are SD
# across runs, not confidence intervals.
panel_b_labels <- subsampling_partition |>
  mutate(
    y = 7L - partition,
    variable_low = pmax(
      0,
      mean_best_match_jaccard_variable_seed -
        sd_best_match_jaccard_variable_seed
    ),
    variable_high = pmin(
      1,
      mean_best_match_jaccard_variable_seed +
        sd_best_match_jaccard_variable_seed
    ),
    n_label = paste0("n=", n_cells_in_submitted_partition),
    n_label_x = if_else(partition == 6L, 0.895, 1.035)
  )

p_b_top <- ggplot(panel_b_labels, aes(y = y)) +
  geom_vline(
    xintercept = 0.75,
    color = neutral_mid,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_segment(
    aes(
      x = mean_best_match_jaccard_fixed_seed,
      xend = mean_best_match_jaccard_variable_seed,
      yend = y
    ),
    color = neutral_light,
    linewidth = 1.0,
    lineend = "round"
  ) +
  geom_segment(
    aes(x = variable_low, xend = variable_high, yend = y),
    color = purple,
    linewidth = 0.65,
    lineend = "round"
  ) +
  geom_segment(
    aes(x = variable_low, xend = variable_low, y = y - 0.10, yend = y + 0.10),
    color = purple,
    linewidth = 0.55
  ) +
  geom_segment(
    aes(x = variable_high, xend = variable_high, y = y - 0.10, yend = y + 0.10),
    color = purple,
    linewidth = 0.55
  ) +
  geom_point(
    aes(x = mean_best_match_jaccard_fixed_seed),
    shape = 21,
    size = 2.5,
    stroke = 0.65,
    fill = "white",
    color = neutral_dark
  ) +
  geom_point(
    aes(x = mean_best_match_jaccard_variable_seed),
    shape = 21,
    size = 2.7,
    stroke = 0.65,
    fill = purple,
    color = neutral_dark
  ) +
  geom_text(
    aes(x = n_label_x, label = n_label),
    family = font_family,
    size = 1.95,
    hjust = 1,
    color = neutral_mid
  ) +
  annotate(
    "point", x = 0.23, y = 7.63, shape = 21, size = 2.3,
    stroke = 0.6, fill = purple, color = neutral_dark
  ) +
  annotate(
    "text", x = 0.26, y = 7.63, label = "Seeds varied (100)",
    hjust = 0, family = font_family, size = 1.95,
    color = neutral_dark
  ) +
  annotate(
    "point", x = 0.57, y = 7.63, shape = 21, size = 2.3,
    stroke = 0.6, fill = "white", color = neutral_dark
  ) +
  annotate(
    "text", x = 0.60, y = 7.63, label = "Seeds fixed (50)",
    hjust = 0, family = font_family, size = 1.95,
    color = neutral_dark
  ) +
  annotate(
    "text", x = 0.75, y = 0.35, label = "0.75",
    family = font_family, size = 1.85, hjust = 0.5,
    color = neutral_mid
  ) +
  scale_x_continuous(
    breaks = c(0.25, 0.50, 0.75, 1.00),
    limits = c(0.20, 1.05),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = panel_b_labels$y,
    labels = panel_b_labels$display_partition,
    limits = c(0.25, 7.78),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "B",
    x = "Best-match Jaccard at resolution 0.4",
    y = NULL
  ) +
  base_theme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    plot.margin = margin(10, 7, 2, 9)
  )

# Panel B, lower: recovered partition-count distribution for the same
# reconstructed-graph subsampling arms.
arm_levels <- c(
  "Graph and Leiden seeds varied",
  "Graph and Leiden seeds fixed"
)
panel_b_counts <- cluster_count_distribution |>
  mutate(
    arm = factor(arm, levels = arm_levels),
    value_label = if_else(
      percentage_of_runs > 0,
      sprintf("%.0f", percentage_of_runs),
      ""
    )
  )

p_b_bottom <- ggplot(
  panel_b_counts,
  aes(x = factor(n_partitions_recovered), y = percentage_of_runs)
) +
  geom_col(
    aes(fill = arm, color = arm),
    position = position_dodge(width = 0.72),
    width = 0.66,
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = value_label, group = arm),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    family = font_family,
    size = 1.85,
    color = neutral_dark
  ) +
  scale_fill_manual(
    values = c(
      "Graph and Leiden seeds varied" = purple,
      "Graph and Leiden seeds fixed" = "white"
    ),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(
      "Graph and Leiden seeds varied" = purple_dark,
      "Graph and Leiden seeds fixed" = neutral_mid
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = c(0, 25, 50),
    limits = c(0, 60),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Partitions recovered per 80% subsample",
    y = "Runs (%)"
  ) +
  base_theme +
  theme(
    panel.grid.major.y = element_line(color = neutral_faint, linewidth = 0.35),
    plot.margin = margin(2, 7, 5, 9)
  )

# Panel C: across-partition summary at eight Leiden resolutions. Vertical
# ranges are the minimum-to-maximum per-partition mean Jaccard, not uncertainty
# intervals. Per-partition values remain available in Table S2F.
panel_c <- resolution_summary |>
  mutate(
    x = seq_len(n()),
    x_label = paste0(
      format(resolution, trim = TRUE, nsmall = if_else(resolution == 1, 1L, 1L)),
      "\n",
      n_partitions_full_cohort,
      "/",
      modal_n_partitions_across_subsamples
    ),
    stable_label = paste0(
      n_partitions_ge_0.75,
      "/",
      n_partitions_full_cohort
    )
  )

panel_c_long <- panel_c |>
  select(x, mean_jaccard, median_jaccard) |>
  pivot_longer(
    cols = c(mean_jaccard, median_jaccard),
    names_to = "summary",
    values_to = "jaccard"
  ) |>
  mutate(
    summary = factor(
      summary,
      levels = c("mean_jaccard", "median_jaccard"),
      labels = c("Mean", "Median")
    )
  )

p_c <- ggplot() +
  annotate(
    "rect", xmin = 3.63, xmax = 4.37, ymin = 0, ymax = 1.14,
    fill = purple_pale, color = NA
  ) +
  geom_hline(
    yintercept = 0.75,
    color = neutral_mid,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_segment(
    data = panel_c,
    aes(x = x, xend = x, y = minimum_jaccard, yend = maximum_jaccard),
    color = neutral_light,
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_line(
    data = panel_c_long,
    aes(x = x, y = jaccard, color = summary, linetype = summary),
    linewidth = 0.65
  ) +
  geom_point(
    data = panel_c_long,
    aes(x = x, y = jaccard, fill = summary, shape = summary),
    size = 2.7,
    stroke = 0.65,
    color = neutral_dark
  ) +
  geom_text(
    data = panel_c,
    aes(x = x, y = 1.055, label = stable_label),
    family = font_family,
    size = 2.0,
    color = neutral_dark
  ) +
  annotate(
    "text", x = 0.65, y = 1.105,
    label = "Partitions with mean Jaccard ≥0.75 / full-cohort k",
    hjust = 0, family = font_family, size = 2.05,
    fontface = "bold", color = neutral_dark
  ) +
  annotate(
    "text", x = 4, y = 0.055,
    label = "submitted\nresolution",
    family = font_family, size = 1.85, lineheight = 0.86,
    color = purple_dark
  ) +
  scale_color_manual(
    values = c(Mean = purple, Median = neutral_dark),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(Mean = purple, Median = "white"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(Mean = 21, Median = 21),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(Mean = "solid", Median = "22"),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = panel_c$x,
    labels = panel_c$x_label,
    limits = c(0.55, 8.45),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    limits = c(0, 1.14),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "C",
    x = paste0(
      "Leiden resolution\n",
      "full-cohort k / modal subsample k"
    ),
    y = "Across-partition best-match Jaccard"
  ) +
  base_theme +
  theme(
    panel.grid.major.y = element_line(color = neutral_faint, linewidth = 0.35),
    legend.position = c(0.93, 0.91),
    legend.direction = "horizontal",
    legend.key.width = grid::unit(3.2, "mm"),
    legend.spacing.x = grid::unit(1.0, "mm"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 7, 5, 8)
  )

p_b <- p_b_top / p_b_bottom + plot_layout(heights = c(1.52, 0.72))
top_row <- p_a + p_b + plot_layout(widths = c(0.96, 1.04))
figure_s2 <- top_row / p_c + plot_layout(heights = c(1.22, 0.78))

# Export and execution records.
base_name <- "FigureS2_partition_stability_diagnostics"
pdf_path <- file.path(output_dir, paste0(base_name, ".pdf"))
tiff_path <- file.path(output_dir, paste0(base_name, ".tiff"))
png_path <- file.path(
  output_dir,
  paste0(base_name, "_preview", preview_suffix, ".png")
)
rds_path <- file.path(output_dir, paste0(base_name, ".rds"))

width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
render_warnings <- character()
capture_render_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      render_warnings <<- unique(c(render_warnings, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
}

save_quartz_pdf <- function(plot, filename, width, height, family) {
  assert_true(capabilities("aqua"), "Quartz PDF export requires macOS Aqua.")
  grDevices::quartz(
    file = filename,
    type = "pdf",
    width = width,
    height = height,
    pointsize = 8,
    family = family,
    bg = "white",
    antialias = TRUE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  capture_render_warnings(print(plot))
}

save_quartz_pdf(figure_s2, pdf_path, width_in, height_in, font_family)

ragg::agg_tiff(
  filename = tiff_path,
  width = width_mm,
  height = height_mm,
  units = "mm",
  pointsize = 8,
  background = "white",
  res = dpi,
  scaling = 1,
  compression = "lzw",
  bitsize = 8
)
capture_render_warnings(print(figure_s2))
grDevices::dev.off()

ragg::agg_png(
  filename = png_path,
  width = width_mm,
  height = height_mm,
  units = "mm",
  pointsize = 8,
  background = "white",
  res = dpi,
  scaling = 1
)
capture_render_warnings(print(figure_s2))
grDevices::dev.off()

saveRDS(figure_s2, rds_path, compress = "xz")
output_paths <- c(
  vector_pdf_master = pdf_path,
  glia_submission_tiff = tiff_path,
  png_preview = png_path,
  plot_object = rds_path
)

writeLines(
  if (length(render_warnings)) render_warnings else "None",
  file.path(manifest_dir, "figS02_render_warnings.txt")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "figS02_plot_sessionInfo.txt")
)
write_csv(
  tibble(
    input = names(input_paths),
    relative_path = file.path("data", "panel_ready", basename(input_paths)),
    sha256 = vapply(input_paths, sha256_file, character(1))
  ),
  file.path(manifest_dir, "figS02_plot_input_manifest.csv")
)
write_csv(
  tibble(
    output = names(output_paths),
    relative_path = file.path("outputs", basename(output_paths)),
    sha256 = vapply(output_paths, sha256_file, character(1))
  ),
  file.path(manifest_dir, "figS02_plot_output_manifest.csv")
)

execution <- list(
  figure = "Figure S2",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  dimensions = list(width_mm = width_mm, height_mm = height_mm),
  panels = list(
    A = paste(
      "Correct fixed stored-graph pairwise seed ARI, with ARI versus",
      "submitted labels and recovered partition count"
    ),
    B = paste(
      "Conditional resolution-0.4 80% cell-subsampling Jaccard summaries",
      "and recovered partition-count distributions after graph reconstruction"
    ),
    C = paste(
      "Conditional reconstructed-graph resolution sweep, with across-partition",
      "mean, median, and min-max Jaccard summaries"
    )
  ),
  fixed_stored_graph_panels = "A only",
  reconstructed_graph_panels = "B-C",
  biological_replication = FALSE,
  end_to_end_pipeline_stability = FALSE,
  state_validation = FALSE,
  claim_ceiling = metadata$claim_ceiling[[1L]],
  submission_raster = list(
    format = "TIFF",
    dpi = dpi,
    compression = "LZW",
    bitsize = 8,
    background = "white"
  ),
  vector_master = list(
    format = "PDF",
    device = "macOS Quartz PDF",
    font_family = font_family
  ),
  plotting = "R only; no manual artwork editing"
)
write_json(
  execution,
  file.path(manifest_dir, "figS02_plot_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Figure S2 rendering complete.")
message("PDF: ", pdf_path)
message("TIFF: ", tiff_path)
message("PNG: ", png_path)
