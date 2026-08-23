#!/usr/bin/env Rscript

# GLIA major revision — canonical Figure 1 plotting script
#
# All panels, annotations, assembly, and exports are generated in R.
# Input: validated panel-ready CSV files from script 78.
# Output: embedded-font vector PDF master, 300-dpi LZW TIFF submission file,
#         PNG preview, RDS plot object, hashes, warnings, and session record.
#
# Example:
# Rscript scripts/79_plot_glia_R1_figure1.R \
#   --figure-root /path/to/figure1

suppressPackageStartupMessages({
  library(colorspace)
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(ragg)
  library(readr)
  library(scales)
  library(systemfonts)
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
    if (idx == length(args)) {
      stop("Missing value after ", flag)
    }
    return(args[[idx + 1L]])
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
figure_root <- normalizePath(
  get_arg("--figure-root", default_figure_root),
  mustWork = TRUE
)

panel_dir <- file.path(figure_root, "data", "panel_ready")
manifest_dir <- file.path(figure_root, "manifests")
output_dir <- file.path(figure_root, "outputs")
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

width_mm <- as.numeric(get_arg("--width-mm", "180"))
height_mm <- as.numeric(get_arg("--height-mm", "215"))
dpi <- as.integer(get_arg("--dpi", "300"))
preview_suffix <- get_arg("--preview-suffix", "")
assert_true(is.finite(width_mm) && width_mm > 0, "Invalid --width-mm.")
assert_true(is.finite(height_mm) && height_mm > 0, "Invalid --height-mm.")
assert_true(dpi == 300L, "GLIA revision submission TIFF is frozen at 300 dpi.")
assert_true(
  grepl("^(_[A-Za-z0-9][A-Za-z0-9._-]*)?$", preview_suffix),
  "Invalid --preview-suffix; use an empty value or a leading underscore followed by letters, numbers, dots, underscores, or hyphens."
)

input_paths <- c(
  design = file.path(panel_dir, "fig01_design.csv"),
  analysis_sets = file.path(panel_dir, "fig01_analysis_sets.csv"),
  sample_qc = file.path(panel_dir, "fig01_sample_qc.csv"),
  umap = file.path(panel_dir, "fig01_umap.csv"),
  seed_stability = file.path(panel_dir, "fig01_seed_stability.csv"),
  composition = file.path(panel_dir, "fig01_composition.csv"),
  composition_effects = file.path(panel_dir, "fig01_composition_effects.csv")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing panel-ready inputs:", paste(basename(missing_inputs), collapse = ", "))
)

design <- read_csv(input_paths[["design"]], show_col_types = FALSE)
analysis_sets <- read_csv(input_paths[["analysis_sets"]], show_col_types = FALSE)
sample_qc <- read_csv(input_paths[["sample_qc"]], show_col_types = FALSE)
umap <- read_csv(input_paths[["umap"]], show_col_types = FALSE)
seed_stability <- read_csv(input_paths[["seed_stability"]], show_col_types = FALSE)
composition <- read_csv(input_paths[["composition"]], show_col_types = FALSE)
composition_effects <- read_csv(
  input_paths[["composition_effects"]],
  show_col_types = FALSE
)

# Revalidate the critical anchors before plotting.
assert_true(nrow(umap) == 7461L, "Figure 1B must show 7,461 cells.")
assert_true(sum(sample_qc$predicted_doublets_removed) == 90L, "Expected 90 predicted doublets.")
assert_true(sum(sample_qc$cells_retained) == 7371L, "Expected 7,371 retained cells.")
assert_true(sum(composition$n_cells) == 7367L, "Figure 1D must use 7,367 non-Rare cells.")
assert_true(
  identical(as.integer(seed_stability$seed), c(0L, 1L, 7L, 42L, 123L, 2024L)),
  "Unexpected Figure 1C seed sequence."
)
assert_true(
  max(abs(
    seed_stability$ari_vs_original -
      c(1, 0.49092414553881564, 0.5858529047117915,
        0.5127939875974998, 0.4033111527516049,
        0.5165138553705313)
  )) < 1e-12,
  "Figure 1C ARI values differ from the verified anchors."
)

font_match <- systemfonts::match_fonts("Arial")
font_family <- if (nrow(font_match) > 0L && nzchar(font_match$path[[1L]])) {
  "Arial"
} else {
  "Helvetica"
}
if (font_family != "Arial") {
  warning("Arial was not found; using Helvetica.")
}

condition_colors <- c(
  "Oxygen control" = "#0072B2",
  "Combined exposure" = "#D55E00"
)
partition_colors <- c(
  "0" = "#EE6677",
  "1" = "#AA3377",
  "2" = "#228833",
  "3" = "#66CCEE",
  "4" = "#4477AA",
  "5" = "#CCBB44",
  "6" = "#BBBBBB"
)
neutral_dark <- "#333333"
neutral_mid <- "#6B6B6B"
neutral_light <- "#D6D6D6"
seed_color <- "#6A51A3"

base_theme <- theme_classic(base_family = font_family, base_size = 8.2) +
  theme(
    text = element_text(color = neutral_dark),
    plot.title = element_text(
      size = 9, face = "bold", hjust = 0,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 7.5, color = neutral_mid, hjust = 0,
      margin = margin(b = 4)
    ),
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7.5, color = neutral_dark),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7.2),
    legend.key.height = grid::unit(3.2, "mm"),
    legend.key.width = grid::unit(7, "mm"),
    plot.margin = margin(5.5, 5.5, 4.5, 5.5)
  )

# ---- Panel A: study design and analysis-set branching ----------------------

control_design <- design |> filter(source_group == "Control")
exposure_design <- design |> filter(source_group == "Surgery")
set_n <- setNames(analysis_sets$n_cells, analysis_sets$analysis_set)

panel_a_boxes <- tibble::tribble(
  ~box, ~xmin, ~xmax, ~ymin, ~ymax, ~fill, ~border, ~label,
  "control", 1, 19, 5.2, 9.9, "#E7F2F8", condition_colors[["Oxygen control"]],
  paste0(
    "Oxygen control\n50% O", "\u2082", ", 30 min\n",
    control_design$samples, " (n = ", control_design$n_animals, " animals)"
  ),
  "exposure", 1, 19, 0.1, 4.8, "#FBEDE5", condition_colors[["Combined exposure"]],
  paste0(
    "Combined exposure\n2.5% sevoflurane\nin 50% O", "\u2082",
    ", 30 min\n+ laparotomy\n",
    exposure_design$samples, " (n = ", exposure_design$n_animals, " animals)"
  ),
  "cohort", 25, 41, 2.5, 7.5, "#F5F5F5", "#8A8A8A",
  "18-month-old\nmale mice\nHippocampus\ncollected 24 h",
  "assay", 46, 59, 2.5, 7.5, "#F5F5F5", "#8A8A8A",
  paste0("10x Chromium\nSingle Cell\n3", "\u2032", " v2"),
  "microglia", 64, 76, 2.5, 7.5, "#F5F5F5", "#555555",
  paste0("Frozen\nsubmitted\nmicroglia\nn = ", comma(set_n[["Submitted-partition audit"]])),
  "partition", 82, 99, 5.2, 9.9, "#F2EDF8", seed_color,
  paste0(
    "Submitted-partition\naudit\nB/C: n = ",
    comma(set_n[["Submitted-partition audit"]]),
    "\nD: n = ",
    comma(set_n[["Conditional composition"]]),
    "\npartition 6 excluded"
  ),
  "molecular", 82, 99, 0.1, 4.8, "#EEF5EC", "#4D7A4B",
  paste0(
    "Molecular analyses\n90 Scrublet-\npredicted doublets\nexcluded\nn = ",
    comma(set_n[["Primary molecular analyses"]])
  )
)

panel_a_arrows <- tibble::tribble(
  ~x, ~xend, ~y, ~yend,
  19.5, 24.4, 7.4, 5.8,
  19.5, 24.4, 2.6, 4.2,
  41.5, 45.4, 5.0, 5.0,
  59.5, 63.4, 5.0, 5.0,
  76.5, 81.4, 5.8, 7.4,
  76.5, 81.4, 4.2, 2.6
)

p_a_body <- ggplot() +
  geom_rect(
    data = panel_a_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = panel_a_boxes$fill,
    color = panel_a_boxes$border,
    linewidth = 0.55
  ) +
  geom_segment(
    data = panel_a_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = neutral_mid,
    linewidth = 0.5,
    lineend = "round",
    arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed")
  ) +
  geom_text(
    data = panel_a_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    family = font_family,
    size = 2.65,
    lineheight = 0.86,
    color = neutral_dark
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 10), clip = "off") +
  theme_void(base_family = font_family, base_size = 8.2) +
  theme(
    plot.margin = margin(0, 5.5, 2, 5.5)
  )

p_a_header <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(tag = "A") +
  theme_void(base_family = font_family) +
  theme(
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(4, 5.5, 1, 5.5)
  )

p_a <- p_a_header / p_a_body +
  plot_layout(heights = c(0.07, 0.93))
p_a_composite <- wrap_elements(full = p_a)

# ---- Panel B: submitted numeric partition ---------------------------------

set.seed(20260728)
umap_plot <- umap[sample.int(nrow(umap)), ]
partition_centers <- umap |>
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    .by = partition
  )

p_b <- ggplot(
  umap_plot,
  aes(x = UMAP1, y = UMAP2, color = factor(partition))
) +
  geom_point(size = 0.27, alpha = 0.68, stroke = 0) +
  geom_label(
    data = partition_centers,
    aes(x = UMAP1, y = UMAP2, label = partition),
    inherit.aes = FALSE,
    family = font_family,
    fontface = "bold",
    size = 3.1,
    color = "#111111",
    fill = alpha("white", 0.88),
    linewidth = 0.25,
    label.padding = grid::unit(1.3, "mm"),
    label.r = grid::unit(1.7, "mm")
  ) +
  scale_color_manual(values = partition_colors, guide = "none") +
  coord_equal() +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  labs(
    tag = "B",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  base_theme +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 0.4, color = neutral_dark),
    plot.margin = margin(5.5, 8, 5, 7)
  )

# ---- Panel C: fixed stored-graph seed sensitivity --------------------------

seed_plot <- seed_stability |>
  mutate(
    seed_label = paste0("Seed ", seed, "  (k=", n_clusters, ")"),
    seed_label = factor(
      seed_label,
      levels = rev(paste0("Seed ", seed, "  (k=", n_clusters, ")"))
    ),
    point_fill = if_else(submitted_seed, "white", seed_color),
    point_color = if_else(submitted_seed, "#111111", seed_color)
  )

p_c <- ggplot(seed_plot, aes(y = seed_label)) +
  geom_segment(
    aes(x = 0, xend = ari_vs_original, yend = seed_label),
    linewidth = 0.65,
    color = neutral_light,
    lineend = "round"
  ) +
  geom_point(
    aes(x = ari_vs_original, fill = point_fill, color = point_color),
    shape = 21,
    size = 2.6,
    stroke = 0.7,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = ari_vs_original, label = sprintf("%.3f", ari_vs_original)),
    family = font_family,
    size = 2.7,
    hjust = -0.22,
    color = neutral_dark
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(
    limits = c(0, 1.18),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "C",
    x = "Adjusted Rand index vs submitted partition",
    y = NULL
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = "#ECECEC", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 15, 5, 9)
  )

# ---- Panel D: animal-level conditional composition -------------------------

sample_levels <- c("C1", "C2", "C3", "S1", "S2", "S3")
display_group_levels <- c("Oxygen control", "Combined exposure")
partition_levels_for_y <- as.character(5:0)

composition_plot <- composition |>
  mutate(
    sample = factor(sample, levels = sample_levels),
    display_group = factor(display_group, levels = display_group_levels),
    partition_factor = factor(as.character(partition), levels = partition_levels_for_y),
    text_color = if_else(pct_of_animal >= 27, "white", "#222222")
  )

heat_colors <- c("#FBF9FC", "#E7D4EC", "#C6A0D2", "#9866A8", "#683B7B")

p_d_heat <- ggplot(
  composition_plot,
  aes(x = sample, y = partition_factor, fill = pct_of_animal)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.1f", pct_of_animal), color = text_color),
    family = font_family,
    size = 2.55,
    show.legend = FALSE
  ) +
  facet_grid(
    cols = vars(display_group),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradientn(
    colors = heat_colors,
    limits = c(0, 45),
    oob = squish,
    breaks = c(0, 15, 30, 45),
    name = "Within-animal\npercentage"
  ) +
  scale_color_identity() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(
    labels = function(x) paste("Partition", x),
    expand = expansion(add = 0.55)
  ) +
  guides(fill = "none") +
  labs(
    x = NULL,
    y = NULL
  ) +
  base_theme +
  theme(
    panel.spacing.x = grid::unit(2, "mm"),
    strip.background = element_rect(fill = "#F1F1F1", color = NA),
    strip.text = element_text(size = 7.5, face = "bold", color = neutral_dark),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.5, face = "bold", color = neutral_dark),
    legend.position = "none",
    plot.margin = margin(4, 5, 5, 7)
  )

effects_plot <- composition_effects |>
  mutate(
    partition_factor = factor(
      as.character(partition),
      levels = partition_levels_for_y
    ),
    p_label = paste0("P = ", sprintf("%.2f", exact_permutation_p))
  )

p_d_forest <- ggplot(effects_plot, aes(y = partition_factor)) +
  geom_vline(
    xintercept = 0,
    color = "#777777",
    linewidth = 0.45,
    linetype = "22"
  ) +
  geom_segment(
    aes(x = ci95_low, xend = ci95_high, yend = partition_factor),
    linewidth = 0.75,
    color = neutral_dark,
    lineend = "round"
  ) +
  geom_point(
    aes(x = difference_percentage_points),
    shape = 21,
    size = 2.5,
    stroke = 0.65,
    fill = "white",
    color = neutral_dark
  ) +
  geom_text(
    aes(x = 37, label = p_label),
    family = font_family,
    size = 2.5,
    hjust = 0,
    color = neutral_mid
  ) +
  scale_x_continuous(
    limits = c(-35, 52),
    breaks = c(-30, -15, 0, 15, 30),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = 0.55)) +
  labs(
    x = "Combined exposure \u2212 oxygen control (percentage points)",
    y = NULL
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = "#ECECEC", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(4, 18, 5, 5)
  )

p_d_header <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(tag = "D") +
  theme_void(base_family = font_family) +
  theme(
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(4, 5.5, 1, 5.5)
  )

panel_d_body <- p_d_heat + p_d_forest +
  plot_layout(widths = c(1.28, 1))

panel_d <- p_d_header / panel_d_body +
  plot_layout(heights = c(0.05, 0.95))

figure_1 <- p_a_composite /
  (p_b + p_c + plot_layout(widths = c(1, 1))) /
  panel_d +
  plot_layout(heights = c(0.28, 0.32, 0.40))

# ---- Export and execution record ------------------------------------------

base_name <- "Figure1_study_design_partition_audit"
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
  assert_true(capabilities("aqua"), "Quartz PDF export requires macOS Aqua support.")
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

save_quartz_pdf(figure_1, pdf_path, width_in, height_in, font_family)

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
capture_render_warnings(print(figure_1))
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
capture_render_warnings(print(figure_1))
grDevices::dev.off()

saveRDS(figure_1, rds_path, compress = "xz")

output_paths <- c(
  vector_pdf_master = pdf_path,
  glia_submission_tiff = tiff_path,
  png_preview = png_path,
  plot_object = rds_path
)

writeLines(
  if (length(render_warnings)) render_warnings else "None",
  file.path(manifest_dir, "fig01_render_warnings.txt")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "fig01_plot_sessionInfo.txt")
)

input_manifest <- data.frame(
  input = names(input_paths),
  relative_path = file.path("data", "panel_ready", basename(input_paths)),
  sha256 = vapply(input_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
output_manifest <- data.frame(
  output = names(output_paths),
  relative_path = file.path("outputs", basename(output_paths)),
  sha256 = vapply(output_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv(input_manifest, file.path(manifest_dir, "fig01_plot_input_manifest.csv"))
write_csv(output_manifest, file.path(manifest_dir, "fig01_plot_output_manifest.csv"))

execution_manifest <- list(
  figure = "Figure 1",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  dimensions = list(width_mm = width_mm, height_mm = height_mm),
  preview_suffix = preview_suffix,
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
    font_family = font_family,
    fonts_requested_embedded = TRUE
  ),
  random_seed = 20260728L,
  render_warnings = as.list(render_warnings),
  r_version = R.version.string,
  packages = as.list(c(
    colorspace = as.character(packageVersion("colorspace")),
    digest = as.character(packageVersion("digest")),
    dplyr = as.character(packageVersion("dplyr")),
    ggplot2 = as.character(packageVersion("ggplot2")),
    jsonlite = as.character(packageVersion("jsonlite")),
    patchwork = as.character(packageVersion("patchwork")),
    ragg = as.character(packageVersion("ragg")),
    readr = as.character(packageVersion("readr")),
    scales = as.character(packageVersion("scales")),
    systemfonts = as.character(packageVersion("systemfonts"))
  )),
  outputs = split(output_manifest, seq_len(nrow(output_manifest)))
)
write_json(
  execution_manifest,
  file.path(manifest_dir, "fig01_plot_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Figure 1 written to: ", output_dir)
message("PDF master: ", basename(pdf_path))
message("TIFF submission file: ", basename(tiff_path), " (", dpi, " dpi, LZW)")
message("PNG preview: ", basename(png_path))
