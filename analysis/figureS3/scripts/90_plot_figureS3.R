#!/usr/bin/env Rscript

# GLIA major revision — canonical Figure S3 plotting script
#
# GSE289098 is displayed only as a same-six-library processed-count
# sensitivity analysis.  No clustering, subtype, per-cell inference, or
# independent-replication display is generated here.

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(patchwork)
  library(ragg)
  library(readr)
  library(systemfonts)
})

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg) != 1L) stop("Unable to resolve script path.")
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
height_mm <- as.numeric(get_arg("--height-mm", "180"))
dpi <- as.integer(get_arg("--dpi", "300"))
preview_suffix <- get_arg("--preview-suffix", "")
assert_true(
  width_mm == 180 && height_mm == 180,
  "Canonical Figure S3 dimensions are frozen at 180 x 180 mm."
)
assert_true(dpi == 300L, "Submission TIFF is frozen at 300 dpi.")
assert_true(
  grepl("^(_[A-Za-z0-9][A-Za-z0-9._-]*)?$", preview_suffix),
  "Invalid preview suffix."
)

input_paths <- c(
  metadata = file.path(panel_dir, "figS03_metadata.csv"),
  library_audit = file.path(panel_dir, "figS03_library_audit.csv"),
  animal_panel_values = file.path(panel_dir, "figS03_animal_panel_values.csv"),
  gene_delta_heatmap = file.path(panel_dir, "figS03_gene_delta_heatmap.csv"),
  panel_effects_full_loo = file.path(
    panel_dir,
    "figS03_panel_effects_full_loo.csv"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing panel-ready input(s):", paste(missing_inputs, collapse = ", "))
)

metadata <- read_csv(input_paths[["metadata"]], show_col_types = FALSE)
library_audit <- read_csv(input_paths[["library_audit"]], show_col_types = FALSE)
animal_panel <- read_csv(
  input_paths[["animal_panel_values"]],
  show_col_types = FALSE
)
gene_delta <- read_csv(
  input_paths[["gene_delta_heatmap"]],
  show_col_types = FALSE
)
panel_effects <- read_csv(
  input_paths[["panel_effects_full_loo"]],
  show_col_types = FALSE
)

assert_columns(
  metadata,
  c(
    "accession", "role", "common_estimand_build_gate",
    "panel_concordance_gate", "n_common_all_cells",
    "n_common_scrublet_negative_microglia", "n_common_features",
    "full_primary_panel_difference", "full_alternative_panel_difference",
    "alternative_minus_primary_panel_difference", "claim_ceiling"
  ),
  "figS03_metadata.csv"
)
assert_columns(
  library_audit,
  c(
    "sample", "source_group", "gse267933_gsm",
    "gse289098_barcode_suffix", "primary_all_cells",
    "integrated_mapped_cells", "n_common_scrublet_negative_cells",
    "primary_total_umi_common_cells", "alternative_total_umi_common_cells",
    "alternative_to_primary_umi_ratio_common_cells"
  ),
  "figS03_library_audit.csv"
)
assert_columns(
  animal_panel,
  c(
    "sample", "source_group", "primary_value", "alternative_value",
    "alternative_minus_primary"
  ),
  "figS03_animal_panel_values.csv"
)
assert_columns(
  gene_delta,
  c(
    "gene", "gene_order", "sample", "source_group",
    "alternative_minus_primary_log2_cpm"
  ),
  "figS03_gene_delta_heatmap.csv"
)
assert_columns(
  panel_effects,
  c(
    "scenario", "scenario_order", "primary_mean_difference",
    "gse289098_mean_difference", "alternative_minus_primary_difference",
    "sign_concordant"
  ),
  "figS03_panel_effects_full_loo.csv"
)

samples <- c("C1", "C2", "C3", "S1", "S2", "S3")
genes <- c("Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3")
assert_true(
  nrow(metadata) == 1L &&
    metadata$accession == "GSE289098" &&
    metadata$common_estimand_build_gate == "PASS" &&
    metadata$panel_concordance_gate == "PASS" &&
    metadata$n_common_all_cells == 20684L &&
    metadata$n_common_features == 27998L &&
    metadata$n_common_scrublet_negative_microglia == 7371L,
  "Figure S3 metadata fail the frozen identity/gate anchors."
)
assert_true(
  nrow(library_audit) == 6L &&
    identical(library_audit$sample, samples) &&
    all(library_audit$primary_all_cells == library_audit$integrated_mapped_cells),
  "Figure S3A library mapping is incomplete."
)
assert_true(
  nrow(animal_panel) == 6L &&
    identical(animal_panel$sample, samples) &&
    abs(
      mean(animal_panel$alternative_value[animal_panel$source_group == "Surgery"]) -
        mean(animal_panel$alternative_value[animal_panel$source_group == "Control"]) -
        1.0982960515733415
    ) < 1e-12,
  "Figure S3B animal values fail the verified effect anchor."
)
assert_true(
  nrow(gene_delta) == 42L &&
    nrow(distinct(gene_delta, gene, sample)) == 42L &&
    setequal(gene_delta$gene, genes) &&
    max(abs(gene_delta$alternative_minus_primary_log2_cpm)) < 0.22,
  "Figure S3C gene-by-animal sensitivity matrix is invalid."
)
assert_true(
  nrow(panel_effects) == 7L &&
    all(panel_effects$sign_concordant) &&
    all(panel_effects$primary_mean_difference > 0) &&
    all(panel_effects$gse289098_mean_difference > 0) &&
    abs(
      panel_effects$alternative_minus_primary_difference[
        panel_effects$scenario == "Full cohort"
      ] - 0.02972237625009333
    ) < 1e-12,
  "Figure S3D full/LOO sensitivity anchors are invalid."
)

# Shared visual specification.
font_match <- systemfonts::match_fonts("Arial")
font_family <- if (nrow(font_match) > 0L && nzchar(font_match$path[[1L]])) {
  "Arial"
} else {
  "Helvetica"
}
if (font_family != "Arial") warning("Arial was not found; using Helvetica.")

control_color <- "#0072B2"
exposure_color <- "#D55E00"
alternative_color <- "#6A51A3"
primary_color <- "#3F3F3F"
neutral_dark <- "#333333"
neutral_mid <- "#737373"
neutral_light <- "#D4D4D4"
neutral_faint <- "#ECECEC"
pale_neutral <- "#F5F5F5"

group_colors <- c(
  Control = control_color,
  Surgery = exposure_color
)
group_labels <- c(
  Control = "Oxygen control",
  Surgery = "Combined exposure"
)

base_theme <- theme_classic(base_family = font_family, base_size = 8.2) +
  theme(
    text = element_text(color = neutral_dark),
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 7.8),
    axis.text = element_text(size = 7.2, color = neutral_dark),
    strip.text = element_text(size = 7.2, face = "bold"),
    legend.title = element_text(size = 7.0),
    legend.text = element_text(size = 6.9),
    plot.margin = margin(6, 7, 5, 7)
  )

# Panel A: exact mapping and retained count payload by animal.
panel_a_data <- library_audit |>
  mutate(
    y = rev(seq_along(samples)),
    umi_retained_percent = 100 * alternative_to_primary_umi_ratio_common_cells,
    label = sprintf("%.1f", umi_retained_percent),
    label_x = if_else(
      umi_retained_percent > 99.65,
      umi_retained_percent - 0.10,
      umi_retained_percent + 0.10
    ),
    label_hjust = if_else(umi_retained_percent > 99.65, 1, 0)
  )

p_a <- ggplot(panel_a_data, aes(y = y)) +
  annotate(
    "rect",
    xmin = 95.55, xmax = 100.45, ymin = 6.72, ymax = 7.48,
    fill = pale_neutral, color = neutral_light, linewidth = 0.4
  ) +
  annotate(
    "text",
    x = 95.72, y = 7.10, hjust = 0,
    label = "One-to-one map: 20,684 cells  |  27,998 features  |  7,371 analysis cells",
    family = font_family, size = 2.25, color = neutral_dark
  ) +
  geom_segment(
    aes(x = 95.6, xend = umi_retained_percent, yend = y),
    color = neutral_light, linewidth = 0.65, lineend = "round"
  ) +
  geom_point(
    aes(x = umi_retained_percent, fill = source_group),
    shape = 21, size = 2.8, stroke = 0.65, color = neutral_dark
  ) +
  geom_text(
    aes(x = label_x, label = label, hjust = label_hjust),
    family = font_family, size = 2.2, color = neutral_dark
  ) +
  geom_vline(
    xintercept = 100,
    color = neutral_mid,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_fill_manual(values = group_colors, labels = group_labels, guide = "none") +
  scale_x_continuous(
    breaks = c(96, 97, 98, 99, 100),
    limits = c(95.5, 100.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = rev(seq_along(samples)),
    labels = samples,
    limits = c(0.5, 7.55),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "A",
    x = "UMIs retained in GSE289098 on the fixed cell set (%)",
    y = NULL
  ) +
  base_theme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    plot.margin = margin(11, 8, 4, 7)
  )

# Panel B: paired animal-level seven-transcript panel values.
identity_limits <- c(1.95, 4.85)
panel_b_data <- animal_panel |>
  mutate(
    label_x = recode(
      sample,
      C1 = 2.06, C2 = 2.00, C3 = 2.70,
      S1 = 3.27, S2 = 2.55, S3 = 4.42
    ),
    label_y = recode(
      sample,
      C1 = 2.35, C2 = 2.03, C3 = 2.62,
      S1 = 3.25, S2 = 2.25, S3 = 4.40
    )
  )
p_b <- ggplot(
  panel_b_data,
  aes(x = primary_value, y = alternative_value)
) +
  geom_abline(
    slope = 1, intercept = 0,
    color = neutral_mid, linetype = "dashed", linewidth = 0.45
  ) +
  geom_segment(
    aes(xend = label_x, yend = label_y),
    color = neutral_mid,
    linewidth = 0.35,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = source_group),
    shape = 21, size = 3.0, stroke = 0.7, color = neutral_dark
  ) +
  geom_text(
    aes(x = label_x, y = label_y, label = sample),
    family = font_family,
    size = 2.35,
    fontface = "bold",
    color = neutral_dark,
    hjust = 0.5,
    vjust = 0.5
  ) +
  scale_fill_manual(
    values = group_colors,
    labels = group_labels,
    name = NULL
  ) +
  scale_x_continuous(limits = identity_limits, breaks = 2:4) +
  scale_y_continuous(limits = identity_limits, breaks = 2:4) +
  coord_equal(clip = "off") +
  labs(
    tag = "B",
    x = "Seven-transcript mean, GSE267933 primary counts\n(log2 counts per million)",
    y = "Seven-transcript mean, GSE289098 counts\n(log2 counts per million)"
  ) +
  base_theme +
  theme(
    panel.grid.major = element_line(color = neutral_faint, linewidth = 0.35),
    legend.position = "bottom",
    legend.key.width = grid::unit(3.2, "mm"),
    legend.margin = margin(t = -2),
    plot.margin = margin(11, 7, 4, 7)
  )

# Panel C: animal-by-gene payload changes.
panel_c_data <- gene_delta |>
  mutate(
    sample = factor(sample, levels = samples),
    gene = factor(gene, levels = rev(genes)),
    display_group = factor(
      display_group,
      levels = c("Oxygen control", "Combined exposure")
    ),
    value_label = sprintf("%+.2f", alternative_minus_primary_log2_cpm),
    text_color = if_else(
      abs(alternative_minus_primary_log2_cpm) >= 0.13,
      "white",
      neutral_dark
    )
  )

p_c <- ggplot(
  panel_c_data,
  aes(x = sample, y = gene, fill = alternative_minus_primary_log2_cpm)
) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(
    aes(label = value_label, color = text_color),
    family = font_family,
    size = 2.05
  ) +
  facet_grid(
    cols = vars(display_group),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.22, 0.22),
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
    name = expression(Delta*" log2 counts per million")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(43, "mm"),
      barheight = grid::unit(2.2, "mm")
    )
  ) +
  scale_color_identity() +
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic(", x, ")")),
    expand = expansion(add = 0)
  ) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(tag = "C", x = NULL, y = NULL) +
  theme_minimal(base_family = font_family, base_size = 8.2) +
  theme(
    text = element_text(color = neutral_dark),
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7.1, face = "bold"),
    axis.text.y = element_text(size = 7.2),
    strip.text = element_text(size = 7.2, face = "bold"),
    strip.background = element_rect(
      fill = pale_neutral, color = neutral_light, linewidth = 0.4
    ),
    panel.spacing.x = grid::unit(2.3, "mm"),
    legend.position = "bottom",
    legend.key.width = grid::unit(19, "mm"),
    legend.key.height = grid::unit(2.2, "mm"),
    legend.title = element_text(size = 7.0),
    legend.text = element_text(size = 6.8),
    plot.margin = margin(5, 6, 4, 7)
  )

# Panel D: full and matched one-animal-omission panel contrasts.
panel_d_wide <- panel_effects |>
  mutate(
    scenario = factor(
      scenario,
      levels = rev(c("Full cohort", paste("Drop", samples)))
    ),
    delta_label = sprintf("%+.3f", alternative_minus_primary_difference)
  )
panel_d_long <- bind_rows(
  panel_d_wide |>
    transmute(
      scenario,
      scenario_order,
      payload = "GSE267933 primary counts",
      mean_difference = primary_mean_difference
    ),
  panel_d_wide |>
    transmute(
      scenario,
      scenario_order,
      payload = "GSE289098 aggregated counts",
      mean_difference = gse289098_mean_difference
    )
) |>
  mutate(
    payload = factor(
      payload,
      levels = c(
        "GSE267933 primary counts",
        "GSE289098 aggregated counts"
      )
    )
  )

p_d <- ggplot(panel_d_long, aes(x = mean_difference, y = scenario)) +
  geom_hline(
    yintercept = seq_len(7),
    color = neutral_faint,
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    color = neutral_mid,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_segment(
    data = panel_d_wide,
    aes(
      x = primary_mean_difference,
      xend = gse289098_mean_difference,
      y = scenario,
      yend = scenario
    ),
    inherit.aes = FALSE,
    color = neutral_light,
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = payload, shape = payload),
    size = 2.8,
    stroke = 0.75,
    color = primary_color
  ) +
  geom_text(
    data = panel_d_wide,
    aes(x = 1.65, y = scenario, label = delta_label),
    inherit.aes = FALSE,
    hjust = 0,
    family = font_family,
    size = 2.05,
    color = neutral_mid
  ) +
  annotate(
    "text", x = 1.82, y = 7.42, label = "GSE289098 − primary",
    hjust = 1, family = font_family, size = 2.1,
    fontface = "bold", color = neutral_mid
  ) +
  annotate(
    "point", x = 0.12, y = 7.42,
    shape = 21, size = 2.5, stroke = 0.7,
    fill = "white", color = primary_color
  ) +
  annotate(
    "text", x = 0.18, y = 7.42, label = "Primary",
    hjust = 0, family = font_family, size = 2.05,
    color = neutral_dark
  ) +
  annotate(
    "point", x = 0.58, y = 7.42,
    shape = 21, size = 2.5, stroke = 0.7,
    fill = alternative_color, color = primary_color
  ) +
  annotate(
    "text", x = 0.64, y = 7.42, label = "GSE289098",
    hjust = 0, family = font_family, size = 2.05,
    color = neutral_dark
  ) +
  scale_fill_manual(
    values = c(
      "GSE267933 primary counts" = "white",
      "GSE289098 aggregated counts" = alternative_color
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "GSE267933 primary counts" = 21,
      "GSE289098 aggregated counts" = 21
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1.0, 1.5),
    limits = c(-0.05, 1.85),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    labels = function(x) sub("^Drop ", "Omit ", x),
    expand = expansion(add = c(0.45, 0.70))
  ) +
  labs(
    tag = "D",
    x = paste0(
      "Difference in seven-transcript mean\n",
      "(combined exposure − oxygen control)\n",
      "(log2 counts per million)"
    ),
    y = NULL
  ) +
  base_theme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    legend.position = "none",
    plot.margin = margin(5, 7, 4, 7)
  )

top_row <- p_a + p_b + plot_layout(widths = c(1.12, 0.88))
bottom_row <- p_c + p_d + plot_layout(widths = c(1.06, 0.94))
figure_s6 <- top_row / bottom_row + plot_layout(heights = c(0.92, 1.08))

# Export and execution record.
base_name <- "FigureS3_same_cohort_count_payload_sensitivity"
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

save_quartz_pdf(figure_s6, pdf_path, width_in, height_in, font_family)

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
capture_render_warnings(print(figure_s6))
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
capture_render_warnings(print(figure_s6))
grDevices::dev.off()

saveRDS(figure_s6, rds_path, compress = "xz")
output_paths <- c(
  vector_pdf_master = pdf_path,
  glia_submission_tiff = tiff_path,
  png_preview = png_path,
  plot_object = rds_path
)

writeLines(
  if (length(render_warnings)) render_warnings else "None",
  file.path(manifest_dir, "figS03_render_warnings.txt")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "figS03_plot_sessionInfo.txt")
)
write_csv(
  tibble(
    input = names(input_paths),
    relative_path = file.path("data", "panel_ready", basename(input_paths)),
    sha256 = vapply(input_paths, sha256_file, character(1))
  ),
  file.path(manifest_dir, "figS03_plot_input_manifest.csv")
)
write_csv(
  tibble(
    output = names(output_paths),
    relative_path = file.path("outputs", basename(output_paths)),
    sha256 = vapply(output_paths, sha256_file, character(1))
  ),
  file.path(manifest_dir, "figS03_plot_output_manifest.csv")
)

execution <- list(
  figure = "Figure S3",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  dimensions = list(width_mm = width_mm, height_mm = height_mm),
  accession = "GSE289098",
  role = "same-cohort processed-count sensitivity only",
  biological_unit = "animal/library",
  fixed_cell_set = 7371L,
  panels = list(
    A = "Per-animal retained UMI percentage and exact mapping audit",
    B = "Paired animal-level seven-transcript panel values",
    C = "GSE289098-minus-primary gene-level log2 counts-per-million values",
    D = "Paired full and six leave-one-animal-out panel contrasts"
  ),
  selected_panel_status =
    "selected during exploratory inspection; secondary descriptive analysis",
  common_estimand_build_gate = "PASS",
  panel_full_and_loo_direction_gate = "PASS",
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
  file.path(manifest_dir, "figS03_plot_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Figure S3 rendering complete.")
message("PDF: ", pdf_path)
message("TIFF: ", tiff_path)
message("PNG: ", png_path)
