#!/usr/bin/env Rscript

# GLIA major revision — canonical Figure 4 plotting script
#
# Figure 4 places the primary single-cell findings in external cross-experiment
# context using GSE283401. It does not present GSE283401 as replication,
# validation, or a within-animal trajectory.
#
# Output: embedded-font Quartz PDF master, 300-dpi 8-bit LZW TIFF submission
#         file, PNG preview, retained R plot object, hashes, warnings, and
#         session record.
#
# Example:
# Rscript scripts/87_plot_glia_R1_figure4.R \
#   --figure-root /path/to/figure4

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(ggplot2)
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

assert_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  assert_true(
    length(missing_columns) == 0L,
    paste0(
      label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

near <- function(x, target, tolerance = 1e-10) {
  length(x) == 1L &&
    is.finite(x) &&
    abs(x - target) <= tolerance
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
height_mm <- as.numeric(get_arg("--height-mm", "165"))
dpi <- as.integer(get_arg("--dpi", "300"))
preview_suffix <- get_arg("--preview-suffix", "")
assert_true(is.finite(width_mm) && width_mm > 0, "Invalid --width-mm.")
assert_true(is.finite(height_mm) && height_mm > 0, "Invalid --height-mm.")
assert_true(
  width_mm == 180 && height_mm == 165,
  "Canonical Figure 4 dimensions are frozen at 180 x 165 mm."
)
assert_true(dpi == 300L, "GLIA revision submission TIFF is frozen at 300 dpi.")
assert_true(
  grepl("^(_[A-Za-z0-9][A-Za-z0-9._-]*)?$", preview_suffix),
  paste(
    "Invalid --preview-suffix; use an empty value or a leading underscore",
    "followed by letters, numbers, dots, underscores, or hyphens."
  )
)

input_paths <- c(
  metadata = file.path(panel_dir, "fig04_metadata.csv"),
  design_counts = file.path(panel_dir, "fig04_design_counts.csv"),
  study_comparison = file.path(panel_dir, "fig04_study_comparison.csv"),
  selected_gene_timepoint = file.path(
    panel_dir, "fig04_selected_gene_timepoint.csv"
  ),
  hallmark_interaction_focus = file.path(
    panel_dir, "fig04_hallmark_interaction_focus.csv"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing panel-ready inputs:", paste(basename(missing_inputs), collapse = ", "))
)

metadata <- read_csv(input_paths[["metadata"]], show_col_types = FALSE)
design_counts <- read_csv(
  input_paths[["design_counts"]],
  show_col_types = FALSE
)
study_comparison <- read_csv(
  input_paths[["study_comparison"]],
  show_col_types = FALSE
)
selected_gene_timepoint <- read_csv(
  input_paths[["selected_gene_timepoint"]],
  show_col_types = FALSE
)
hallmark_interaction_focus <- read_csv(
  input_paths[["hallmark_interaction_focus"]],
  show_col_types = FALSE
)

# ---- Frozen schema and numerical anchors ----------------------------------

assert_columns(
  metadata,
  c(
    "figure_id", "external_accession", "biological_unit", "old_model_n",
    "old_6h_control_n", "old_6h_exposure_n", "old_48h_control_n",
    "old_48h_exposure_n", "interaction", "database", "claim_ceiling",
    "limitation"
  ),
  "fig04_metadata.csv"
)
assert_columns(
  design_counts,
  c(
    "time", "time_order", "age", "age_order", "treatment",
    "treatment_order", "n_animals"
  ),
  "fig04_design_counts.csv"
)
assert_columns(
  study_comparison,
  c(
    "dataset", "role", "age", "assay", "exposure", "sampling",
    "biological_n", "animal_relation", "claim_role"
  ),
  "fig04_study_comparison.csv"
)
assert_columns(
  selected_gene_timepoint,
  c(
    "estimand", "symbol", "log2FoldChange", "ci95_low", "ci95_high",
    "padj", "transcript_order", "time", "time_order", "n_control",
    "n_combined_exposure", "panel_status"
  ),
  "fig04_selected_gene_timepoint.csv"
)
assert_columns(
  hallmark_interaction_focus,
  c(
    "estimand", "ID", "NES", "p.adjust", "fdr_lt_0_05",
    "hallmark_display", "hallmark_order", "interpretation", "limitation"
  ),
  "fig04_hallmark_interaction_focus.csv"
)

assert_true(
  nrow(metadata) == 1L &&
    metadata$figure_id[[1L]] == "Figure 4" &&
    metadata$external_accession[[1L]] == "GSE283401" &&
    metadata$biological_unit[[1L]] == "Animal/library" &&
    metadata$old_model_n[[1L]] == 32L,
  "Figure 4 metadata do not match the frozen external-context design."
)
assert_true(
  identical(
    as.integer(c(
      metadata$old_6h_control_n,
      metadata$old_6h_exposure_n,
      metadata$old_48h_control_n,
      metadata$old_48h_exposure_n
    )),
    c(8L, 7L, 8L, 9L)
  ),
  "Old-animal sample counts must remain 8/7 at 6 h and 8/9 at 48 h."
)
assert_true(
  nrow(design_counts) == 8L &&
    sum(design_counts$n_animals) == 48L &&
    sum(
      design_counts$n_animals[
        design_counts$age == "Old (20–22 months)"
      ]
    ) == 32L,
  "Figure 4A must retain the full 48-library design and 32-library old model."
)
assert_true(
  nrow(study_comparison) == 3L &&
    setequal(study_comparison$dataset, c("GSE267933", "GSE283401")),
  "Figure 4A study-comparison rows are incomplete."
)

gene_order <- c("Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3")
assert_true(
  nrow(selected_gene_timepoint) == 14L &&
    setequal(selected_gene_timepoint$symbol, gene_order) &&
    setequal(selected_gene_timepoint$time, c("6 h", "48 h")) &&
    nrow(distinct(selected_gene_timepoint, symbol, time)) == 14L &&
    all(is.finite(selected_gene_timepoint$log2FoldChange)) &&
    all(is.finite(selected_gene_timepoint$ci95_low)) &&
    all(is.finite(selected_gene_timepoint$ci95_high)),
  "Figure 4B must contain seven unique transcripts at each of two timepoints."
)

gene_anchor <- function(symbol, time) {
  selected_gene_timepoint |>
    filter(.data$symbol == .env$symbol, .data$time == .env$time) |>
    pull(log2FoldChange)
}
assert_true(
  near(gene_anchor("Irf7", "6 h"), -0.5030606386528876) &&
    near(gene_anchor("Irf7", "48 h"), 0.4220254975908904) &&
    near(gene_anchor("Ifit2", "6 h"), -0.30933761811634319) &&
    near(gene_anchor("Ifit2", "48 h"), 1.3515671136934) &&
    near(gene_anchor("Ifit3", "6 h"), -1.2460914380118211) &&
    near(gene_anchor("Ifit3", "48 h"), 0.7602468852222085),
  "Figure 4B gene-level numerical anchors are not verified."
)

hallmark_ids <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
assert_true(
  nrow(hallmark_interaction_focus) == 4L &&
    setequal(hallmark_interaction_focus$ID, hallmark_ids) &&
    all(hallmark_interaction_focus$estimand == "time_by_treatment") &&
    nrow(distinct(hallmark_interaction_focus, ID)) == 4L,
  "Figure 4C must contain four unique focused Hallmark interactions."
)

hallmark_anchor <- function(id, column) {
  hallmark_interaction_focus |>
    filter(.data$ID == .env$id) |>
    pull(all_of(column))
}
assert_true(
  near(
    hallmark_anchor(hallmark_ids[[1L]], "NES"),
    2.624158343520711
  ) &&
    near(
      hallmark_anchor(hallmark_ids[[1L]], "p.adjust"),
      7.313829969354819e-11,
      1e-16
    ) &&
    near(
      hallmark_anchor(hallmark_ids[[2L]], "NES"),
      2.365226167247039
    ) &&
    near(
      hallmark_anchor(hallmark_ids[[3L]], "NES"),
      1.3432884692591072
    ) &&
    near(
      hallmark_anchor(hallmark_ids[[3L]], "p.adjust"),
      0.11355003162542474
    ) &&
    near(
      hallmark_anchor(hallmark_ids[[4L]], "NES"),
      1.708914471375484
    ),
  "Figure 4C Hallmark interaction anchors are not verified."
)
assert_true(
  identical(
    hallmark_interaction_focus$fdr_lt_0_05[
      order(hallmark_interaction_focus$hallmark_order)
    ],
    c(TRUE, TRUE, FALSE, TRUE)
  ),
  "Figure 4C focused Hallmark BH-FDR status is not verified."
)

# ---- Shared visual specification ------------------------------------------

font_match <- systemfonts::match_fonts("Arial")
font_family <- if (nrow(font_match) > 0L && nzchar(font_match$path[[1L]])) {
  "Arial"
} else {
  "Helvetica"
}
if (font_family != "Arial") {
  warning("Arial was not found; using Helvetica.")
}

control_color <- "#0072B2"
exposure_color <- "#D55E00"
interaction_color <- "#6A51A3"
pale_blue <- "#E7F2F8"
pale_orange <- "#FBEDE5"
pale_purple <- "#F2EDF8"
neutral_dark <- "#333333"
neutral_mid <- "#707070"
neutral_light <- "#D2D2D2"
neutral_faint <- "#ECECEC"

base_theme <- theme_classic(base_family = font_family, base_size = 8.2) +
  theme(
    text = element_text(color = neutral_dark),
    plot.title = element_text(
      size = 9, face = "bold", hjust = 0,
      margin = margin(b = 1.5)
    ),
    plot.subtitle = element_text(
      size = 7.2, color = neutral_mid, hjust = 0,
      lineheight = 0.96, margin = margin(b = 3)
    ),
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 7.8),
    axis.text = element_text(size = 7.2, color = neutral_dark),
    strip.text = element_text(size = 7.4, face = "bold"),
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 7.0),
    plot.margin = margin(5, 5, 4, 7)
  )

# ---- Panel A: independent dataset and design audit ------------------------

panel_a_cards <- data.frame(
  xmin = c(2, 51),
  xmax = c(49, 98),
  ymin = c(0.8, 0.8),
  ymax = c(19.2, 19.2),
  header_fill = c(pale_blue, pale_purple),
  header_border = c(control_color, interaction_color)
)

panel_a_row_lines <- data.frame(
  x = rep(c(2, 51), each = 3),
  xend = rep(c(49, 98), each = 3),
  y = rep(c(13.8, 10.6, 6.8), times = 2)
)

panel_a_row_labels <- data.frame(
  x = rep(c(4, 53), each = 4),
  y = rep(c(15.1, 12.1, 8.8, 4.7), times = 2),
  label = rep(c("Age", "Assay", "Exposure", "Samples"), times = 2)
)

panel_a_values <- data.frame(
  x = c(13, 13, 13, 62, 62, 62),
  y = c(15.1, 12.1, 8.8, 15.1, 12.1, 8.8),
  label = c(
    "18-month male mice",
    "10x scRNA-seq (hippocampus)",
    "2.5% sevoflurane, 30 min + laparotomy",
    "20–22-month male mice (old subset)",
    "FACS microglial bulk RNA-seq (hippocampus)",
    "1.2% isoflurane in 30% O2, 2 h + laparotomy"
  )
)

panel_a_sample_values <- data.frame(
  x = c(13, 19, 36, 62, 73, 87, 62, 73, 87),
  y = c(4.7, 4.7, 4.7, 5.7, 5.7, 5.7, 3.6, 3.6, 3.6),
  label = c(
    "24 h", "Oxygen control n=3", "Exposed n=3",
    "6 h from start", "Carrier gas n=8", "Exposed n=7",
    "48 h from start", "Carrier gas n=8", "Exposed n=9"
  ),
  color = c(
    neutral_dark, control_color, exposure_color,
    neutral_dark, control_color, exposure_color,
    neutral_dark, control_color, exposure_color
  ),
  face = c("bold", rep("plain", 2), "bold", rep("plain", 2), "bold", rep("plain", 2))
)

p_a <- ggplot() +
  geom_rect(
    data = panel_a_cards,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    color = neutral_light,
    linewidth = 0.5
  ) +
  geom_rect(
    data = panel_a_cards,
    aes(xmin = xmin, xmax = xmax, ymin = 16.3, ymax = ymax),
    fill = panel_a_cards$header_fill,
    color = panel_a_cards$header_border,
    linewidth = 0.55
  ) +
  geom_segment(
    data = panel_a_row_lines,
    aes(x = x, xend = xend, y = y, yend = y),
    color = neutral_faint,
    linewidth = 0.4
  ) +
  annotate(
    "text", x = 4, y = 17.75, hjust = 0,
    label = "GSE267933  |  Primary cohort",
    family = font_family, fontface = "bold", size = 3.0,
    color = control_color
  ) +
  annotate(
    "text", x = 53, y = 17.75, hjust = 0,
    label = "GSE283401  |  External context",
    family = font_family, fontface = "bold", size = 3.0,
    color = interaction_color
  ) +
  geom_text(
    data = panel_a_row_labels,
    aes(x = x, y = y, label = label),
    family = font_family, fontface = "bold", size = 2.25,
    hjust = 0, color = neutral_mid
  ) +
  geom_text(
    data = panel_a_values,
    aes(x = x, y = y, label = label),
    family = font_family, size = 2.35,
    hjust = 0, color = neutral_dark
  ) +
  geom_text(
    data = panel_a_sample_values,
    aes(x = x, y = y, label = label),
    family = font_family, fontface = panel_a_sample_values$face,
    size = 2.25, hjust = 0,
    color = panel_a_sample_values$color
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 20), clip = "off") +
  labs(tag = "A") +
  theme_void(base_family = font_family, base_size = 8.2) +
  theme(
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(4, 5, 1, 7)
  )

# ---- Panel B: selected transcript contrasts at each deposited time --------

time_labels <- c(
  "6 h" = "6 h\n8 control; 7 exposed",
  "48 h" = "48 h\n8 control; 9 exposed"
)

selected_gene_plot <- selected_gene_timepoint |>
  mutate(
    symbol = factor(symbol, levels = rev(gene_order)),
    time_label = factor(
      unname(time_labels[time]),
      levels = unname(time_labels[c("6 h", "48 h")])
    ),
    favored_group = if_else(
      log2FoldChange >= 0,
      "Combined exposure",
      "Carrier-gas control"
    )
  )

p_b <- ggplot(
  selected_gene_plot,
  aes(x = log2FoldChange, y = symbol)
) +
  geom_hline(
    yintercept = seq_along(gene_order),
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
    aes(x = ci95_low, xend = ci95_high, yend = symbol),
    color = neutral_dark,
    linewidth = 0.65,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = favored_group),
    shape = 21,
    size = 2.55,
    stroke = 0.65,
    color = neutral_dark
  ) +
  facet_grid(cols = vars(time_label)) +
  scale_fill_manual(
    values = c(
      "Carrier-gas control" = control_color,
      "Combined exposure" = exposure_color
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = seq(-3, 3, by = 1),
    limits = c(-3.5, 3.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    labels = function(x) parse(text = paste0("italic(", x, ")")),
    expand = expansion(add = 0.45)
  ) +
  labs(
    tag = "B",
    x = paste0(
      "DESeq2 log2 fold change (exposure \u2212 carrier-gas control)\n",
      "\u2190 more positive in carrier-gas control          ",
      "more positive in exposure \u2192"
    ),
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "#F5F5F5", color = neutral_light, linewidth = 0.4
    ),
    panel.spacing.x = grid::unit(4.5, "mm"),
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.4),
    plot.margin = margin(5, 7, 5, 7)
  )

# ---- Panel C: formal Hallmark time-by-treatment interaction ---------------

format_fdr <- function(x) {
  ifelse(
    x < 1e-3,
    gsub("-", "\u2212", formatC(x, format = "e", digits = 2), fixed = TRUE),
    ifelse(x < 0.01, sprintf("%.4f", x), sprintf("%.3f", x))
  )
}

hallmark_plot <- hallmark_interaction_focus |>
  arrange(hallmark_order) |>
  mutate(
    hallmark_label = factor(
      c(
        "Interferon alpha\nresponse",
        "Interferon gamma\nresponse",
        "TNFA signaling\nvia NF-kB",
        "Hallmark inflammatory-\nresponse gene set"
      ),
      levels = rev(c(
        "Interferon alpha\nresponse",
        "Interferon gamma\nresponse",
        "TNFA signaling\nvia NF-kB",
        "Hallmark inflammatory-\nresponse gene set"
      ))
    ),
    significance = if_else(
      fdr_lt_0_05,
      "BH FDR < 0.05",
      "BH FDR \u2265 0.05"
    ),
    direct_label = paste0(
      "NES ", sprintf("%.2f", NES), "\nFDR ", format_fdr(p.adjust)
    )
  )

p_c <- ggplot(hallmark_plot, aes(x = NES, y = hallmark_label)) +
  geom_hline(
    yintercept = seq_len(4),
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
    aes(x = 0, xend = NES, yend = hallmark_label),
    color = interaction_color,
    linewidth = 0.8,
    lineend = "round"
  ) +
  geom_point(
    aes(fill = significance),
    shape = 21,
    size = 3.0,
    stroke = 0.75,
    color = interaction_color
  ) +
  geom_text(
    aes(label = direct_label),
    family = font_family,
    size = 2.25,
    hjust = 0,
    x = -2.85,
    lineheight = 0.9,
    color = neutral_dark
  ) +
  scale_fill_manual(
    values = c(
      "BH FDR < 0.05" = interaction_color,
      "BH FDR \u2265 0.05" = "white"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = -3:3,
    limits = c(-3.0, 3.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = 0.55)) +
  labs(
    tag = "C",
    x = paste0(
      "NES from time × treatment-ranked genes\n",
      "(48-h exposure contrast) \u2212 (6-h exposure contrast)\n",
      "\u2190 genes rank more positive at 6 h       at 48 h \u2192"
    ),
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.8, lineheight = 0.92),
    plot.margin = margin(5, 7, 5, 7)
  )

# ---- Assembly --------------------------------------------------------------

bottom_row <- p_b + p_c + plot_layout(widths = c(1.16, 0.84))
figure_4 <- p_a / bottom_row +
  plot_layout(heights = c(0.52, 1.48))

# ---- Export and execution record ------------------------------------------

base_name <- "Figure4_independent_temporal_context"
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

save_quartz_pdf(figure_4, pdf_path, width_in, height_in, font_family)

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
capture_render_warnings(print(figure_4))
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
capture_render_warnings(print(figure_4))
grDevices::dev.off()

saveRDS(figure_4, rds_path, compress = "xz")

output_paths <- c(
  vector_pdf_master = pdf_path,
  glia_submission_tiff = tiff_path,
  png_preview = png_path,
  plot_object = rds_path
)

writeLines(
  if (length(render_warnings)) render_warnings else "None",
  file.path(manifest_dir, "fig04_render_warnings.txt")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "fig04_plot_sessionInfo.txt")
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
write_csv(input_manifest, file.path(manifest_dir, "fig04_plot_input_manifest.csv"))
write_csv(output_manifest, file.path(manifest_dir, "fig04_plot_output_manifest.csv"))

execution_manifest <- list(
  figure = "Figure 4",
  purpose = paste(
    "Independent temporal context for selected transcript contrasts and",
    "focused Hallmark time-by-treatment interactions"
  ),
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  dimensions = list(width_mm = width_mm, height_mm = height_mm),
  preview_suffix = preview_suffix,
  external_accession = "GSE283401",
  external_role = "independent context; not replication or validation",
  biological_unit = "animal/library",
  old_model_n = 32L,
  old_design = list(
    six_h_control_n = 8L,
    six_h_exposure_n = 7L,
    forty_eight_h_control_n = 8L,
    forty_eight_h_exposure_n = 9L
  ),
  gene_display = list(
    selection_status =
      "selected during exploratory inspection; secondary descriptive analysis",
    estimator = "unshrunk DESeq2 log2 fold change with 95% Wald CI",
    fits = "separate old-animal treatment models at 6 h and 48 h"
  ),
  hallmark_display = list(
    estimator = paste(
      "mouse MSigDB 2026.1.Mm Hallmark NES ranked by signed",
      "DESeq2 time-by-treatment Wald statistic"
    ),
    interaction = paste(
      "(48-h exposure minus control) minus",
      "(6-h exposure minus control)"
    ),
    multiplicity = "BH FDR across all 50 Hallmark gene sets"
  ),
  claim_ceiling = metadata$claim_ceiling[[1L]],
  limitation = metadata$limitation[[1L]],
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
  random_seed = NULL,
  render_warnings = as.list(render_warnings),
  r_version = R.version.string,
  packages = as.list(c(
    digest = as.character(packageVersion("digest")),
    dplyr = as.character(packageVersion("dplyr")),
    ggplot2 = as.character(packageVersion("ggplot2")),
    jsonlite = as.character(packageVersion("jsonlite")),
    patchwork = as.character(packageVersion("patchwork")),
    ragg = as.character(packageVersion("ragg")),
    readr = as.character(packageVersion("readr")),
    systemfonts = as.character(packageVersion("systemfonts"))
  )),
  outputs = split(output_manifest, seq_len(nrow(output_manifest)))
)
write_json(
  execution_manifest,
  file.path(manifest_dir, "fig04_plot_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Figure 4 written to: ", output_dir)
message("PDF master: ", basename(pdf_path))
message("TIFF submission file: ", basename(tiff_path), " (", dpi, " dpi, LZW)")
message("PNG preview: ", basename(png_path))
