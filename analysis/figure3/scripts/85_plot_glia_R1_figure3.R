#!/usr/bin/env Rscript

# GLIA major revision — canonical Figure 3 plotting script
#
# Figure 3 is an influence and estimator-sensitivity display, not a validation
# figure. All panels, assembly, and exports are generated in R from the
# panel-ready CSV files frozen by script 84.
#
# Output: embedded-font Quartz PDF master, 300-dpi 8-bit LZW TIFF submission
#         file, PNG preview, retained R plot object, hashes, warnings, and
#         session record.
#
# Example:
# Rscript scripts/85_plot_glia_R1_figure3.R \
#   --figure-root /path/to/figure3

suppressPackageStartupMessages({
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
height_mm <- as.numeric(get_arg("--height-mm", "155"))
dpi <- as.integer(get_arg("--dpi", "300"))
preview_suffix <- get_arg("--preview-suffix", "")
assert_true(is.finite(width_mm) && width_mm > 0, "Invalid --width-mm.")
assert_true(is.finite(height_mm) && height_mm > 0, "Invalid --height-mm.")
assert_true(
  width_mm == 180 && height_mm == 155,
  "Canonical Figure 3 dimensions are frozen at 180 x 155 mm."
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
  metadata = file.path(panel_dir, "fig03_metadata.csv"),
  sample_metadata = file.path(panel_dir, "fig03_sample_metadata.csv"),
  isg_score_animal_values = file.path(
    panel_dir, "fig03_isg_score_animal_values.csv"
  ),
  isg_score_effect_loo = file.path(
    panel_dir, "fig03_isg_score_effect_loo.csv"
  ),
  hallmark_focus_loo = file.path(
    panel_dir, "fig03_hallmark_focus_loo.csv"
  ),
  inclusion_panel_full = file.path(
    panel_dir, "fig03_inclusion_panel_full.csv"
  ),
  inclusion_hallmark_full = file.path(
    panel_dir, "fig03_inclusion_hallmark_full.csv"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing panel-ready inputs:", paste(basename(missing_inputs), collapse = ", "))
)

metadata <- read_csv(input_paths[["metadata"]], show_col_types = FALSE)
sample_metadata <- read_csv(
  input_paths[["sample_metadata"]],
  show_col_types = FALSE
)
isg_score_animal_values <- read_csv(
  input_paths[["isg_score_animal_values"]],
  show_col_types = FALSE
)
isg_score_effect_loo <- read_csv(
  input_paths[["isg_score_effect_loo"]],
  show_col_types = FALSE
)
hallmark_focus_loo <- read_csv(
  input_paths[["hallmark_focus_loo"]],
  show_col_types = FALSE
)
inclusion_panel_full <- read_csv(
  input_paths[["inclusion_panel_full"]],
  show_col_types = FALSE
)
inclusion_hallmark_full <- read_csv(
  input_paths[["inclusion_hallmark_full"]],
  show_col_types = FALSE
)

samples <- c("C1", "C2", "C3", "S1", "S2", "S3")
group_levels <- c("Oxygen control", "Combined exposure")
scenario_levels <- c(
  "Full cohort",
  "Omit C1", "Omit C2", "Omit C3",
  "Omit S1", "Omit S2", "Omit S3"
)
hallmark_ids <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
hallmark_labels <- c(
  "Interferon alpha\nresponse",
  "Interferon gamma\nresponse",
  "TNFA signaling\nvia NF-kB",
  "Hallmark inflammatory-\nresponse gene set"
)
names(hallmark_labels) <- hallmark_ids
analysis_set_levels <- c(
  "Primary: 7,371 Scrublet-negative microglia",
  "Sensitivity: all 7,461 frozen microglia"
)

# ---- Frozen schema and numerical anchors ----------------------------------

assert_columns(
  metadata,
  c(
    "figure_id", "primary_accession", "analysis_set", "biological_unit",
    "contrast", "score_estimator", "hallmark_estimator", "scale_rule"
  ),
  "fig03_metadata.csv"
)
assert_columns(
  sample_metadata,
  c(
    "sample", "sample_order", "source_group", "display_group",
    "condition_order", "biological_unit"
  ),
  "fig03_sample_metadata.csv"
)
animal_required <- c(
  "sample", "sample_order", "source_group", "display_group",
  "condition_order", "outcome", "estimator", "value"
)
effect_required <- c(
  "dropped_animal", "dropped_source_group", "mean_difference",
  "mean_difference_ci95_low", "mean_difference_ci95_high", "cohens_d",
  "scenario", "scenario_order", "is_full_cohort"
)
assert_columns(
  isg_score_animal_values,
  animal_required,
  "fig03_isg_score_animal_values.csv"
)
assert_columns(
  isg_score_effect_loo,
  effect_required,
  "fig03_isg_score_effect_loo.csv"
)
assert_columns(
  hallmark_focus_loo,
  c(
    "dropped_animal", "ID", "NES", "p.adjust", "hallmark_order",
    "hallmark_display", "scenario_order", "scenario", "is_full_cohort",
    "fdr_lt_0_05"
  ),
  "fig03_hallmark_focus_loo.csv"
)
assert_columns(
  inclusion_panel_full,
  c(
    "analysis_set", "analysis_set_order", "outcome", "estimator",
    "value", "delta_scrublet_negative_minus_allcell"
  ),
  "fig03_inclusion_panel_full.csv"
)
assert_columns(
  inclusion_hallmark_full,
  c(
    "ID", "hallmark_order", "analysis_set", "analysis_set_order",
    "estimator", "value", "adjusted_p_value",
    "delta_scrublet_negative_minus_allcell"
  ),
  "fig03_inclusion_hallmark_full.csv"
)

assert_true(
  nrow(metadata) == 1L &&
    metadata$figure_id[[1L]] == "Figure 3" &&
    metadata$primary_accession[[1L]] == "GSE267933" &&
    grepl("7,371", metadata$analysis_set[[1L]], fixed = TRUE) &&
    grepl("separate axes", metadata$scale_rule[[1L]], fixed = TRUE),
  "Figure 3 metadata differ from the frozen analysis and scale rules."
)
assert_true(
  nrow(sample_metadata) == 6L &&
    identical(sample_metadata$sample, samples) &&
    identical(as.integer(sample_metadata$sample_order), seq_along(samples)) &&
    identical(
      sample_metadata$display_group,
      c(rep("Oxygen control", 3), rep("Combined exposure", 3))
    ),
  "Figure 3 sample metadata differ from the six-animal design."
)

for (animal_data in list(isg_score_animal_values)) {
  assert_true(
    nrow(animal_data) == 6L &&
      identical(animal_data$sample, samples) &&
      !anyDuplicated(animal_data$sample) &&
      all(is.finite(animal_data$value)),
    "Each Figure 3 animal-value panel must contain six unique finite values."
  )
}

for (effect_data in list(isg_score_effect_loo)) {
  assert_true(
    nrow(effect_data) == 7L &&
      identical(effect_data$scenario, scenario_levels) &&
      identical(
        as.integer(effect_data$scenario_order),
        seq_along(scenario_levels)
      ) &&
      sum(effect_data$is_full_cohort) == 1L &&
      !anyDuplicated(effect_data$scenario),
    paste(
      "Each Figure 3 effect diagnostic must contain the full cohort and",
      "six ordered one-animal omissions."
    )
  )
}

isg_full <- isg_score_effect_loo |>
  filter(is_full_cohort)
isg_omit_s3 <- isg_score_effect_loo |>
  filter(dropped_animal == "S3")
assert_true(
  nrow(isg_full) == 1L &&
    near(isg_full$mean_difference, 0.04649471117813341) &&
    near(isg_full$mean_difference_ci95_low, -0.17662484192462355) &&
    near(isg_full$mean_difference_ci95_high, 0.26961426428089036) &&
    near(isg_full$cohens_d, 0.7137166150596256),
  "Figure 3A full-cohort score_genes anchors are not verified."
)
assert_true(
  nrow(isg_omit_s3) == 1L &&
    near(isg_omit_s3$mean_difference, -0.00608966251040171) &&
    isg_omit_s3$mean_difference < 0 &&
    all(
      isg_score_effect_loo$mean_difference[
        isg_score_effect_loo$dropped_animal != "S3"
      ] > 0
    ),
  "Figure 3A must retain the verified omit-S3 sign change."
)
assert_true(
  near(
    isg_score_animal_values$value[
      isg_score_animal_values$sample == "S3"
    ],
    0.10990842655042898
  ),
  "Figure 3A S3 animal score differs from the verified anchor."
)

assert_true(
  nrow(hallmark_focus_loo) == 28L &&
    setequal(hallmark_focus_loo$ID, hallmark_ids) &&
    setequal(hallmark_focus_loo$scenario, scenario_levels) &&
    nrow(distinct(hallmark_focus_loo, ID, scenario)) == 28L &&
    all(is.finite(hallmark_focus_loo$NES)) &&
    all(hallmark_focus_loo$NES >= -3 & hallmark_focus_loo$NES <= 3),
  "Figure 3B must contain four Hallmarks x seven unique fitted scenarios."
)

hallmark_anchor <- function(id, scenario) {
  hallmark_focus_loo |>
    filter(.data$ID == .env$id, .data$scenario == .env$scenario)
}

ifna_full <- hallmark_anchor(hallmark_ids[[1L]], "Full cohort")
ifng_full <- hallmark_anchor(hallmark_ids[[2L]], "Full cohort")
tnfa_full <- hallmark_anchor(hallmark_ids[[3L]], "Full cohort")
inflam_full <- hallmark_anchor(hallmark_ids[[4L]], "Full cohort")
ifna_omit_s3 <- hallmark_anchor(hallmark_ids[[1L]], "Omit S3")
ifng_omit_s3 <- hallmark_anchor(hallmark_ids[[2L]], "Omit S3")
assert_true(
  nrow(ifna_full) == 1L &&
    near(ifna_full$NES, 2.1714562753577713) &&
    near(ifna_full$p.adjust, 3.159876028891314e-6, 1e-12) &&
    nrow(ifng_full) == 1L &&
    near(ifng_full$NES, 1.5803088067385394) &&
    near(ifng_full$p.adjust, 0.00468642985416463, 1e-12) &&
    nrow(tnfa_full) == 1L &&
    near(tnfa_full$NES, -2.470813545455537) &&
    nrow(inflam_full) == 1L &&
    near(inflam_full$NES, -1.8461972784026033),
  "Figure 3B full-cohort Hallmark NES/FDR anchors are not verified."
)
assert_true(
  nrow(ifna_omit_s3) == 1L &&
    near(ifna_omit_s3$NES, -1.098079647814092) &&
    near(ifna_omit_s3$p.adjust, 0.343873813261867, 1e-10) &&
    nrow(ifng_omit_s3) == 1L &&
    near(ifng_omit_s3$NES, -1.568752656952539) &&
    near(ifng_omit_s3$p.adjust, 0.00334246714259, 1e-10),
  "Figure 3B omit-S3 interferon Hallmark anchors are not verified."
)
assert_true(
  all(
    hallmark_focus_loo$NES[
      hallmark_focus_loo$ID %in% hallmark_ids[3:4]
    ] < 0
  ),
  "TNFA/NF-kB and inflammatory-response Hallmark NES must remain negative."
)

assert_true(
  nrow(inclusion_panel_full) == 2L &&
    identical(inclusion_panel_full$analysis_set, analysis_set_levels) &&
    near(
      inclusion_panel_full$value[
        inclusion_panel_full$analysis_set == analysis_set_levels[[1L]]
      ],
      1.0685736753232482
    ) &&
    near(
      inclusion_panel_full$value[
        inclusion_panel_full$analysis_set == analysis_set_levels[[2L]]
      ],
      1.066003929153014
    ),
  "Figure 3C selected-panel inclusion anchors are not verified."
)
assert_true(
  nrow(inclusion_hallmark_full) == 8L &&
    setequal(inclusion_hallmark_full$ID, hallmark_ids) &&
    nrow(
      distinct(inclusion_hallmark_full, ID, analysis_set)
    ) == 8L,
  "Figure 3C must contain four Hallmarks x two cell-inclusion rules."
)
inclusion_hallmark_anchors <- c(
  2.1714562753577713, 2.2125531313154747,
  1.5803088067385394, 1.5666702727671056,
  -2.470813545455537, -2.496471789740173,
  -1.8461972784026033, -1.7866049673279405
)
inclusion_hallmark_ordered <- inclusion_hallmark_full |>
  arrange(hallmark_order, analysis_set_order)
assert_true(
  max(
    abs(inclusion_hallmark_ordered$value - inclusion_hallmark_anchors)
  ) < 1e-10,
  "Figure 3C Hallmark inclusion anchors are not verified."
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

condition_colors <- c(
  "Oxygen control" = "#0072B2",
  "Combined exposure" = "#D55E00"
)
condition_shapes <- c(
  "Oxygen control" = 16,
  "Combined exposure" = 17
)
source_group_colors <- c(
  "Control" = "#0072B2",
  "Surgery" = "#D55E00"
)
analysis_set_shapes <- c(
  "Primary: 7,371 Scrublet-negative microglia" = 16,
  "Sensitivity: all 7,461 frozen microglia" = 1
)
neutral_dark <- "#333333"
neutral_mid <- "#707070"
neutral_light <- "#D2D2D2"
neutral_faint <- "#ECECEC"

base_theme <- theme_classic(base_family = font_family, base_size = 8.2) +
  theme(
    text = element_text(color = neutral_dark),
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7.4, color = neutral_dark),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.title = element_text(size = 7.3),
    legend.text = element_text(size = 7.1),
    legend.key.height = grid::unit(3.0, "mm"),
    legend.key.width = grid::unit(4.6, "mm"),
    legend.spacing.x = grid::unit(1.2, "mm"),
    legend.margin = margin(0, 0, 1, 0),
    plot.margin = margin(5, 5, 4, 7)
  )

sample_x_offsets <- c(
  C1 = -0.13, C2 = 0, C3 = 0.13,
  S1 = -0.13, S2 = 0, S3 = 0.13
)

prepare_animal_values <- function(data) {
  data |>
    mutate(
      display_group = factor(display_group, levels = group_levels),
      group_x = condition_order,
      x_plot = group_x + unname(sample_x_offsets[sample])
    )
}

group_means <- function(data) {
  data |>
    summarise(
      group_mean = mean(value),
      group_x = first(group_x),
      .by = display_group
    )
}

# ---- Panel A: fixed score_genes and animal influence -----------------------

isg_animal_plot <- prepare_animal_values(isg_score_animal_values) |>
  mutate(
    label_x = x_plot + case_when(
      sample %in% c("C1", "S1") ~ -0.025,
      sample %in% c("C3", "S2") ~ 0.025,
      TRUE ~ 0
    ),
    label_y = value + case_when(
      sample %in% c("C1", "S2") ~ -0.010,
      sample == "S3" ~ 0.012,
      TRUE ~ 0.010
    )
  )
isg_group_means <- group_means(isg_animal_plot)

p_a_animal <- ggplot(
  isg_animal_plot,
  aes(
    x = x_plot,
    y = value,
    color = display_group,
    shape = display_group
  )
) +
  geom_hline(
    yintercept = 0,
    color = neutral_mid,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_segment(
    data = isg_group_means,
    aes(
      x = group_x - 0.18,
      xend = group_x + 0.18,
      y = group_mean,
      yend = group_mean
    ),
    inherit.aes = FALSE,
    color = neutral_dark,
    linewidth = 0.75,
    lineend = "round"
  ) +
  geom_point(size = 2.5, stroke = 0.45) +
  geom_text(
    aes(x = label_x, y = label_y, label = sample),
    family = font_family,
    size = 2.25,
    color = neutral_dark,
    show.legend = FALSE
  ) +
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_shape_manual(values = condition_shapes, guide = "none") +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Oxygen\ncontrol", "Combined\nexposure"),
    limits = c(0.63, 2.37),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(-0.05, 0.10, by = 0.05),
    labels = label_number(accuracy = 0.01),
    limits = c(-0.075, 0.130),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "A",
    title = "Animal-level scores",
    x = NULL,
    y = "Background-adjusted seven-transcript score\n(score units)"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = neutral_faint, linewidth = 0.35),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
    plot.tag.position = c(-0.015, 1.08),
    plot.margin = margin(9, 5, 1, 10)
  )

isg_loo_plot <- isg_score_effect_loo |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    y_plot = 8 - scenario_order,
    value_label = sprintf("%+.3f", mean_difference)
  )

p_a_loo <- ggplot(
  isg_loo_plot,
  aes(x = mean_difference, y = y_plot)
) +
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
  geom_vline(
    xintercept = isg_full$mean_difference,
    color = neutral_light,
    linetype = "dotted",
    linewidth = 0.55
  ) +
  geom_segment(
    data = isg_loo_plot |> filter(is_full_cohort),
    aes(
      x = mean_difference_ci95_low,
      xend = mean_difference_ci95_high,
      y = y_plot,
      yend = y_plot
    ),
    inherit.aes = FALSE,
    color = neutral_dark,
    linewidth = 0.70,
    lineend = "round"
  ) +
  geom_point(
    data = isg_loo_plot |> filter(is_full_cohort),
    shape = 18,
    size = 3.0,
    color = neutral_dark
  ) +
  geom_point(
    data = isg_loo_plot |> filter(!is_full_cohort),
    aes(color = dropped_source_group),
    shape = 21,
    size = 2.35,
    stroke = 0.65,
    fill = "white"
  ) +
  geom_text(
    data = isg_loo_plot |> filter(!is_full_cohort),
    aes(label = value_label),
    family = font_family,
    size = 2.25,
    color = neutral_dark,
    hjust = 0,
    nudge_x = 0.009
  ) +
  geom_text(
    data = isg_loo_plot |> filter(is_full_cohort),
    aes(label = value_label),
    family = font_family,
    size = 2.25,
    color = neutral_dark,
    hjust = 0,
    nudge_x = 0.009,
    nudge_y = 0.26
  ) +
  scale_color_manual(values = source_group_colors, guide = "none") +
  scale_x_continuous(
    breaks = seq(-0.2, 0.3, by = 0.1),
    labels = label_number(accuracy = 0.1),
    limits = c(-0.20, 0.30),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = 7:1,
    labels = scenario_levels,
    limits = c(0.55, 7.45),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Full cohort and one-animal-omission estimates",
    x = "Combined exposure \u2212 oxygen control (score units)",
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
    plot.margin = margin(2, 5, 4, 7)
  )

# ---- Panel B: focused Hallmark full and omission refits --------------------

hallmark_plot <- hallmark_focus_loo |>
  mutate(
    scenario = factor(scenario, levels = scenario_levels),
    hallmark_label = factor(
      unname(hallmark_labels[ID]),
      levels = rev(unname(hallmark_labels))
    ),
    nes_label = sprintf("%.2f", NES),
    text_color = if_else(abs(NES) >= 1.75, "white", neutral_dark)
  )

p_b <- ggplot(
  hallmark_plot,
  aes(x = scenario, y = hallmark_label, fill = NES)
) +
  geom_tile(color = "white", linewidth = 0.65) +
  geom_text(
    aes(label = nes_label, color = text_color),
    family = font_family,
    size = 2.30,
    fontface = "plain"
  ) +
  scale_color_identity(guide = "none") +
  scale_fill_gradient2(
    low = condition_colors[["Oxygen control"]],
    mid = "white",
    high = condition_colors[["Combined exposure"]],
    midpoint = 0,
    limits = c(-3, 3),
    oob = squish,
    breaks = c(-3, 0, 3),
    name = "NES: oxygen control \u2190 0 \u2192 combined exposure",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      barwidth = grid::unit(43, "mm"),
      barheight = grid::unit(2.7, "mm"),
      ticks = TRUE
    )
  ) +
  scale_x_discrete(
    labels = c("Full", "\u2212C1", "\u2212C2", "\u2212C3", "\u2212S1", "\u2212S2", "\u2212S3"),
    expand = expansion(add = 0)
  ) +
  scale_y_discrete(expand = expansion(add = 0)) +
  labs(
    tag = "B",
    x = "Full cohort and one-animal-omission refits",
    y = NULL
  ) +
  base_theme +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 7.0),
    axis.text.y = element_text(size = 7.1, lineheight = 0.92),
    legend.position = "top",
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 6.8),
    panel.background = element_blank(),
    plot.margin = margin(5, 4, 4, 8)
  )

# ---- Panel C: predicted-doublet inclusion sensitivity ---------------------

analysis_set_legend_labels <- c(
  "Primary: 7,371 Scrublet-negative microglia" =
    "Excluded (7,371 cells)",
  "Sensitivity: all 7,461 frozen microglia" =
    "Retained (7,461 cells)"
)
inclusion_panel_plot <- inclusion_panel_full |>
  mutate(
    analysis_set = factor(
      analysis_set,
      levels = analysis_set_levels
    ),
    y_plot = if_else(analysis_set_order == 1, 1.06, 0.94),
    value_label = sprintf("%.3f", value)
  )
inclusion_panel_range <- inclusion_panel_plot |>
  summarise(xmin = min(value), xmax = max(value))

p_c_panel <- ggplot(
  inclusion_panel_plot,
  aes(x = value, y = y_plot, shape = analysis_set)
) +
  geom_segment(
    data = inclusion_panel_range,
    aes(x = xmin, xend = xmax, y = 1, yend = 1),
    inherit.aes = FALSE,
    color = neutral_light,
    linewidth = 0.8
  ) +
  geom_point(size = 2.45, color = neutral_dark, stroke = 0.7) +
  geom_text(
    aes(label = value_label),
    family = font_family,
    size = 2.25,
    color = neutral_dark,
    hjust = 0,
    nudge_x = 0.025
  ) +
  scale_shape_manual(
    values = analysis_set_shapes,
    breaks = analysis_set_levels,
    labels = unname(analysis_set_legend_labels[analysis_set_levels]),
    name = "Cell set (applies to both subplots)"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1.2, by = 0.3),
    limits = c(0, 1.32),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = 1,
    labels = "Seven-transcript\nexpression",
    limits = c(0.68, 1.30),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "C",
    title = "Selected-transcript difference",
    x = paste(
      "Combined exposure − oxygen control",
      "(log2 counts per million)",
      sep = "\n"
    ),
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 8.0, face = "bold", hjust = 0),
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 5, 1, 8)
  )

inclusion_hallmark_plot <- inclusion_hallmark_full |>
  mutate(
    hallmark_label = factor(
      unname(hallmark_labels[ID]),
      levels = rev(unname(hallmark_labels))
    ),
    analysis_set = factor(
      analysis_set,
      levels = analysis_set_levels
    )
  )
inclusion_hallmark_ranges <- inclusion_hallmark_plot |>
  summarise(
    xmin = min(value),
    xmax = max(value),
    .by = hallmark_label
  )

p_c_hallmark <- ggplot(
  inclusion_hallmark_plot,
  aes(x = value, y = hallmark_label, shape = analysis_set)
) +
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
    data = inclusion_hallmark_ranges,
    aes(
      x = xmin,
      xend = xmax,
      y = hallmark_label,
      yend = hallmark_label
    ),
    inherit.aes = FALSE,
    color = neutral_light,
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.45,
    color = neutral_dark,
    stroke = 0.7,
    position = position_dodge(width = 0.18)
  ) +
  scale_shape_manual(
    values = analysis_set_shapes,
    breaks = analysis_set_levels,
    labels = unname(analysis_set_legend_labels[analysis_set_levels]),
    name = "Cell set (applies to both subplots)",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      nrow = 1,
      byrow = TRUE
    )
  ) +
  scale_x_continuous(
    breaks = -3:3,
    limits = c(-3, 3),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(expand = expansion(add = 0.45)) +
  labs(
    title = "Focused Hallmark rankings",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    shape = "Cell set (applies to both subplots)"
  ) +
  base_theme +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 6.8),
    legend.text = element_text(size = 6.6),
    legend.key.width = grid::unit(3.5, "mm"),
    legend.spacing.x = grid::unit(0.8, "mm"),
    plot.title = element_text(size = 8.0, face = "bold", hjust = 0),
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.0, lineheight = 0.92),
    plot.margin = margin(1, 5, 4, 8)
  )

# ---- Assembly --------------------------------------------------------------

figure_design <- "
AABB
CCBB
DDEE
"
figure_3 <-
  p_a_animal + p_b + p_a_loo +
  p_c_panel + p_c_hallmark +
  plot_layout(
    design = figure_design,
    widths = c(0.49, 0.49, 0.51, 0.51),
    heights = c(0.47, 0.57, 0.58)
  )

# ---- Export and execution record ------------------------------------------

base_name <- "Figure3_animal_influence_estimator_sensitivity"
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

save_quartz_pdf(figure_3, pdf_path, width_in, height_in, font_family)

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
capture_render_warnings(print(figure_3))
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
capture_render_warnings(print(figure_3))
grDevices::dev.off()

saveRDS(figure_3, rds_path, compress = "xz")

output_paths <- c(
  vector_pdf_master = pdf_path,
  glia_submission_tiff = tiff_path,
  png_preview = png_path,
  plot_object = rds_path
)

writeLines(
  if (length(render_warnings)) render_warnings else "None",
  file.path(manifest_dir, "fig03_render_warnings.txt")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(manifest_dir, "fig03_plot_sessionInfo.txt")
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
write_csv(input_manifest, file.path(manifest_dir, "fig03_plot_input_manifest.csv"))
write_csv(output_manifest, file.path(manifest_dir, "fig03_plot_output_manifest.csv"))

execution_manifest <- list(
  figure = "Figure 3",
  purpose = paste(
    "Animal influence and estimator sensitivity for selected",
    "interferon-responsive summaries"
  ),
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = basename(script_path),
  script_sha256 = sha256_file(script_path),
  dimensions = list(width_mm = width_mm, height_mm = height_mm),
  preview_suffix = preview_suffix,
  primary_accession = "GSE267933",
  primary_analysis_set_cells = 7371L,
  inclusion_sensitivity_cells = 7461L,
  biological_unit = "animal/library; n=3 per group",
  contrast = "combined exposure minus oxygen control",
  panel_selection_status =
    "selected during exploratory inspection; secondary descriptive analysis",
  leave_one_animal_out_status = paste(
    "descriptive influence diagnostics; not independent hypothesis tests",
    "or biological replication"
  ),
  scales = list(
    direct_expression =
      "animal-level equal-gene-weight mean log2 counts per million",
    score_genes = paste(
      "Scanpy score_genes; seed=42, ctrl_size=50, n_bins=25,",
      "use_raw=FALSE, ctrl_as_ref=TRUE"
    ),
    hallmark = paste(
      "NES from mouse MSigDB 2026.1.Mm Hallmarks ranked by signed",
      "animal-pseudobulk DESeq2 Wald statistics"
    )
  ),
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
    scales = as.character(packageVersion("scales")),
    systemfonts = as.character(packageVersion("systemfonts"))
  )),
  outputs = split(output_manifest, seq_len(nrow(output_manifest)))
)
write_json(
  execution_manifest,
  file.path(manifest_dir, "fig03_plot_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

message("Figure 3 written to: ", output_dir)
message("PDF master: ", basename(pdf_path))
message("TIFF submission file: ", basename(tiff_path), " (", dpi, " dpi, LZW)")
message("PNG preview: ", basename(png_path))
