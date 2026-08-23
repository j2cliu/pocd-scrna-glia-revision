#!/usr/bin/env Rscript

# GLIA major revision — canonical Figure S1 plotting script
#
# Retrospective audit of the submitted seven-way partition.  Traceability
# labels are not treated as validated biological states.  All plotting and
# multipanel assembly are performed in R.

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

italic_gene_labels <- function(values) {
  lapply(values, function(value) bquote(italic(.(value))))
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
height_mm <- as.numeric(get_arg("--height-mm", "150"))
dpi <- as.integer(get_arg("--dpi", "300"))
assert_true(
  width_mm == 180 && height_mm == 150,
  "Canonical Figure S1 dimensions are frozen at 180 x 150 mm."
)
assert_true(dpi == 300L, "Submission TIFF is frozen at 300 dpi.")

input_paths <- c(
  metadata = file.path(panel_dir, "figS01_metadata.csv"),
  selected = file.path(panel_dir, "figS01_selected_gene_dotplot.csv"),
  doublets = file.path(panel_dir, "figS01_doublet_summary.csv"),
  detection = file.path(panel_dir, "figS01_detection_heatmap.csv")
)
missing_inputs <- input_paths[!file.exists(input_paths)]
assert_true(
  length(missing_inputs) == 0L,
  paste("Missing panel-ready input(s):", paste(missing_inputs, collapse = ", "))
)

metadata <- read_csv(input_paths[["metadata"]], show_col_types = FALSE)
selected <- read_csv(input_paths[["selected"]], show_col_types = FALSE)
doublets <- read_csv(input_paths[["doublets"]], show_col_types = FALSE)
detection <- read_csv(input_paths[["detection"]], show_col_types = FALSE)

assert_columns(
  metadata,
  c(
    "n_submitted_cells", "n_submitted_partitions",
    "n_detected_symbol_deduplicated_genes", "claim_ceiling",
    "published_signature_comparison"
  ),
  "figS01_metadata.csv"
)
assert_columns(
  selected,
  c(
    "partition", "traceability_name", "gene", "gene_order", "n_cells",
    "mean_log1p_cp10k", "pct_detected", "selection_status"
  ),
  "figS01_selected_gene_dotplot.csv"
)
assert_columns(
  doublets,
  c(
    "partition", "traceability_name", "n_cells", "n_doublets",
    "doublet_rate_percent", "interpretation_ceiling"
  ),
  "figS01_doublet_summary.csv"
)
assert_columns(
  detection,
  c(
    "partition", "traceability_name", "gene_group", "gene_group_order",
    "gene", "within_group_order", "n_cells", "pct_detected"
  ),
  "figS01_detection_heatmap.csv"
)

partitions <- 0:6
selected_genes <- c(
  "P2ry12", "Tmem119", "Cx3cr1", "Tnf", "Il1b", "Irf7", "Ifitm3",
  "Nr1d1", "Dbp"
)
assert_true(
  nrow(metadata) == 1L &&
    metadata$n_submitted_cells == 7461L &&
    metadata$n_submitted_partitions == 7L &&
    metadata$n_detected_symbol_deduplicated_genes == 17878L &&
    metadata$published_signature_comparison == "withdrawn",
  "Figure S1 metadata fail the frozen audit anchors."
)
assert_true(
  nrow(selected) == 63L &&
    setequal(selected$partition, partitions) &&
    setequal(selected$gene, selected_genes),
  "Figure S1A selected-gene matrix is incomplete."
)
assert_true(
  nrow(doublets) == 7L &&
    sum(doublets$n_cells) == 7461L &&
    sum(doublets$n_doublets) == 90L &&
    doublets$n_doublets[doublets$partition == 6L] == 82L &&
    doublets$n_cells[doublets$partition == 6L] == 94L,
  "Figure S1B doublet audit fails the verified anchors."
)
assert_true(
  nrow(detection) == 112L &&
    nrow(distinct(detection, partition, gene)) == 112L &&
    detection$pct_detected[detection$partition == 6L & detection$gene == "Plp1"] > 99 &&
    max(detection$pct_detected[detection$gene_group == "BAM/CAM"]) < 2,
  "Figure S1C detection audit fails identity/contamination anchors."
)

font_match <- systemfonts::match_fonts("Arial")
font_family <- if (nrow(font_match) > 0L && nzchar(font_match$path[[1L]])) {
  "Arial"
} else {
  "Helvetica"
}
if (font_family != "Arial") warning("Arial was not found; using Helvetica.")

neutral_dark <- "#303030"
neutral_mid <- "#707070"
neutral_light <- "#D0D0D0"
neutral_faint <- "#EFEFEF"
doublet_color <- "#6A51A3"
rare_outline <- "#A63603"
sequential_low <- "#F7FBFF"
sequential_mid <- "#6BAED6"
sequential_high <- "#08306B"

base_theme <- theme_classic(base_family = font_family, base_size = 8.2) +
  theme(
    text = element_text(color = neutral_dark),
    plot.tag = element_text(size = 11, face = "bold", color = "#111111"),
    plot.tag.position = c(0, 1),
    axis.title = element_text(size = 7.8),
    axis.text = element_text(size = 7.1, color = neutral_dark),
    strip.text = element_text(size = 7.0, face = "bold"),
    legend.title = element_text(size = 7.0),
    legend.text = element_text(size = 6.8),
    plot.margin = margin(5, 6, 5, 6)
  )

# Panel A: the exact nine-gene panel used retrospectively to describe the
# submitted partition.  These genes were author-selected rather than selected
# by the marker test.
selected <- selected |>
  mutate(
    partition_position = partition + 1,
    gene = factor(gene, levels = rev(selected_genes)),
    gene_family = case_when(
      as.character(gene) %in% c("P2ry12", "Tmem119", "Cx3cr1") ~ "Microglia core",
      as.character(gene) %in% c("Tnf", "Il1b", "Irf7", "Ifitm3") ~ "Inflammatory / IFN",
      TRUE ~ "Clock-related"
    ),
    gene_family = factor(
      gene_family,
      levels = c("Microglia core", "Inflammatory / IFN", "Clock-related")
    )
  )

p_a <- ggplot(
  selected,
  aes(
    x = partition_position,
    y = gene,
    size = pct_detected,
    color = mean_log1p_cp10k
  )
) +
  annotate(
    "rect",
    xmin = 6.55, xmax = 7.45, ymin = -Inf, ymax = Inf,
    fill = NA, color = rare_outline, linewidth = 0.45, linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(2.5, 6.5),
    color = neutral_light,
    linewidth = 0.35
  ) +
  geom_point(alpha = 0.95) +
  scale_x_continuous(
    breaks = 1:7,
    labels = c("0", "1", "2", "3", "4", "5", "6†"),
    limits = c(0.55, 7.45),
    expand = expansion(mult = 0)
  ) +
  scale_y_discrete(labels = italic_gene_labels) +
  scale_size_continuous(
    name = "Cells detected (%)",
    range = c(0.5, 4.2),
    breaks = c(25, 50, 75, 100),
    limits = c(0, 100)
  ) +
  scale_color_gradientn(
    name = "Mean log1p CP10K",
    colors = c(sequential_low, sequential_mid, sequential_high),
    values = scales::rescale(c(0, 0.8, 3.4)),
    limits = c(0, 3.4),
    oob = scales::squish
  ) +
  labs(x = "Submitted numeric partition", y = NULL, tag = "A") +
  base_theme +
  theme(
    panel.grid.major = element_line(color = neutral_faint, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(9, 6, 5, 9)
  )

# Panel B: Scrublet calls, emphasizing that partition 6 is overwhelmingly
# doublet-enriched while the other submitted partitions are not.
doublets <- doublets |>
  mutate(
    partition_factor = factor(partition, levels = rev(0:6)),
    label = paste0(n_doublets, "/", n_cells),
    label_x = if_else(doublet_rate_percent > 75, doublet_rate_percent - 3, pmax(3, doublet_rate_percent + 3)),
    label_hjust = if_else(doublet_rate_percent > 75, 1, 0)
  )

p_b <- ggplot(doublets, aes(y = partition_factor)) +
  geom_segment(
    aes(x = 0, xend = doublet_rate_percent, yend = partition_factor),
    color = neutral_light, linewidth = 0.8, lineend = "round"
  ) +
  geom_point(
    aes(x = doublet_rate_percent),
    shape = 21, size = 2.7, stroke = 0.55,
    fill = doublet_color, color = neutral_dark
  ) +
  geom_text(
    aes(x = label_x, label = label, hjust = label_hjust),
    family = font_family, size = 2.35, color = neutral_dark
  ) +
  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(labels = function(values) ifelse(values == "6", "6†", values)) +
  labs(
    x = "Scrublet-predicted doublets (%)",
    y = "Submitted numeric partition",
    tag = "B"
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = neutral_faint, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Panel C: raw-count detection audit.  Ttr and oligodendrocyte-lineage signal
# are shown as contamination/ambient diagnostics, not partition identities.
detection <- detection |>
  arrange(gene_group_order, within_group_order) |>
  mutate(
    partition_position = partition + 1,
    gene_group = factor(
      gene_group,
      levels = c(
        "Microglia core", "BAM/CAM", "Choroid-plexus-associated",
        "Oligodendrocyte lineage"
      )
    )
  )
gene_levels <- detection |>
  distinct(gene, gene_group_order, within_group_order) |>
  arrange(gene_group_order, within_group_order) |>
  pull(gene)
detection$gene <- factor(detection$gene, levels = rev(gene_levels))

p_c <- ggplot(detection, aes(x = partition_position, y = gene, fill = pct_detected)) +
  geom_tile(color = "white", linewidth = 0.35) +
  annotate(
    "rect",
    xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
    fill = NA, color = rare_outline, linewidth = 0.55, linetype = "dashed"
  ) +
  facet_grid(
    gene_group ~ ., scales = "free_y", space = "free_y", switch = "y"
  ) +
  scale_x_continuous(
    breaks = 1:7,
    labels = c("0", "1", "2", "3", "4", "5", "6†"),
    limits = c(0.5, 7.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_discrete(labels = italic_gene_labels) +
  scale_fill_gradientn(
    name = "Cells detected (%)",
    colors = c(sequential_low, sequential_mid, sequential_high),
    values = scales::rescale(c(0, 25, 100)),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  labs(x = "Submitted numeric partition", y = NULL, tag = "C") +
  base_theme +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, color = neutral_mid),
    legend.position = "right"
  )

figure_plot <- (
  p_a + p_b + plot_layout(widths = c(1.35, 0.80), guides = "keep")
) / p_c +
  plot_layout(heights = c(1.03, 1.12))

output_paths <- c(
  pdf = file.path(output_dir, "FigureS1_retrospective_partition_annotation_audit.pdf"),
  tiff = file.path(output_dir, "FigureS1_retrospective_partition_annotation_audit.tiff"),
  preview = file.path(output_dir, "FigureS1_retrospective_partition_annotation_audit_preview.png"),
  rds = file.path(output_dir, "FigureS1_retrospective_partition_annotation_audit.rds")
)

render_warnings <- character()
capture_render_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(warning_condition) {
      render_warnings <<- unique(c(
        render_warnings,
        conditionMessage(warning_condition)
      ))
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

withCallingHandlers(
  {
    save_quartz_pdf(
      figure_plot,
      output_paths[["pdf"]],
      width_mm / 25.4,
      height_mm / 25.4,
      font_family
    )

    agg_tiff(
      filename = output_paths[["tiff"]],
      width = width_mm,
      height = height_mm,
      units = "mm",
      res = dpi,
      compression = "lzw",
      background = "white"
    )
    print(figure_plot)
    dev.off()

    agg_png(
      filename = output_paths[["preview"]],
      width = width_mm,
      height = height_mm,
      units = "mm",
      res = 150,
      background = "white"
    )
    print(figure_plot)
    dev.off()
  },
  warning = function(warning_condition) {
    render_warnings <<- c(render_warnings, conditionMessage(warning_condition))
    invokeRestart("muffleWarning")
  }
)
saveRDS(figure_plot, output_paths[["rds"]])

writeLines(
  if (length(render_warnings)) unique(render_warnings) else "No render warnings.",
  file.path(manifest_dir, "figS01_render_warnings.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(manifest_dir, "figS01_plot_sessionInfo.txt")
)

input_manifest <- tibble(
  role = names(input_paths),
  path = unname(input_paths),
  sha256 = vapply(unname(input_paths), sha256_file, character(1)),
  size_bytes = file.info(unname(input_paths))$size
)
write_csv(input_manifest, file.path(manifest_dir, "figS01_plot_input_manifest.csv"))

output_manifest <- tibble(
  role = names(output_paths),
  path = unname(output_paths),
  sha256 = vapply(unname(output_paths), sha256_file, character(1)),
  size_bytes = file.info(unname(output_paths))$size
)
write_csv(output_manifest, file.path(manifest_dir, "figS01_plot_output_manifest.csv"))

execution <- list(
  analysis = "Figure S1 retrospective submitted-partition annotation audit",
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = script_path,
  script_sha256 = sha256_file(script_path),
  dimensions_mm = list(width = width_mm, height = height_mm),
  submission_tiff_dpi = dpi,
  submission_tiff_compression = "LZW",
  font_family = font_family,
  anchors = list(
    n_cells = metadata$n_submitted_cells[[1]],
    n_partitions = metadata$n_submitted_partitions[[1]],
    n_scrublet_calls = sum(doublets$n_doublets),
    partition_6_scrublet_calls = doublets$n_doublets[doublets$partition == 6L],
    partition_6_cells = doublets$n_cells[doublets$partition == 6L],
    max_bam_cam_detection_percent = max(
      detection$pct_detected[detection$gene_group == "BAM/CAM"]
    )
  ),
  claim_ceiling = metadata$claim_ceiling[[1]],
  render_warnings = unique(render_warnings)
)
write_json(
  execution,
  file.path(manifest_dir, "figS01_plot_execution.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

message("Figure S1 written to: ", output_dir)
