#!/usr/bin/env Rscript

# Convert the frozen internal evidence ledger into a clean repository-facing
# map. The historical retired IEG row is omitted, final Table S4 numbering is
# used, and every listed code/data path must exist in this repository.

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Cannot resolve script path.")
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

repository_root <- normalizePath(
  file.path(dirname(get_script_path()), ".."),
  mustWork = TRUE
)
setwd(repository_root)

ledger <- read.csv(
  "provenance/RESULTS_EVIDENCE_LEDGER_INTERNAL.csv",
  check.names = FALSE
)
ledger <- ledger[ledger$ledger_id != "R4.P1", , drop = FALSE]

code_map <- list(
  "R1.P1" = c(
    "analysis/figure1/scripts/78_prepare_glia_R1_figure1_data.R",
    "analysis/tables1_2/scripts/80_prepare_glia_R1_tables1_2.R"
  ),
  "R1.P2" = c(
    "analysis/figureS1/scripts/91_prepare_glia_R1_figureS1_tableS1.py",
    "analysis/figureS1/scripts/92_plot_glia_R1_figureS1.R"
  ),
  "R1.P3" = c(
    "analysis/figureS2/scripts/91_prepare_glia_R1_figureS2_tableS2.R",
    "analysis/figureS2/scripts/92_plot_glia_R1_figureS2.R"
  ),
  "R1.P4" = c(
    "analysis/figure1/scripts/78_prepare_glia_R1_figure1_data.R",
    "analysis/tables1_2/scripts/80_prepare_glia_R1_tables1_2.R"
  ),
  "R2.P0" = c(
    "analysis/figure2/scripts/95_build_glia_R1_tableS3_primary_transcriptome.R",
    "analysis/figureS1/scripts/91_prepare_glia_R1_figureS1_tableS1.py"
  ),
  "R2.P1" = c(
    "analysis/figure2/scripts/82_prepare_glia_R1_figure2_data.R",
    "analysis/figure2/scripts/84_plot_glia_R1_figure2_rebalanced.R"
  ),
  "R2.P2" = c(
    "analysis/figure2/scripts/82_prepare_glia_R1_figure2_data.R",
    "analysis/figure2/scripts/84_plot_glia_R1_figure2_rebalanced.R"
  ),
  "R2.P3" = c(
    "analysis/figure2/scripts/82_prepare_glia_R1_figure2_data.R",
    "analysis/figure2/scripts/84_plot_glia_R1_figure2_rebalanced.R"
  ),
  "R3.P1" = c(
    "analysis/source_audit/scripts/audit_score_genes_fixed_component_decomposition.py",
    "analysis/source_audit/scripts/audit_score_genes_doublet_negative.py",
    "analysis/figure3/scripts/85_plot_glia_R1_figure3.R"
  ),
  "R3.P2" = c(
    "analysis/source_audit/scripts/audit_scrublet_negative_deseq_gsea.R",
    "analysis/figure3/scripts/85_plot_glia_R1_figure3.R"
  ),
  "R3.P3" = c(
    "analysis/source_audit/scripts/audit_scrublet_negative_deseq_gsea.R",
    "analysis/figure3/scripts/85_plot_glia_R1_figure3.R"
  ),
  "R4.P2" = c(
    "analysis/source_audit/scripts/audit_scrublet_negative_deseq_gsea.R",
    "analysis/figure3/scripts/85_plot_glia_R1_figure3.R"
  ),
  "R4.P3" = c(
    "analysis/figureS3/scripts/88_preflight_gse289098_common_payload.py",
    "analysis/figureS3/scripts/89_prepare_figureS3_tableS6.R",
    "analysis/figureS3/scripts/90_plot_figureS3.R"
  ),
  "R5.P1" = c(
    "analysis/figure4/scripts/86_prepare_glia_R1_figure4_data.R",
    "analysis/figure4/scripts/87_plot_glia_R1_figure4.R"
  ),
  "R5.P2" = c(
    "analysis/figure4/scripts/86_prepare_glia_R1_figure4_data.R",
    "analysis/figure4/scripts/87_plot_glia_R1_figure4.R"
  ),
  "R5.P3" = c(
    "analysis/figure4/scripts/86_prepare_glia_R1_figure4_data.R",
    "analysis/figure4/scripts/87_plot_glia_R1_figure4.R"
  )
)

table_s4_rename <- c(
  TableS4B_score_parameter_sensitivity.csv =
    "TableS4B_selected_panel_score_parameter_sensitivity.csv",
  TableS4C_score_parameter_summary.csv =
    "TableS4C_selected_panel_score_parameter_summary.csv",
  TableS4D_score_leave_one_gene_out.csv =
    "TableS4D_score_leave_one_transcript_out.csv",
  TableS4H_cell_inclusion_selected_panel.csv =
    "TableS4E_cell_inclusion_selected_panel.csv",
  TableS4I_cell_inclusion_hallmark.csv =
    "TableS4F_cell_inclusion_hallmark.csv",
  TableS4J_fixed_isg_score_animal_values.csv =
    "TableS4G_fixed_selected_panel_score_animal_values.csv",
  TableS4K_fixed_isg_score_effects_full_and_animal_omissions.csv =
    "TableS4H_fixed_selected_panel_score_effects_full_and_animal_omissions.csv",
  TableS4L_fixed_selected_panel_score_control_gene_manifest.csv =
    "TableS4I_fixed_selected_panel_score_control_gene_manifest.csv",
  TableS4M_fixed_selected_panel_score_component_animal_values.csv =
    "TableS4J_fixed_selected_panel_score_component_animal_values.csv",
  TableS4N_fixed_selected_panel_score_component_effects_full_and_animal_omissions.csv =
    "TableS4K_fixed_selected_panel_score_component_effects_full_and_animal_omissions.csv"
)

all_repository_files <- list.files(
  ".",
  recursive = TRUE,
  full.names = FALSE,
  all.files = FALSE,
  include.dirs = FALSE
)

locate_source <- function(path) {
  base <- basename(trimws(path))
  if (base %in% names(table_s4_rename)) {
    base <- unname(table_s4_rename[[base]])
  }
  candidates <- all_repository_files[basename(all_repository_files) == base]
  candidates <- candidates[!grepl("^provenance/RESULTS_EVIDENCE", candidates)]
  if (!length(candidates)) {
    stop("No repository file found for evidence source: ", base)
  }
  preferred <- candidates[
    grepl("/outputs/|/data/panel_ready/", candidates)
  ]
  if (length(preferred)) candidates <- preferred
  if (length(candidates) > 1L) {
    candidates <- candidates[order(nchar(candidates), candidates)]
  }
  candidates[[1L]]
}

map_source_field <- function(field) {
  pieces <- trimws(unlist(strsplit(field, "\\|")))
  paste(vapply(pieces, locate_source, character(1)), collapse = " | ")
}

public <- data.frame(
  evidence_id = ledger$ledger_id,
  results_section = ledger$results_section,
  paragraph_heading = ledger$paragraph_heading,
  claim_level = ledger$claim_level,
  permitted_claim = ledger$permitted_claim,
  exact_numeric_anchors = ledger$exact_numeric_anchors,
  repository_data = vapply(ledger$source_data, map_source_field, character(1)),
  repository_code = vapply(ledger$ledger_id, function(id) {
    paths <- code_map[[id]]
    if (is.null(paths) || !length(paths)) stop("Missing code map for ", id)
    if (any(!file.exists(paths))) {
      stop("Mapped code file missing for ", id)
    }
    paste(paths, collapse = " | ")
  }, character(1)),
  display_cross_reference = ledger$display_cross_reference,
  required_caveat = ledger$required_caveat,
  prohibited_wording = ledger$prohibited_wording,
  package_status = "FROZEN_FOR_PREPUBLICATION_REVIEW",
  stringsAsFactors = FALSE
)

write.csv(
  public,
  "provenance/RESULTS_EVIDENCE_MAP.csv",
  row.names = FALSE,
  na = ""
)

cat("Wrote public evidence map with ", nrow(public), " active rows.\n", sep = "")

