#!/usr/bin/env Rscript

# Re-render all seven final manuscript displays in a temporary directory and
# compare TIFF/PNG hashes with the frozen repository artifacts. PDF hashes are
# not required because device metadata can vary by render time.

suppressPackageStartupMessages(library(digest))

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Cannot resolve script path.")
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

copy_tree <- function(source, destination) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(
    source,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE,
    include.dirs = FALSE
  )
  relative <- substring(entries, nchar(source) + 2L)
  for (i in seq_along(entries)) {
    target <- file.path(destination, relative[[i]])
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(entries[[i]], target, overwrite = FALSE)) {
      stop("Failed to copy QA input: ", entries[[i]])
    }
  }
}

script_path <- get_script_path()
repository_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(repository_root)

qa_root <- tempfile("glia_display_reproduction_")
dir.create(qa_root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(qa_root, recursive = TRUE, force = TRUE), add = TRUE)
file.copy(
  file.path(repository_root, "analysis", "RESULTS_DISPLAY_MAP_V2_20260820.md"),
  qa_root,
  overwrite = FALSE
)

jobs <- data.frame(
  analysis_id = c(
    "figure1", "figure2", "figure3", "figure4",
    "figureS1", "figureS2", "figureS3"
  ),
  script = c(
    "79_plot_glia_R1_figure1.R",
    "84_plot_glia_R1_figure2_rebalanced.R",
    "85_plot_glia_R1_figure3.R",
    "87_plot_glia_R1_figure4.R",
    "92_plot_glia_R1_figureS1.R",
    "92_plot_glia_R1_figureS2.R",
    "90_plot_figureS3.R"
  ),
  output_base = c(
    "Figure1_study_design_partition_audit",
    "Figure2_transcriptome_and_selected_transcripts",
    "Figure3_animal_influence_estimator_sensitivity",
    "Figure4_independent_temporal_context",
    "FigureS1_retrospective_partition_annotation_audit",
    "FigureS2_partition_stability_diagnostics",
    "FigureS3_same_cohort_count_payload_sensitivity"
  ),
  stringsAsFactors = FALSE
)

results <- list()

for (i in seq_len(nrow(jobs))) {
  id <- jobs$analysis_id[[i]]
  source_root <- file.path(repository_root, "analysis", id)
  target_root <- file.path(qa_root, id)

  copy_tree(file.path(source_root, "data"), file.path(target_root, "data"))
  copy_tree(file.path(source_root, "scripts"), file.path(target_root, "scripts"))
  dir.create(file.path(target_root, "outputs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(target_root, "manifests"), recursive = TRUE, showWarnings = FALSE)

  if (id == "figure2") {
    table_dir <- file.path(target_root, "outputs", "tableS3")
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(
      file.path(
        source_root,
        "outputs",
        "tableS3",
        "TableS3D_primary_transcriptome_deseq2.csv"
      ),
      table_dir,
      overwrite = FALSE
    )
  }

  plot_script <- file.path(target_root, "scripts", jobs$script[[i]])
  log_path <- file.path(qa_root, paste0(id, ".log"))
  status <- system2(
    "Rscript",
    plot_script,
    stdout = log_path,
    stderr = log_path
  )
  if (!identical(as.integer(status), 0L)) {
    stop(
      "Display reproduction failed for ", id, ". Log: ",
      paste(readLines(log_path, warn = FALSE), collapse = " | ")
    )
  }

  for (extension in c("tiff", "png")) {
    suffix <- if (extension == "png") "_preview.png" else ".tiff"
    rendered <- file.path(target_root, "outputs", paste0(jobs$output_base[[i]], suffix))
    frozen <- file.path(source_root, "outputs", paste0(jobs$output_base[[i]], suffix))
    if (!file.exists(rendered) || !file.exists(frozen)) {
      stop("Missing rendered or frozen display for ", id, " (", extension, ")")
    }
    rendered_hash <- sha256_file(rendered)
    frozen_hash <- sha256_file(frozen)
    results[[length(results) + 1L]] <- data.frame(
      analysis_id = id,
      format = toupper(extension),
      rendered_sha256 = rendered_hash,
      frozen_sha256 = frozen_hash,
      identical = identical(rendered_hash, frozen_hash),
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, results)
write.csv(
  results,
  "provenance/DISPLAY_REPRODUCTION_HASHES.csv",
  row.names = FALSE
)

failed_tiff <- results[results$format == "TIFF" & !results$identical, , drop = FALSE]
if (nrow(failed_tiff)) {
  stop(
    "One or more submission TIFF hashes differed: ",
    paste(failed_tiff$analysis_id, collapse = ", ")
  )
}

n_tiff_identical <- sum(results$format == "TIFF" & results$identical)
n_png_identical <- sum(results$format == "PNG" & results$identical)

report <- c(
  "# Frozen-display reproduction audit",
  "",
  "- Status: **PASS**",
  "- Final displays re-rendered: 7",
  paste0("- Submission TIFF hashes: ", n_tiff_identical, "/7 identical"),
  paste0("- Convenience PNG preview hashes: ", n_png_identical, "/7 byte-identical"),
  "- Execution location: temporary directory removed after audit",
  "",
  "The 300-dpi TIFF files are the journal-submission raster artifacts and are",
  "the release-blocking comparison. PNG files are convenience previews; two",
  "frozen previews were resized or re-exported after the canonical TIFF render",
  "and are therefore reported but not used as a release criterion. PDF hashes",
  "were not used because PDF device metadata can vary by render time.",
  "",
  "This audit verifies deterministic reconstruction of the final submission",
  "rasters from repository-contained derived inputs; it does not rerun the",
  "upstream biological analyses or add independent validation."
)
writeLines(report, "provenance/DISPLAY_REPRODUCTION_AUDIT.md", useBytes = TRUE)

cat(
  "PASS: 7 displays re-rendered; ", n_tiff_identical,
  "/7 submission TIFF hashes identical; ", n_png_identical,
  "/7 PNG preview hashes byte-identical.\n",
  sep = ""
)
