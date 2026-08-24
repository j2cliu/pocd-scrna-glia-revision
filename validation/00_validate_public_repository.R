#!/usr/bin/env Rscript

# One bounded prepublication audit for the clean GLIA revision repository.

suppressPackageStartupMessages(library(digest))

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Cannot resolve validator path.")
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

script_path <- get_script_path()
repository_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(repository_root)

required <- c(
  "README.md",
  "CITATION.cff",
  "LICENSE",
  "LICENSE-CODE",
  "LICENSE-DOCUMENTATION",
  "analysis/figure1/outputs/Figure1_study_design_partition_audit.tiff",
  "analysis/figure2/outputs/Figure2_transcriptome_and_selected_transcripts.tiff",
  "analysis/figure3/outputs/Figure3_animal_influence_estimator_sensitivity.tiff",
  "analysis/figure4/outputs/Figure4_independent_temporal_context.tiff",
  "analysis/figureS1/outputs/FigureS1_retrospective_partition_annotation_audit.tiff",
  "analysis/figureS2/outputs/FigureS2_partition_stability_diagnostics.tiff",
  "analysis/figureS3/outputs/FigureS3_same_cohort_count_payload_sensitivity.tiff",
  "analysis/tables1_2/outputs/Table1_microglial_analysis_set_audit.csv",
  "analysis/tables1_2/outputs/Table2_partA_per_animal_composition.csv",
  "analysis/tables1_2/outputs/Table2_partB_animal_level_contrasts.csv",
  "analysis/figureS1/outputs/tableS1/TableS1_partA_partition_summary.csv",
  "analysis/figureS2/outputs/tableS2/TableS2_data_dictionary.csv",
  "analysis/figure2/outputs/tableS3/TableS3D_primary_transcriptome_deseq2.csv",
  "analysis/figure3/outputs/tableS4/TableS4_data_dictionary.csv",
  "analysis/figure4/outputs/tableS5/TableS5_data_dictionary.csv",
  "analysis/figureS3/outputs/tableS6/TableS6_data_dictionary.csv",
  "provenance/INCLUSION_MANIFEST.csv",
  "provenance/EXCLUSION_POLICY.csv"
)

missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing required repository files: ", paste(missing, collapse = ", "))
}

inclusion <- read.csv(
  "provenance/INCLUSION_MANIFEST.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_manifest_fields <- c(
  "repository_path", "repository_bytes", "repository_sha256", "transformed"
)
missing_manifest_fields <- setdiff(required_manifest_fields, names(inclusion))
if (length(missing_manifest_fields)) {
  stop(
    "INCLUSION_MANIFEST.csv is missing required fields: ",
    paste(missing_manifest_fields, collapse = ", ")
  )
}
if (anyDuplicated(inclusion$repository_path)) {
  stop("INCLUSION_MANIFEST.csv contains duplicate repository paths.")
}
missing_included <- inclusion$repository_path[!file.exists(inclusion$repository_path)]
if (length(missing_included)) {
  stop(
    "Files listed in INCLUSION_MANIFEST.csv are missing: ",
    paste(missing_included, collapse = ", ")
  )
}
current_inclusion_bytes <- unname(file.info(inclusion$repository_path)$size)
current_inclusion_sha256 <- vapply(
  inclusion$repository_path,
  sha256_file,
  character(1)
)
stale_inclusion <- inclusion$repository_path[
  current_inclusion_bytes != inclusion$repository_bytes |
    current_inclusion_sha256 != inclusion$repository_sha256
]
if (length(stale_inclusion)) {
  stop(
    "INCLUSION_MANIFEST.csv has stale repository identities: ",
    paste(stale_inclusion, collapse = ", ")
  )
}

all_files <- list.files(
  ".",
  recursive = TRUE,
  full.names = FALSE,
  all.files = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)
all_files <- all_files[!grepl("^\\.git(/|$)", all_files)]
all_files <- all_files[
  !grepl("(^|/)__pycache__/|\\.py[co]$|(^|/)\\.DS_Store$", all_files)
]
all_files <- setdiff(
  all_files,
  c(
    "provenance/REPOSITORY_FILE_MANIFEST.csv",
    "provenance/PREPUBLICATION_AUDIT.md"
  )
)

forbidden_scope <- c(
  "non_ifn_discovery",
  "gse222430",
  "gse303920",
  "tnf_ptprs",
  "173_gene"
)
scope_hits <- unlist(lapply(forbidden_scope, function(term) {
  all_files[grepl(term, tolower(all_files), fixed = TRUE)]
}), use.names = FALSE)
if (length(scope_hits)) {
  stop("Excluded scope leaked into repository: ", paste(unique(scope_hits), collapse = ", "))
}

text_extensions <- c(
  "R", "r", "py", "md", "txt", "csv", "tsv", "json", "cff",
  "gitignore"
)
extension <- sub("^.*\\.", "", all_files)
text_files <- all_files[extension %in% text_extensions | basename(all_files) == ".gitignore"]
text_files_to_scan <- setdiff(
  text_files,
  "validation/00_validate_public_repository.R"
)

forbidden_patterns <- c(
  "/Users/",
  "Documents/Lab",
  "Documents/Codex",
  "jasperliumba",
  "crosswalk.passphrase",
  "github_pat_",
  "ghp_",
  "AKIA[0-9A-Z]{16}",
  "-----BEGIN (RSA |EC |OPENSSH )?PRIVATE KEY-----",
  "sk-[A-Za-z0-9]{20,}"
)

text_hits <- character()
for (path in text_files_to_scan) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  for (pattern in forbidden_patterns) {
    if (any(grepl(pattern, lines, perl = TRUE))) {
      text_hits <- c(text_hits, paste(path, pattern, sep = " :: "))
    }
  }
}
if (length(text_hits)) {
  stop("Private path or secret-pattern hits: ", paste(unique(text_hits), collapse = "; "))
}

r_scripts <- all_files[grepl("\\.[Rr]$", all_files)]
r_parse_failures <- character()
for (path in r_scripts) {
  tryCatch(
    parse(file = path),
    error = function(e) {
      r_parse_failures <<- c(r_parse_failures, paste(path, conditionMessage(e), sep = " :: "))
    }
  )
}
if (length(r_parse_failures)) {
  stop("R parse failures: ", paste(r_parse_failures, collapse = "; "))
}

python_scripts <- all_files[grepl("\\.py$", all_files)]
python_status <- if (length(python_scripts)) {
  status <- system2(
    "python3",
    c("-m", "py_compile", python_scripts),
    stdout = TRUE,
    stderr = TRUE
  )
  status_attr <- attr(status, "status")
  if (is.null(status_attr)) 0L else as.integer(status_attr)
} else {
  0L
}
if (!identical(as.integer(python_status), 0L)) {
  stop("At least one Python script failed byte-compilation.")
}

file_info <- file.info(all_files)
if (any(file_info$size >= 100 * 1024^2)) {
  stop("At least one file is at or above GitHub's 100-MB file limit.")
}

manifest <- data.frame(
  repository_path = all_files,
  bytes = unname(file_info$size),
  sha256 = vapply(all_files, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
manifest <- manifest[order(manifest$repository_path), ]
write.csv(
  manifest,
  "provenance/REPOSITORY_FILE_MANIFEST.csv",
  row.names = FALSE,
  na = ""
)

largest <- manifest[order(-manifest$bytes), ][seq_len(min(10L, nrow(manifest))), ]
report <- c(
  "# Prepublication repository audit",
  "",
  paste0("- Status: **PASS**"),
  paste0("- Files audited: ", nrow(manifest)),
  paste0("- Total bytes: ", format(sum(manifest$bytes), big.mark = ",", scientific = FALSE)),
  paste0("- R scripts parsed: ", length(r_scripts)),
  paste0("- Python scripts byte-compiled: ", length(python_scripts)),
  paste0("- Required artifacts present: ", length(required), "/", length(required)),
  paste0("- Inclusion-manifest identities verified: ", nrow(inclusion), "/", nrow(inclusion)),
  "- Private absolute paths: 0",
  "- Credential/secret-pattern hits: 0",
  "- Files at or above 100 MB: 0",
  "- Excluded-scope path hits: 0",
  "",
  "## Largest files",
  "",
  "| Repository path | Bytes | SHA-256 |",
  "|---|---:|---|",
  apply(largest, 1L, function(row) {
    paste0("| `", row[["repository_path"]], "` | ", row[["bytes"]], " | `", row[["sha256"]], "` |")
  }),
  "",
  "This audit checks repository integrity and public-release hygiene. It does",
  "not constitute biological replication or orthogonal validation."
)
writeLines(report, "provenance/PREPUBLICATION_AUDIT.md", useBytes = TRUE)

cat(
  "PASS: ", nrow(manifest), " files; ", length(r_scripts),
  " R scripts; ", length(python_scripts), " Python scripts; ",
  " no private paths or secret-pattern hits.\n",
  sep = ""
)
