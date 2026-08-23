#!/usr/bin/env Rscript

# GLIA major revision — complete primary pseudobulk result for Table S3D–E
#
# Run from the revision workspace root:
#   Rscript revision_work/results_rebuild/figure2/scripts/95_build_glia_R1_tableS3_primary_transcriptome.R

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
audit_root <- file.path(project_root, "independent_audit")
output_dir <- file.path(
  project_root,
  "revision_work/results_rebuild/figure2/outputs/tableS3"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

full_path <- file.path(audit_root, "scrublet_negative_deseq2_full.csv")
loo_path <- file.path(audit_root, "scrublet_negative_deseq2_loo.csv.gz")
count_path <- file.path(audit_root, "scrublet_negative_pseudobulk_counts.csv.gz")
stopifnot(file.exists(full_path), file.exists(loo_path), file.exists(count_path))

full <- read.csv(full_path, check.names = FALSE, na.strings = "NA")
loo <- read.csv(gzfile(loo_path), check.names = FALSE, na.strings = "NA")
counts <- read.csv(gzfile(count_path), check.names = FALSE)

required_full <- c(
  "gene", "symbol", "baseMean", "log2FoldChange", "lfcSE", "ci95_low",
  "ci95_high", "stat", "pvalue", "padj"
)
required_loo <- c("dropped_animal", "n_surgery", "n_control", required_full)
stopifnot(
  all(required_full %in% names(full)),
  all(required_loo %in% names(loo)),
  identical(names(counts), c("gene", "C1", "C2", "C3", "S1", "S2", "S3")),
  nrow(full) == 13926L,
  !anyDuplicated(full$gene),
  !anyDuplicated(counts$gene)
)

animal_order <- c("C1", "C2", "C3", "S1", "S2", "S3")
stopifnot(
  setequal(unique(loo$dropped_animal), animal_order),
  all(table(loo$dropped_animal) == nrow(full)),
  all(vapply(split(loo$gene, loo$dropped_animal), function(x) !anyDuplicated(x), logical(1)))
)

count_index <- match(full$gene, counts$gene)
stopifnot(!anyNA(count_index))
raw_counts <- counts[count_index, animal_order, drop = FALSE]
names(raw_counts) <- paste0("raw_umi_", animal_order)

full_fdr <- !is.na(full$padj) & full$padj < 0.05
full_sign <- sign(full$log2FoldChange)

loo_fdr_count <- integer(nrow(full))
loo_fdr_same_sign_count <- integer(nrow(full))
for (animal in animal_order) {
  one <- loo[loo$dropped_animal == animal, , drop = FALSE]
  one <- one[match(full$gene, one$gene), , drop = FALSE]
  stopifnot(identical(full$gene, one$gene))
  one_fdr <- !is.na(one$padj) & one$padj < 0.05
  loo_fdr_count <- loo_fdr_count + as.integer(one_fdr)
  loo_fdr_same_sign_count <- loo_fdr_same_sign_count + as.integer(
    one_fdr & sign(one$log2FoldChange) == full_sign
  )
}

part_d <- cbind(
  data.frame(
    nominal_p_rank = rank(full$pvalue, ties.method = "min", na.last = "keep"),
    ensembl_id = full$gene,
    symbol = full$symbol,
    raw_counts,
    base_mean = full$baseMean,
    log2_fold_change = full$log2FoldChange,
    lfc_se = full$lfcSE,
    ci95_low = full$ci95_low,
    ci95_high = full$ci95_high,
    wald_statistic = full$stat,
    nominal_p_value = full$pvalue,
    bh_adjusted_p_value = full$padj,
    bh_fdr_lt_0_05 = full_fdr,
    effect_direction = ifelse(
      full$log2FoldChange > 0,
      "More positive in combined exposure",
      ifelse(
        full$log2FoldChange < 0,
        "More positive in oxygen control",
        "No directional difference"
      )
    ),
    n_leave_one_animal_out_fdr_lt_0_05 = loo_fdr_count,
    n_leave_one_animal_out_fdr_lt_0_05_same_direction = loo_fdr_same_sign_count,
    fdr_lt_0_05_in_all_six_leave_one_animal_out_fits = (
      loo_fdr_count == length(animal_order) &
        loo_fdr_same_sign_count == length(animal_order)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)

part_d <- part_d[order(part_d$nominal_p_rank, part_d$ensembl_id), , drop = FALSE]
rownames(part_d) <- NULL

summarize_fit <- function(data, scenario, dropped_animal = NA_character_) {
  fdr <- !is.na(data$padj) & data$padj < 0.05
  retained_full_hits <- if (is.na(dropped_animal)) {
    sum(full_fdr)
  } else {
    sum(full_fdr[match(data$gene, full$gene)] & fdr)
  }
  data.frame(
    scenario = scenario,
    dropped_animal = dropped_animal,
    n_combined_exposure_animals = if (is.na(dropped_animal)) 3L else unique(data$n_surgery),
    n_oxygen_control_animals = if (is.na(dropped_animal)) 3L else unique(data$n_control),
    n_tested_features = nrow(data),
    n_nonmissing_bh_adjusted_p_values = sum(!is.na(data$padj)),
    n_bh_fdr_lt_0_05 = sum(fdr),
    n_positive_bh_fdr_lt_0_05 = sum(fdr & data$log2FoldChange > 0),
    n_negative_bh_fdr_lt_0_05 = sum(fdr & data$log2FoldChange < 0),
    n_full_fit_fdr_hits_retained = retained_full_hits,
    stringsAsFactors = FALSE
  )
}

part_e <- summarize_fit(full, "Full cohort")
for (animal in animal_order) {
  part_e <- rbind(
    part_e,
    summarize_fit(
      loo[loo$dropped_animal == animal, , drop = FALSE],
      paste("Omit", animal),
      animal
    )
  )
}
rownames(part_e) <- NULL

stopifnot(
  sum(full_fdr) == 40L,
  sum(full_fdr & full$log2FoldChange > 0) == 11L,
  sum(full_fdr & full$log2FoldChange < 0) == 29L,
  identical(part_e$n_bh_fdr_lt_0_05, c(40L, 18L, 18L, 10L, 15L, 58L, 191L)),
  identical(part_e$n_full_fit_fdr_hits_retained, c(40L, 17L, 11L, 8L, 11L, 24L, 34L)),
  sum(part_d$fdr_lt_0_05_in_all_six_leave_one_animal_out_fits) == 5L,
  setequal(
    na.omit(part_d$symbol[part_d$fdr_lt_0_05_in_all_six_leave_one_animal_out_fits]),
    c("Ccl9", "H1f2", "Lst1", "Tnfaip3", "Ttr")
  )
)

write_pair <- function(data, stem) {
  write.csv(
    data,
    file.path(output_dir, paste0(stem, ".csv")),
    row.names = FALSE,
    na = "NA",
    quote = TRUE
  )
  write.table(
    data,
    file.path(output_dir, paste0(stem, ".tsv")),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    quote = FALSE
  )
}

write_pair(part_d, "TableS3D_primary_transcriptome_deseq2")
write_pair(part_e, "TableS3E_primary_transcriptome_loo_summary")

cat(
  sprintf(
    paste0(
      "PASS: Table S3D contains %d primary features (%d BH FDR < 0.05; ",
      "%d positive, %d negative); Table S3E contains the full fit and six ",
      "animal-omission summaries.\n"
    ),
    nrow(part_d),
    sum(full_fdr),
    sum(full_fdr & full$log2FoldChange > 0),
    sum(full_fdr & full$log2FoldChange < 0)
  )
)
