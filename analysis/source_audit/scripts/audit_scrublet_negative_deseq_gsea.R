#!/usr/bin/env Rscript

# Scrublet-negative GSE267933 pseudobulk DESeq2 and mouse-native Hallmark audit.
# The input counts have already been filtered at >=10 total UMI in the full
# six-animal cohort. Every leave-one-animal-out fit retains that frozen universe.
# Mouse-native MSigDB is mandatory: this script has no human-Hallmark fallback.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(clusterProfiler)
  library(msigdbr)
})

project <- Sys.getenv("POCD_SCRNA_PROJECT_ROOT")
if (!nzchar(project)) stop("Set POCD_SCRNA_PROJECT_ROOT before running this audit.")
project <- normalizePath(project, mustWork = TRUE)
out_dir <- normalizePath(dirname(sub(
  "^--file=", "",
  grep("^--file=", commandArgs(FALSE), value=TRUE)[1]
)))
counts_path <- file.path(
  out_dir,
  "scrublet_negative_pseudobulk_counts.csv.gz"
)
annotation_path <- file.path(
  out_dir,
  "scrublet_negative_gene_annotation.csv"
)
meta_path <- file.path(out_dir, "scrublet_negative_sample_meta.csv")

counts_raw <- read_csv(counts_path, show_col_types=FALSE)
gene_ids <- counts_raw[[1]]
counts <- as.data.frame(counts_raw[, -1])
rownames(counts) <- gene_ids
counts <- round(as.matrix(counts))

annotation <- read_csv(annotation_path, show_col_types=FALSE) |>
  distinct(gene, .keep_all=TRUE)
meta <- read_csv(meta_path, show_col_types=FALSE) |>
  as.data.frame()
rownames(meta) <- meta$sample

if (!identical(colnames(counts), meta$sample)) {
  stop("Count columns do not exactly match sample metadata order")
}
if (any(rowSums(counts) < 10)) {
  stop("The frozen full-cohort gene universe contains a feature with <10 UMI")
}
if (!identical(sort(unique(meta$group)), c("Control", "Surgery"))) {
  stop("Unexpected group labels")
}

# Fail loudly if the explicitly requested native-mouse collection is unavailable
# or if its identity/version does not match the frozen Methods specification.
hallmark_raw <- msigdbr(
  db_species="MM",
  species="Mus musculus",
  collection="MH"
)
db_versions <- sort(unique(hallmark_raw$db_version))
if (!identical(db_versions, "2026.1.Mm")) {
  stop(
    paste(
      "Expected MSigDB 2026.1.Mm, obtained",
      paste(db_versions, collapse=", ")
    )
  )
}
hallmark <- hallmark_raw |>
  select(gs_name, gene_symbol) |>
  distinct()
if (n_distinct(hallmark$gs_name) != 50) {
  stop(
    paste(
      "Expected 50 mouse-native Hallmark sets, obtained",
      n_distinct(hallmark$gs_name)
    )
  )
}

focus <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
core <- c("Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3")

cohens_d <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  pooled_sd <- sqrt(
    ((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
      (length(x) + length(y) - 2)
  )
  if (is.finite(pooled_sd) && pooled_sd > 0) {
    (mean(x) - mean(y)) / pooled_sd
  } else {
    NA_real_
  }
}

exact_perm_p <- function(x, y) {
  pooled <- c(as.numeric(x), as.numeric(y))
  observed <- abs(mean(x) - mean(y))
  assignments <- combn(seq_along(pooled), length(x))
  permuted <- apply(assignments, 2, function(i) {
    abs(mean(pooled[i]) - mean(pooled[-i]))
  })
  mean(permuted >= observed - 1e-12)
}

fit_rows <- list()
gsea_rows <- list()

for (drop in c("(none)", meta$sample)) {
  keep <- if (drop == "(none)") meta$sample else setdiff(meta$sample, drop)
  local_meta <- meta[keep, , drop=FALSE]
  local_counts <- counts[, keep, drop=FALSE]
  col_data <- data.frame(
    group=factor(local_meta$group, levels=c("Control", "Surgery")),
    row.names=local_meta$sample
  )
  dds <- DESeqDataSetFromMatrix(
    countData=local_counts,
    colData=col_data,
    design=~group
  )
  dds <- DESeq(dds, quiet=TRUE)
  fit <- results(
    dds,
    contrast=c("group", "Surgery", "Control")
  ) |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    as_tibble() |>
    left_join(annotation, by="gene") |>
    mutate(
      dropped_animal=drop,
      n_surgery=sum(local_meta$group == "Surgery"),
      n_control=sum(local_meta$group == "Control"),
      ci95_low=log2FoldChange - 1.96 * lfcSE,
      ci95_high=log2FoldChange + 1.96 * lfcSE
    )
  fit_rows[[drop]] <- fit

  ranks <- fit |>
    filter(!is.na(symbol), !is.na(stat)) |>
    group_by(symbol) |>
    slice_max(order_by=abs(stat), n=1, with_ties=FALSE) |>
    ungroup() |>
    arrange(desc(stat))
  gene_list <- setNames(ranks$stat, ranks$symbol)
  set.seed(42)
  pathway_fit <- GSEA(
    gene_list,
    TERM2GENE=hallmark,
    pvalueCutoff=1,
    pAdjustMethod="BH",
    minGSSize=10,
    maxGSSize=500,
    eps=1e-30,
    seed=TRUE,
    verbose=FALSE
  ) |>
    as.data.frame() |>
    mutate(
      dropped_animal=drop,
      n_surgery=sum(local_meta$group == "Surgery"),
      n_control=sum(local_meta$group == "Control")
    )
  if (nrow(pathway_fit) != 50) {
    stop(
      paste(
        "Expected 50 tested Hallmark sets for drop =",
        drop,
        "but obtained",
        nrow(pathway_fit)
      )
    )
  }
  gsea_rows[[drop]] <- pathway_fit
}

fits <- bind_rows(fit_rows)
gsea <- bind_rows(gsea_rows)

full_fit <- fits |>
  filter(dropped_animal == "(none)") |>
  select(
    gene,
    symbol,
    baseMean,
    log2FoldChange,
    lfcSE,
    ci95_low,
    ci95_high,
    stat,
    pvalue,
    padj
  ) |>
  arrange(pvalue)
write_csv(
  full_fit,
  file.path(out_dir, "scrublet_negative_deseq2_full.csv")
)

# Preserve every transcriptome-wide LOO coefficient, rather than only the
# selected core genes, so the systematic sensitivity can be re-audited without
# rerunning DESeq2.
deseq_loo <- fits |>
  filter(dropped_animal != "(none)") |>
  select(
    dropped_animal,
    n_surgery,
    n_control,
    gene,
    symbol,
    baseMean,
    log2FoldChange,
    lfcSE,
    ci95_low,
    ci95_high,
    stat,
    pvalue,
    padj
  ) |>
  arrange(dropped_animal, pvalue)
write_csv(
  deseq_loo,
  file.path(out_dir, "scrublet_negative_deseq2_loo.csv.gz")
)

# Transcriptome-wide comparison with the all-7461-cell primary DESeq2 fit.
allcell_full_path <- file.path(
  project,
  "data/results/GSE267933_pseudobulk/deseq2_all_microglia_annotated.csv"
)
allcell_full <- read_csv(allcell_full_path, show_col_types=FALSE) |>
  select(
    gene,
    allcell_symbol=symbol,
    allcell_baseMean=baseMean,
    allcell_log2FoldChange=log2FoldChange,
    allcell_lfcSE=lfcSE,
    allcell_stat=stat,
    allcell_pvalue=pvalue,
    allcell_padj=padj
  )
deseq_compare <- full_join(
  full_fit |>
    rename(
      scrublet_negative_symbol=symbol,
      scrublet_negative_baseMean=baseMean,
      scrublet_negative_log2FoldChange=log2FoldChange,
      scrublet_negative_lfcSE=lfcSE,
      scrublet_negative_stat=stat,
      scrublet_negative_pvalue=pvalue,
      scrublet_negative_padj=padj
    ),
  allcell_full,
  by="gene"
) |>
  mutate(
    in_scrublet_negative_universe=
      !is.na(scrublet_negative_log2FoldChange),
    in_allcell_universe=!is.na(allcell_log2FoldChange),
    delta_log2FoldChange=
      scrublet_negative_log2FoldChange - allcell_log2FoldChange,
    sign_same=if_else(
      in_scrublet_negative_universe & in_allcell_universe,
      sign(scrublet_negative_log2FoldChange) ==
        sign(allcell_log2FoldChange),
      NA
    ),
    fdr05_same=if_else(
      in_scrublet_negative_universe & in_allcell_universe &
        !is.na(scrublet_negative_padj) & !is.na(allcell_padj),
      (scrublet_negative_padj < 0.05) == (allcell_padj < 0.05),
      NA
    )
  )
write_csv(
  deseq_compare,
  file.path(out_dir, "scrublet_negative_vs_allcell_deseq2.csv")
)

core_loo <- fits |>
  filter(symbol %in% core) |>
  select(
    dropped_animal,
    n_surgery,
    n_control,
    gene,
    symbol,
    baseMean,
    log2FoldChange,
    lfcSE,
    ci95_low,
    ci95_high,
    stat,
    pvalue,
    padj
  ) |>
  arrange(dropped_animal, match(symbol, core))
write_csv(
  core_loo,
  file.path(out_dir, "scrublet_negative_core_genes_loo.csv")
)

# Direct equal-gene-weight selected-panel summary from the Scrublet-negative
# pseudobulk counts. Total retained UMI (including features below the DE filter)
# comes from the extraction audit, matching the direct-expression estimator.
count_audit <- read_csv(
  file.path(out_dir, "scrublet_negative_count_audit.csv"),
  show_col_types=FALSE
)
core_counts <- as_tibble(counts, rownames="gene") |>
  left_join(annotation, by="gene") |>
  filter(symbol %in% core) |>
  group_by(symbol) |>
  summarise(across(all_of(meta$sample), sum), .groups="drop")
if (!setequal(core_counts$symbol, core)) {
  stop("Not all seven selected-panel genes were recovered")
}
core_values <- core_counts |>
  pivot_longer(
    cols=all_of(meta$sample),
    names_to="sample",
    values_to="gene_umi"
  ) |>
  left_join(
    count_audit |>
      select(sample, group, total_retained_umi=umi_retained),
    by="sample"
  ) |>
  mutate(
    log2_cpm=log2(
      (gene_umi + 0.5) / (total_retained_umi + 1) * 1e6
    )
  ) |>
  arrange(sample, match(symbol, core))
write_csv(
  core_values,
  file.path(out_dir, "scrublet_negative_selected_panel_animal_values.csv")
)

panel_values <- core_values |>
  group_by(sample, group) |>
  summarise(panel_mean_log2_cpm=mean(log2_cpm), .groups="drop")
panel_effect_rows <- list()
for (drop in c("(none)", meta$sample)) {
  local <- if (drop == "(none)") {
    panel_values
  } else {
    panel_values |> filter(sample != drop)
  }
  x <- local$panel_mean_log2_cpm[local$group == "Surgery"]
  y <- local$panel_mean_log2_cpm[local$group == "Control"]
  panel_effect_rows[[drop]] <- tibble(
    dropped_animal=drop,
    n_surgery=length(x),
    n_control=length(y),
    surgery_mean=mean(x),
    control_mean=mean(y),
    mean_difference=mean(x) - mean(y),
    cohens_d=cohens_d(x, y),
    exact_perm_p=if (
      drop == "(none)"
    ) exact_perm_p(x, y) else NA_real_
  )
}
panel_effects <- bind_rows(panel_effect_rows)
write_csv(
  panel_effects,
  file.path(out_dir, "scrublet_negative_selected_panel_effect_loo.csv")
)

allcell_panel_full <- read_csv(
  file.path(out_dir, "panel_effects.csv"),
  show_col_types=FALSE
) |>
  filter(
    mask == "all_microglia_labels",
    metric == "pseudobulk_log2cpm",
    panel == "core_isg7"
  ) |>
  transmute(
    dropped_animal="(none)",
    allcell_mean_difference=mean_difference,
    allcell_cohens_d=cohens_d,
    allcell_exact_perm_p=exact_perm_p
  )
allcell_panel_loo <- read_csv(
  file.path(out_dir, "panel_leave_one_animal_out.csv"),
  show_col_types=FALSE
) |>
  filter(
    mask == "all_microglia_labels",
    metric == "pseudobulk_log2cpm",
    panel == "core_isg7"
  ) |>
  transmute(
    dropped_animal,
    allcell_mean_difference=mean_difference,
    allcell_cohens_d=cohens_d,
    allcell_exact_perm_p=NA_real_
  )
panel_compare <- panel_effects |>
  left_join(
    bind_rows(allcell_panel_full, allcell_panel_loo),
    by="dropped_animal"
  ) |>
  mutate(
    delta_mean_difference=
      mean_difference - allcell_mean_difference,
    delta_cohens_d=cohens_d - allcell_cohens_d,
    sign_same=
      sign(mean_difference) == sign(allcell_mean_difference)
  )
if (any(is.na(panel_compare$allcell_mean_difference))) {
  stop("Failed to match all-cell selected-panel sensitivity rows")
}
write_csv(
  panel_compare,
  file.path(out_dir, "scrublet_negative_vs_allcell_selected_panel.csv")
)

gsea_full <- gsea |>
  filter(dropped_animal == "(none)") |>
  select(
    ID,
    Description,
    setSize,
    enrichmentScore,
    NES,
    pvalue,
    p.adjust,
    qvalue,
    rank,
    leading_edge,
    core_enrichment
  ) |>
  arrange(pvalue)
write_csv(
  gsea_full,
  file.path(out_dir, "scrublet_negative_hallmark_full.csv")
)

gsea_full_and_loo <- gsea |>
  select(
    dropped_animal,
    n_surgery,
    n_control,
    ID,
    Description,
    setSize,
    enrichmentScore,
    NES,
    pvalue,
    p.adjust,
    qvalue,
    rank,
    leading_edge,
    core_enrichment
  ) |>
  arrange(dropped_animal, pvalue)
write_csv(
  gsea_full_and_loo,
  file.path(out_dir, "scrublet_negative_hallmark_full_and_loo.csv")
)

gsea_focus_loo <- gsea |>
  filter(ID %in% focus) |>
  select(
    dropped_animal,
    n_surgery,
    n_control,
    ID,
    setSize,
    enrichmentScore,
    NES,
    pvalue,
    p.adjust,
    qvalue,
    core_enrichment
  ) |>
  arrange(dropped_animal, match(ID, focus))
write_csv(
  gsea_focus_loo,
  file.path(out_dir, "scrublet_negative_hallmark_focus_loo.csv")
)

# Compare against the previously audited all-7461-cell primary analysis.
allcell_gsea_path <- file.path(
  out_dir,
  "primary_hallmark_database_compare_loo.csv"
)
if (!file.exists(allcell_gsea_path)) {
  stop("Missing harmonized all-cell Hallmark comparison artifact")
}
allcell_gsea <- read_csv(allcell_gsea_path, show_col_types=FALSE) |>
  filter(database == "mouse_MH_native") |>
  select(
    dropped_animal,
    ID,
    allcell_setSize=setSize,
    allcell_NES=NES,
    allcell_pvalue=pvalue,
    allcell_p_adjust=p.adjust
  )

hallmark_compare <- gsea_focus_loo |>
  left_join(allcell_gsea, by=c("dropped_animal", "ID")) |>
  mutate(
    delta_NES_scrublet_negative_minus_allcell=NES - allcell_NES,
    sign_same=sign(NES) == sign(allcell_NES),
    fdr05_same=(p.adjust < 0.05) == (allcell_p_adjust < 0.05)
  )
if (any(is.na(hallmark_compare$allcell_NES))) {
  stop("Failed to match all-cell Hallmark sensitivity rows")
}
write_csv(
  hallmark_compare,
  file.path(out_dir, "scrublet_negative_vs_allcell_hallmark.csv")
)

allcell_core_path <- file.path(
  out_dir,
  "primary_pseudobulk_loo_core_genes.csv"
)
if (!file.exists(allcell_core_path)) {
  stop("Missing all-cell core-gene LOO artifact")
}
allcell_core <- read_csv(allcell_core_path, show_col_types=FALSE) |>
  select(
    dropped_animal,
    symbol,
    allcell_log2FoldChange=log2FoldChange,
    allcell_lfcSE=lfcSE,
    allcell_stat=stat,
    allcell_pvalue=pvalue,
    allcell_padj=padj
  )
core_compare <- core_loo |>
  left_join(allcell_core, by=c("dropped_animal", "symbol")) |>
  mutate(
    delta_log2FoldChange_scrublet_negative_minus_allcell=
      log2FoldChange - allcell_log2FoldChange,
    sign_same=
      sign(log2FoldChange) == sign(allcell_log2FoldChange)
  )
if (any(is.na(core_compare$allcell_log2FoldChange))) {
  stop("Failed to match all-cell core-gene sensitivity rows")
}
write_csv(
  core_compare,
  file.path(out_dir, "scrublet_negative_vs_allcell_core_genes.csv")
)

analysis_manifest <- tibble(
  item=c(
    "R",
    "DESeq2",
    "clusterProfiler",
    "msigdbr",
    "fgsea",
    "MSigDB",
    "Hallmark collection",
    "Hallmark sets",
    "full-cohort filtered features",
    "LOO gene-universe policy",
    "ranking metric",
    "random seed"
  ),
  value=c(
    as.character(getRversion()),
    as.character(packageVersion("DESeq2")),
    as.character(packageVersion("clusterProfiler")),
    as.character(packageVersion("msigdbr")),
    as.character(packageVersion("fgsea")),
    db_versions,
    "MH (mouse-native; no fallback)",
    as.character(n_distinct(hallmark$gs_name)),
    as.character(nrow(counts)),
    "retain full-cohort filtered features in every LOO fit",
    "signed DESeq2 Wald statistic",
    "42"
  )
)
write_csv(
  analysis_manifest,
  file.path(out_dir, "scrublet_negative_analysis_manifest.csv")
)

gate_summary <- tribble(
  ~gate, ~status, ~evidence,
  "Input", "PASS",
  "90/7461 predicted doublets removed; 7371 barcodes matched to integer UMI counts",
  "Replication", "PARTIAL",
  "Animal/library is the unit, but only 3 animals per group",
  "Doublet sensitivity", "PASS",
  paste0(
    "Across 4 focus Hallmarks x 7 full/LOO fits: ",
    sum(!hallmark_compare$sign_same),
    " sign and ",
    sum(!hallmark_compare$fdr05_same),
    " FDR<0.05 classification discordances versus all-cell"
  ),
  "Signal stability", "PARTIAL",
  paste0(
    "Full IFN enrichment persists, but drop-S3 remains ",
    "negative for IFN-alpha and IFN-gamma"
  ),
  "Primary-analysis suitability", "PASS_WITH_CAVEAT",
  paste0(
    "Scrublet-negative pseudobulk can be primary with all-cell sensitivity; ",
    "it does not rescue S3 dependence or support a coherent robust IFN program"
  )
)
write_csv(
  gate_summary,
  file.path(out_dir, "scrublet_negative_gate_summary.csv")
)

cat("Scrublet-negative full-cohort Hallmark focus\n")
print(
  as.data.frame(
    gsea_focus_loo |>
      filter(dropped_animal == "(none)") |>
      select(ID, setSize, NES, p.adjust)
  ),
  row.names=FALSE,
  digits=5
)
cat("\nScrublet-negative drop-S3 Hallmark focus\n")
print(
  as.data.frame(
    gsea_focus_loo |>
      filter(dropped_animal == "S3") |>
      select(ID, setSize, NES, p.adjust)
  ),
  row.names=FALSE,
  digits=5
)
cat("\nScrublet-negative versus all-cell qualitative Hallmark differences\n")
print(
  as.data.frame(
    hallmark_compare |>
      group_by(ID) |>
      summarise(
        sign_discordant=sum(!sign_same),
        fdr05_discordant=sum(!fdr05_same),
        max_abs_delta_NES=max(
          abs(delta_NES_scrublet_negative_minus_allcell)
        ),
        .groups="drop"
      )
  ),
  row.names=FALSE,
  digits=5
)
cat("\nDirect selected-panel effect\n")
print(as.data.frame(panel_compare), row.names=FALSE, digits=5)
