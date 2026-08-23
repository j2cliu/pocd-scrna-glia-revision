#!/usr/bin/env Rscript

# Refit the primary GSE267933 pseudobulk model after dropping each animal.
# This separates compact core-ISG directionality from a broad Hallmark program.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(msigdbr)
})

project <- Sys.getenv("POCD_SCRNA_PROJECT_ROOT")
if (!nzchar(project)) stop("Set POCD_SCRNA_PROJECT_ROOT before running this audit.")
project <- normalizePath(project, mustWork = TRUE)
out_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])))
pb_dir <- file.path(project, "data/results/GSE267933_pseudobulk")

counts_raw <- read_csv(
  file.path(pb_dir, "pseudobulk_all_microglia_counts.csv"),
  show_col_types=FALSE
)
gene_ids <- counts_raw[[1]]
counts <- as.data.frame(counts_raw[, -1])
rownames(counts) <- gene_ids

meta <- read_csv(
  file.path(pb_dir, "pseudobulk_all_microglia_meta.csv"),
  show_col_types=FALSE
) |>
  as.data.frame()
rownames(meta) <- meta$sample

annotation <- read_csv(
  file.path(pb_dir, "deseq2_all_microglia_annotated.csv"),
  show_col_types=FALSE
) |>
  select(gene, symbol) |>
  distinct(gene, .keep_all=TRUE)

hallmark <- tryCatch(
  msigdbr(species="Mus musculus", collection="H") |>
    select(gs_name, gene_symbol),
  error=function(e) {
    msigdbr(species="Mus musculus", category="H") |>
      select(gs_name, gene_symbol)
  }
)

focus <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
core <- c("Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3")

gsea_rows <- list()
core_rows <- list()

for (drop in c("(none)", meta$sample)) {
  keep <- if (drop == "(none)") meta$sample else setdiff(meta$sample, drop)
  local_meta <- meta[keep, , drop=FALSE]
  local_counts <- round(as.matrix(counts[, keep, drop=FALSE]))
  col_data <- data.frame(
    group=factor(local_meta$group, levels=c("Control", "Surgery")),
    row.names=local_meta$sample
  )
  dds <- DESeqDataSetFromMatrix(local_counts, col_data, design=~group)
  dds <- DESeq(dds, quiet=TRUE)
  fit <- results(dds, contrast=c("group", "Surgery", "Control")) |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    as_tibble() |>
    left_join(annotation, by="gene") |>
    filter(!is.na(symbol))

  core_rows[[drop]] <- fit |>
    filter(symbol %in% core) |>
    transmute(
      dropped_animal=drop,
      symbol,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    )

  ranks <- fit |>
    filter(!is.na(stat)) |>
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
    minGSSize=10,
    maxGSSize=500,
    eps=1e-30,
    seed=TRUE,
    verbose=FALSE
  ) |>
    as.data.frame() |>
    filter(ID %in% focus) |>
    transmute(
      dropped_animal=drop,
      ID,
      setSize,
      NES,
      pvalue,
      p.adjust,
      core_enrichment
    )
  gsea_rows[[drop]] <- pathway_fit
}

gsea <- bind_rows(gsea_rows)
core_result <- bind_rows(core_rows)
write_csv(gsea, file.path(out_dir, "primary_pseudobulk_loo_hallmark.csv"))
write_csv(core_result, file.path(out_dir, "primary_pseudobulk_loo_core_genes.csv"))

cat("Hallmark leave-one-animal-out\n")
print(
  as.data.frame(
    gsea |>
      filter(ID %in% c(
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE"
      )) |>
      select(dropped_animal, ID, NES, p.adjust)
  ),
  row.names=FALSE,
  digits=4
)

cat("\nCore ISG after dropping S3\n")
print(
  as.data.frame(
    core_result |>
      filter(dropped_animal == "S3") |>
      select(symbol, log2FoldChange, lfcSE, stat, pvalue, padj)
  ),
  row.names=FALSE,
  digits=4
)
