#!/usr/bin/env Rscript

# Independent pathway-level audit of the primary GSE267933 pseudobulk result.
# The prespecified primary ranking is the DESeq2 signed Wald statistic. The
# manuscript's historical log2FC * -log10(p) rank is retained as a sensitivity.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
})

project <- Sys.getenv("POCD_SCRNA_PROJECT_ROOT")
if (!nzchar(project)) stop("Set POCD_SCRNA_PROJECT_ROOT before running this audit.")
project <- normalizePath(project, mustWork = TRUE)
out_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])))
input <- file.path(
  project,
  "data/results/GSE267933_pseudobulk/deseq2_all_microglia_annotated.csv"
)

d <- read_csv(input, show_col_types=FALSE) |>
  filter(!is.na(symbol), !is.na(stat), !is.na(log2FoldChange), !is.na(pvalue)) |>
  group_by(symbol) |>
  slice_max(order_by=abs(stat), n=1, with_ties=FALSE) |>
  ungroup()

hallmark <- tryCatch(
  msigdbr(species="Mus musculus", collection="H") |>
    select(gs_name, gene_symbol),
  error=function(e) {
    msigdbr(species="Mus musculus", category="H") |>
      select(gs_name, gene_symbol)
  }
)

run_gsea <- function(metric_name, values) {
  ranks <- setNames(values, d$symbol)
  ranks <- sort(ranks[is.finite(ranks)], decreasing=TRUE)
  set.seed(42)
  fit <- GSEA(
    ranks,
    TERM2GENE=hallmark,
    pvalueCutoff=1,
    minGSSize=10,
    maxGSSize=500,
    eps=1e-30,
    seed=TRUE,
    verbose=FALSE
  )
  result <- as.data.frame(fit) |>
    mutate(rank_metric=metric_name)
  write_csv(result, file.path(out_dir, paste0("primary_hallmark_gsea_", metric_name, ".csv")))
  result
}

wald <- run_gsea("signed_wald", d$stat)
legacy <- run_gsea(
  "lfc_times_neglog10p",
  d$log2FoldChange * -log10(pmax(d$pvalue, .Machine$double.xmin))
)

focus <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

combined <- bind_rows(wald, legacy) |>
  filter(ID %in% focus) |>
  select(rank_metric, ID, setSize, NES, pvalue, p.adjust, core_enrichment)
write_csv(combined, file.path(out_dir, "primary_hallmark_gsea_focus.csv"))

cat("Primary GSE267933 pseudobulk Hallmark GSEA\n")
print(as.data.frame(combined[, c("rank_metric", "ID", "setSize", "NES", "pvalue", "p.adjust")]),
      row.names=FALSE, digits=4)
cat("\nR session:\n")
print(sessionInfo())
