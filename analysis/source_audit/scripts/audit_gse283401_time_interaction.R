#!/usr/bin/env Rscript

# Formal old-animal time x treatment model for GSE283401.
# Positive interaction = the surgery-control contrast is more positive at 48 h
# than at 6 h. Time remains potentially confounded with sequencing batch.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

project <- Sys.getenv("POCD_SCRNA_PROJECT_ROOT")
if (!nzchar(project)) stop("Set POCD_SCRNA_PROJECT_ROOT before running this audit.")
project <- normalizePath(project, mustWork = TRUE)
data_dir <- file.path(project, "data/raw/GSE283401")
out_dir <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])))

meta <- tribble(
  ~sample, ~time, ~treatment,
  "FACS40_HP", "6h", "control",
  "FACS43_HP", "6h", "surgery",
  "FACS44_HP", "6h", "control",
  "FACS47_HP", "6h", "surgery",
  "FACS48_HP", "6h", "control",
  "FACS51_HP", "6h", "surgery",
  "FACS52_HP", "6h", "control",
  "FACS55_HP", "6h", "surgery",
  "FACS56_HP", "6h", "control",
  "FACS59_HP", "6h", "surgery",
  "FACS60_HP", "6h", "control",
  "FACS63_HP", "6h", "surgery",
  "FACS64_HP", "6h", "control",
  "FACS71_HP", "6h", "surgery",
  "FACS72_HP", "6h", "control",
  "FACS102_HP", "48h", "control",
  "FACS105_HP", "48h", "surgery",
  "FACS106_HP", "48h", "control",
  "FACS107_HP", "48h", "surgery",
  "FACS109_HP", "48h", "surgery",
  "FACS76_HP", "48h", "control",
  "FACS77_HP", "48h", "surgery",
  "FACS80_HP", "48h", "control",
  "FACS83_HP", "48h", "surgery",
  "FACS86_HP", "48h", "control",
  "FACS88_HP", "48h", "control",
  "FACS89_HP", "48h", "surgery",
  "FACS93_HP", "48h", "surgery",
  "FACS95_HP", "48h", "surgery",
  "FACS96_HP", "48h", "control",
  "FACS98_HP", "48h", "control",
  "FACS99_HP", "48h", "surgery"
) |>
  mutate(
    time=factor(time, levels=c("6h", "48h")),
    treatment=factor(treatment, levels=c("control", "surgery"))
  ) |>
  as.data.frame()
rownames(meta) <- meta$sample

read_counts <- function(filename) {
  d <- read_csv(file.path(data_dir, filename), show_col_types=FALSE)
  ids <- d$gene
  x <- as.data.frame(d[, -which(names(d) == "gene")])
  rownames(x) <- ids
  x
}

c6 <- read_counts("GSE283401_hp_6h_counts.csv.gz")
c48 <- read_counts("GSE283401_hp_48h_counts.csv.gz")
common <- intersect(rownames(c6), rownames(c48))
counts <- cbind(
  c6[common, meta$sample[meta$time == "6h"], drop=FALSE],
  c48[common, meta$sample[meta$time == "48h"], drop=FALSE]
)
counts <- round(as.matrix(counts[, meta$sample, drop=FALSE]))

dds <- DESeqDataSetFromMatrix(counts, meta, design=~time + treatment + time:treatment)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, quiet=TRUE)
interaction_name <- grep("time48h.*treatmentsurgery", resultsNames(dds), value=TRUE)
if (length(interaction_name) != 1) {
  stop(paste("Could not identify interaction coefficient:", paste(resultsNames(dds), collapse=", ")))
}
fit <- results(dds, name=interaction_name) |>
  as.data.frame() |>
  rownames_to_column("ensembl") |>
  as_tibble()

clean_ids <- sub("\\.\\d+$", "", fit$ensembl)
fit$symbol <- unname(
  mapIds(
    org.Mm.eg.db,
    keys=clean_ids,
    column="SYMBOL",
    keytype="ENSEMBL",
    multiVals="first"
  )
)
write_csv(fit, file.path(out_dir, "gse283401_old_time_by_treatment_deseq2.csv"))

hallmark <- tryCatch(
  msigdbr(db_species="MM", species="Mus musculus", collection="MH") |>
    dplyr::select(gs_name, gene_symbol),
  error=function(e) {
    msigdbr(species="Mus musculus", collection="H") |>
      dplyr::select(gs_name, gene_symbol)
  }
)

ranks <- fit |>
  filter(!is.na(symbol), !is.na(stat)) |>
  group_by(symbol) |>
  slice_max(abs(stat), n=1, with_ties=FALSE) |>
  ungroup() |>
  arrange(desc(stat))
gene_list <- setNames(ranks$stat, ranks$symbol)
set.seed(42)
gsea <- GSEA(
  gene_list,
  TERM2GENE=hallmark,
  pvalueCutoff=1,
  minGSSize=10,
  maxGSSize=500,
  eps=1e-30,
  seed=TRUE,
  verbose=FALSE
) |>
  as.data.frame()
write_csv(gsea, file.path(out_dir, "gse283401_old_time_by_treatment_hallmark.csv"))

focus <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)
cat("Interaction coefficient:", interaction_name, "\n")
print(
  as.data.frame(
    gsea |>
      filter(ID %in% focus) |>
      dplyr::select(ID, setSize, NES, pvalue, p.adjust)
  ),
  row.names=FALSE,
  digits=4
)

core <- c("Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3")
cat("\nCore ISG interaction effects\n")
print(
  as.data.frame(
    fit |>
      filter(symbol %in% core) |>
      dplyr::select(symbol, log2FoldChange, lfcSE, stat, pvalue, padj)
  ),
  row.names=FALSE,
  digits=4
)
