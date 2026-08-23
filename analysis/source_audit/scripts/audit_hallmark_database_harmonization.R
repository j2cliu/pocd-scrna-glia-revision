#!/usr/bin/env Rscript

# Compare the human Hallmark collection ortholog-mapped to mouse with the
# mouse-native Hallmark collection, holding the primary GSE267933 DESeq2
# ranking and all GSEA settings fixed. This is an audit artifact only.

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

human_raw <- msigdbr(
  db_species="HS",
  species="Mus musculus",
  collection="H"
)
mouse_raw <- msigdbr(
  db_species="MM",
  species="Mus musculus",
  collection="MH"
)

hallmark <- list(
  human_H_ortholog_mapped=human_raw |>
    select(gs_name, gene_symbol) |>
    distinct(),
  mouse_MH_native=mouse_raw |>
    select(gs_name, gene_symbol) |>
    distinct()
)

focus <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

membership_rows <- lapply(focus, function(pathway) {
  human_genes <- sort(unique(
    hallmark$human_H_ortholog_mapped$gene_symbol[
      hallmark$human_H_ortholog_mapped$gs_name == pathway
    ]
  ))
  mouse_genes <- sort(unique(
    hallmark$mouse_MH_native$gene_symbol[
      hallmark$mouse_MH_native$gs_name == pathway
    ]
  ))
  shared <- intersect(human_genes, mouse_genes)
  union_genes <- union(human_genes, mouse_genes)
  tibble(
    ID=pathway,
    human_H_n=length(human_genes),
    mouse_MH_n=length(mouse_genes),
    shared_n=length(shared),
    union_n=length(union_genes),
    jaccard=length(shared) / length(union_genes),
    human_H_only=paste(setdiff(human_genes, mouse_genes), collapse=";"),
    mouse_MH_only=paste(setdiff(mouse_genes, human_genes), collapse=";")
  )
})
membership <- bind_rows(membership_rows)
write_csv(
  membership,
  file.path(out_dir, "hallmark_database_focus_membership.csv")
)

db_manifest <- tibble(
  database=c("human_H_ortholog_mapped", "mouse_MH_native"),
  db_species=c("HS", "MM"),
  target_species=c("Mus musculus", "Mus musculus"),
  collection=c("H", "MH"),
  msigdb_version=c(
    paste(sort(unique(human_raw$db_version)), collapse=";"),
    paste(sort(unique(mouse_raw$db_version)), collapse=";")
  ),
  rows=c(nrow(hallmark[[1]]), nrow(hallmark[[2]])),
  gene_sets=c(
    n_distinct(hallmark[[1]]$gs_name),
    n_distinct(hallmark[[2]]$gs_name)
  ),
  unique_genes=c(
    n_distinct(hallmark[[1]]$gene_symbol),
    n_distinct(hallmark[[2]]$gene_symbol)
  ),
  msigdbr_version=as.character(packageVersion("msigdbr"))
)
write_csv(
  db_manifest,
  file.path(out_dir, "hallmark_database_manifest.csv")
)

gsea_rows <- list()

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
    filter(!is.na(symbol), !is.na(stat))

  ranks <- fit |>
    group_by(symbol) |>
    slice_max(order_by=abs(stat), n=1, with_ties=FALSE) |>
    ungroup() |>
    arrange(desc(stat))
  gene_list <- setNames(ranks$stat, ranks$symbol)

  for (database in names(hallmark)) {
    set.seed(42)
    pathway_fit <- GSEA(
      gene_list,
      TERM2GENE=hallmark[[database]],
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
        database,
        ID,
        setSize,
        NES,
        pvalue,
        p.adjust,
        core_enrichment
      )
    gsea_rows[[paste(drop, database, sep="::")]] <- pathway_fit
  }
}

gsea <- bind_rows(gsea_rows)
write_csv(
  gsea,
  file.path(out_dir, "primary_hallmark_database_compare_loo.csv")
)

paired <- gsea |>
  select(dropped_animal, database, ID, setSize, NES, p.adjust) |>
  pivot_wider(
    names_from=database,
    values_from=c(setSize, NES, p.adjust)
  ) |>
  mutate(
    delta_NES_mouse_minus_human=
      NES_mouse_MH_native - NES_human_H_ortholog_mapped,
    sign_same=
      sign(NES_mouse_MH_native) == sign(NES_human_H_ortholog_mapped),
    fdr05_same=
      (p.adjust_mouse_MH_native < 0.05) ==
      (p.adjust_human_H_ortholog_mapped < 0.05)
  )
write_csv(
  paired,
  file.path(out_dir, "primary_hallmark_database_compare_loo_paired.csv")
)

cat("Database manifest\n")
print(as.data.frame(db_manifest), row.names=FALSE, digits=4)
cat("\nFocus-set membership comparison\n")
print(
  as.data.frame(
    membership |>
      select(ID, human_H_n, mouse_MH_n, shared_n, union_n, jaccard)
  ),
  row.names=FALSE,
  digits=4
)
cat("\nPrimary GSE267933: full cohort and drop-S3 database comparison\n")
print(
  as.data.frame(
    gsea |>
      filter(
        dropped_animal %in% c("(none)", "S3"),
        ID %in% c(
          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
          "HALLMARK_INTERFERON_GAMMA_RESPONSE"
        )
      ) |>
      select(dropped_animal, database, ID, setSize, NES, p.adjust)
  ),
  row.names=FALSE,
  digits=5
)
cat("\nQualitative concordance across all leave-one-animal-out fits\n")
print(
  as.data.frame(
    paired |>
      group_by(ID) |>
      summarise(
        comparisons=n(),
        sign_discordant=sum(!sign_same),
        fdr05_discordant=sum(!fdr05_same),
        max_abs_delta_NES=max(abs(delta_NES_mouse_minus_human)),
        .groups="drop"
      )
  ),
  row.names=FALSE,
  digits=4
)

