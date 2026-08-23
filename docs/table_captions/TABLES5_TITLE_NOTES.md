# Table S5 title and notes

## Title

**Table S5. GSE283401 hippocampal sample provenance, old-animal treatment contrasts, and gene-level and Hallmark time-by-treatment analyses.**

## General note

The hippocampal component of GSE283401 contains 48 bulk RNA-sequencing libraries from FACS-isolated microglia of male C57BL/6 mice. The source design included young mice aged 3–5 months and old mice aged 20–22 months. Exposed animals received 1.2% isoflurane in 30% oxygen balanced with nitrogen for 2 h plus laparotomy. Controls received carrier gas for 2 h and the same bupivacaine and meloxicam regimen, but no anesthesia or laparotomy. The 6- and 48-h endpoints were measured from the start of anesthesia or carrier-gas administration.

Formal analyses in Parts B–E were restricted to 32 old-animal libraries: eight carrier-gas controls and seven exposed animals at 6 h, and eight controls and nine exposed animals at 48 h. Each animal contributed one library and constituted one biological unit. Different animals were sampled at the two endpoints. The accession also contains 48 hemisphere libraries, which were not included in the hippocampal manifest or analyses reported here.

The 6- and 48-h counts were deposited as separate matrices. Sampling time is therefore aligned with deposited matrix and with any unrecorded processing or sequencing differences between the experiments. These analyses provide external cross-experiment context but do not constitute a within-animal time course, replication of GSE267933, or orthogonal validation. The deposited treatment label `surgery` is retained in machine-readable provenance fields; manuscript-facing fields describe the treatment as isoflurane plus laparotomy.

## Part A — GSE283401 sample manifest

Part A reports GEO, BioSample, and SRA identifiers; source metadata; matrix order; assay information; deposited matrix; total library counts; and analysis-inclusion status for all 48 hippocampal libraries. Young animals document the hippocampal source design but were not included in the old-animal models. `included_in_old_model = TRUE` identifies the 32 libraries used in Parts B–E. `library_total` is the sum of deposited counts for that library and is descriptive rather than a normalization factor or sample-quality exclusion rule.

## Part B — Old-animal 6-h gene-level treatment contrast

Part B reports the complete DESeq2 `~ treatment` results for the 15 old animals sampled at 6 h. The feature universe contains 26,317 Ensembl rows with at least 10 total counts across these libraries. Positive `log2FoldChange` values denote higher expression under isoflurane plus laparotomy than under carrier-gas control, and negative values denote the reverse. Log2 fold changes are unshrunk; `ci95_low` and `ci95_high` are 95% Wald confidence limits. `padj` is the Benjamini–Hochberg adjusted *P* value returned after DESeq2 independent filtering. An unavailable adjusted value is reported as `NA` and is not evidence of significance or absence of an effect.

Seventy-three Ensembl features met adjusted *P* < 0.05: 40 had positive and 33 had negative exposure-minus-control log2 fold changes. These counts describe the separately deposited 6-h experiment.

## Part C — Old-animal 48-h gene-level treatment contrast

Part C reports the complete DESeq2 `~ treatment` results for the 17 old animals sampled at 48 h. The feature universe contains 23,732 Ensembl rows with at least 10 total counts across these libraries. Contrast orientation, unshrunk log2 fold changes, Wald confidence intervals, and adjusted-*P*-value handling follow Part B. Parts B and C are separate cross-sectional fits and are not paired measurements.

No Ensembl feature met adjusted *P* < 0.05 in this fit. The 73-versus-0 comparison between Parts B and C involves different animals, matrices, feature universes, and multiple-testing families and does not establish disappearance, emergence, or resolution over time.

## Part D — Old-animal gene-level time-by-treatment interaction

Part D reports the complete gene-level results from `~ time + treatment + time:treatment` across all 32 old animals. The feature universe contains 30,280 shared Ensembl rows with at least 10 total counts across the combined libraries. The reported interaction log2 fold change is `(48-h exposed − 48-h control) − (6-h exposed − 6-h control)`. A positive value denotes a more positive treatment-associated contrast in the deposited 48-h experiment than in the deposited 6-h experiment; it does not denote within-animal induction or establish a temporal transition. Gene-level confidence intervals and adjusted *P* values follow Parts B–C.

Twelve Ensembl features met adjusted *P* < 0.05 in the interaction fit: five had positive and seven had negative coefficients. The signs describe between-experiment differences in treatment contrast rather than within-animal induction or suppression.

For Parts B–D, version suffixes were removed from Ensembl identifiers before mapping to mouse gene symbols with `org.Mm.eg.db` 3.22.0 and `multiVals = "first"`. `baseMean` is the DESeq2 mean of normalized counts; `lfcSE` is the standard error of the unshrunk log2 fold change; `stat` is the Wald statistic; `pvalue` is the unadjusted Wald-test *P* value; and `padj` is the Benjamini–Hochberg adjusted *P* value. Adjustment was performed separately within each model after independent filtering. Parts B, C, and D contain 12,525, 23,720, and 6,796 nonmissing adjusted *P* values, respectively.

## Part E — Mouse MSigDB Hallmark time-by-treatment interaction

Part E reports all 50 mouse MSigDB 2026.1.Mm Hallmark gene sets. Genes from Part D were ranked by the signed interaction Wald statistic. When multiple Ensembl rows mapped to one gene symbol, the row with the largest absolute Wald statistic was retained. Enrichment was calculated with `clusterProfiler::GSEA` using `minGSSize = 10`, `maxGSSize = 500`, `pvalueCutoff = 1`, `eps = 1 × 10^−30`, and seed 42. No human-gene-set fallback was used.

`NES` is the normalized enrichment score; positive values denote ranking toward a more positive treatment contrast in the deposited 48-h experiment than in the deposited 6-h experiment. `p.adjust` is the Benjamini–Hochberg false-discovery rate across all 50 sets, not only those displayed in Figure 4. NES is a rank-based transcript-level statistic and does not measure protein abundance, transcription-factor activity, cytokine signaling, or pathway activity.

The interaction ranking contained 21,514 unique mapped mouse symbols after duplicate-symbol handling. `enrichmentScore` is the unnormalized enrichment score; `pvalue`, `p.adjust`, and `qvalue` are values returned by `clusterProfiler::GSEA`; `rank` is the ranked-list position at which the enrichment score was attained; `leading_edge` summarizes tag, list, and signal percentages; and `core_enrichment` lists leading-edge genes. `nes_rank_descending` and `absolute_nes_rank` rank the 50 sets by NES and absolute NES, respectively. `fdr_lt_0_05` indicates whether `p.adjust` is below 0.05. Eleven Hallmark sets met this threshold.

## Data dictionary

The accompanying data dictionary describes the contents, unit or scale, and interpretation limit for Parts A–E.
