# Table S3 title and notes

## Title

**Table S3. Complete primary animal-level pseudobulk differential expression, animal-influence summaries, and detailed results for the selected seven-transcript panel.**

## General note

The seven interferon-responsive transcripts (*Irf7*, *Ifitm3*, *Isg15*, *Mx1*, *Ifit1*, *Ifit2*, and *Ifit3*) were selected during exploratory inspection of the primary dataset and used to summarize animal-level expression patterns. Primary analyses used the 7,371 Scrublet-negative microglia, with the animal/library as the biological unit. C1–C3 denote oxygen-control animals; S1–S3 denote animals receiving combined anesthesia-plus-laparotomy exposure. Parts D–E provide transcriptome-wide context for the selected-transcript results. Positive log2 fold changes denote higher expression under combined exposure, and negative values denote higher expression in oxygen controls.

## Part A — Animal values

`TableS3_animal_values` reports each transcript’s UMI count, retained microglial UMI denominator, and direct animal-level expression, calculated as `log2[((gene UMI + 0.5)/(total retained microglial UMI + 1)) × 10^6]`. The panel value is the arithmetic mean of the seven transcript-level values, giving each transcript equal weight.

## Part B — Full-cohort effects

`TableS3_direct_and_deseq2_effects` reports combined-exposure and oxygen-control means, the combined-exposure-minus-control mean difference, Welch–Satterthwaite 95% confidence interval, Cohen’s *d*, and two-sided exact label-permutation *P* value across all 20 three-versus-three assignments for each direct-expression outcome. These confidence intervals and exact permutation tests are descriptive and are not multiplicity-adjusted.

Transcript rows also report unshrunk DESeq2 log2 fold changes, 95% Wald intervals, Wald statistics, nominal *P* values, and Benjamini–Hochberg adjusted *P* values from animal-level pseudobulk UMI counts. Adjustment was performed across the complete 13,926-feature DESeq2 analysis, not only the seven displayed transcripts; Wald intervals are not multiplicity-adjusted. Adjusted *P* values are `NA` for five selected transcripts after DESeq2 independent filtering and were not replaced by nonsignificant numerical values. DESeq2 and direct-expression columns represent distinct estimands derived from the same six libraries.

## Part C — Animal-influence summaries

`TableS3_leave_one_animal_out` reports the full-cohort direct-expression estimate and six systematic one-animal-omission estimates for every transcript and for the equal-weight panel. Omission rows are influence summaries rather than independent hypothesis tests; exact permutation *P* values are therefore reported only for the full three-versus-three comparison. `sign_matches_full` indicates whether the raw mean-difference sign matches the corresponding full-cohort sign. `perfect_group_separation` is a descriptive ordering indicator for the retained animal values. No animal was excluded from the full-cohort estimate.

## Part D — Complete primary transcriptome-wide DESeq2 results

`TableS3D_primary_transcriptome_deseq2` reports all 13,926 Ensembl features tested in the primary `~ group` model. Fields include raw UMI counts for each animal/library, normalized base mean, unshrunk log2 fold change, Wald standard error and 95% interval, Wald statistic, nominal *P* value, and Benjamini–Hochberg adjusted *P* value. DESeq2 independent filtering left 6,632 nonmissing adjusted *P* values. Forty features met adjusted *P* < 0.05: 11 had positive and 29 had negative log2 fold changes. `bh_fdr_lt_0_05` is a reporting indicator rather than an additional test.

The final three columns summarize the six systematic one-animal-omission refits. They report the number of omission fits in which a feature met the false-discovery-rate threshold, the number that also retained the full-fit direction, and whether both conditions held in all six refits. These are influence summaries rather than independent replications. No animal was excluded from the full fit.

## Part E — Transcriptome-wide animal-omission summary

`TableS3E_primary_transcriptome_loo_summary` reports, for the full fit and each one-animal omission, the group sizes, tested-feature count, number of nonmissing adjusted *P* values, number of false-discovery-rate-threshold features in each direction, and number of the 40 full-fit hits retained. Threshold-hit counts ranged from 10 to 191 across omissions. Five full-fit hits—*Ccl9*, *H1f2*, *Lst1*, *Tnfaip3*, and *Ttr*—met adjusted *P* < 0.05 with the same direction in all six omission refits.

The raw-count detection results for *Ttr* are provided in Figure S1 and Table S1. Because unfiltered droplet matrices were unavailable, biological *Ttr* expression could not be separated from ambient RNA or other sample-processing contributions. Neither the false-discovery-rate-hit counts nor their omission patterns define a coherent exposure-induced program.
