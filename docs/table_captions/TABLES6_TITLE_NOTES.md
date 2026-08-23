# Table S6 title and notes

## Title

**Table S6. GSE289098 same-cohort processed-count correspondence and animal-level sensitivity analyses on the fixed 7,371-cell microglial population.**

## Parts

- **Part A, global payload and analysis summary:** matrix dimensions; one-to-one cell-barcode and feature-identifier correspondence between the primary GSE267933 and alternative GSE289098 processed-count matrices; matrix-wide and fixed-analysis-population count differences; common-estimand checks; and fixed feature universes.
- **Part B, library mapping and counts:** animal/library identifiers, GSE267933 GSM/BioSample/SRA Experiment accessions, GSE289098 barcode-suffix mapping, all-cell and microglial cell counts, and per-library count differences.
- **Part C, paired animal values:** transcript UMI counts, all-feature UMI denominators, transcript-level log2 counts-per-million values, equal-weight seven-transcript values, and GSE289098-minus-primary paired differences for each animal.
- **Part D, paired full-cohort direct effects:** group means, combined-exposure-minus-oxygen-control differences, 95% Welch–Satterthwaite confidence intervals, Cohen’s *d*, and exact permutation *P* values for the seven transcripts and their equal-weight panel.
- **Part E, panel full-cohort and one-animal-omission results:** matched full-cohort and six systematic animal-omission panel estimates under both count matrices, including paired contrast differences and direction concordance.
- **Part F, selected-transcript DESeq2 results:** unshrunk DESeq2 estimates under both count matrices using the same primary-derived 13,926-feature universe.
- **Part G, DESeq2 size factors:** median-ratio size factors for each animal and count matrix.
- **Part H, transcript-level full-cohort and one-animal-omission results:** matched direct-expression effects for all seven transcripts. Forty-eight of 49 full-cohort or omission scenarios retained the same direction. The only boundary crossing was the near-zero *Isg15* estimate after S3 was omitted (primary, −0.031; GSE289098, +0.052 log2 counts per million).
- **Data dictionary:** row unit and purpose of each part.

## General notes

1. GSE289098 contains a Cell Ranger v3.0.2-aggregated processed-count matrix derived from the same six GSE267933 libraries and 20,684 filtered cell barcodes. It adds no animals and does not provide independent biological replication or orthogonal validation. Machine-readable `source_group` and `dropped_source_group` fields retain the deposited labels `Control` and `Surgery`; manuscript-facing text maps these labels to oxygen control and combined sevoflurane-plus-laparotomy exposure, respectively.
2. The analysis carries forward the exact 7,371 Scrublet-negative microglia by sample-aware barcode matching. No GSE289098-specific quality-control threshold, microglial reclassification, clustering, or subtype definition was applied. Cells remain nested within animals; the biological unit is the animal/library (*n* = 3 per group).
3. Primary counts were obtained from the verified integer-count object corresponding to the six original GSE267933 count matrices. Processed-object `.raw` attributes containing log-normalized expression were not used as counts. GSE289098 counts were obtained from the matched 27,998-feature processed matrix.
4. The seven transcripts selected during exploratory inspection were *Irf7*, *Ifitm3*, *Isg15*, *Mx1*, *Ifit1*, *Ifit2*, and *Ifit3*. For each transcript and animal, expression was calculated as `log2[((transcript UMI + 0.5)/(total UMI across all 27,998 features + 1)) × 10^6]`. The panel value is the arithmetic mean of the seven transcript-level values rather than a value calculated from summed panel counts.
5. Direct effects are combined-exposure-minus-oxygen-control differences in animal-level values. Confidence intervals use Welch–Satterthwaite degrees of freedom, and Cohen’s *d* uses the pooled within-group standard deviation. Descriptive full-cohort exact two-sided permutation *P* values enumerate all 20 possible three-versus-three label assignments without a +1 correction. No permutation *P* values are assigned to unequal-group omission rows.
6. Alternative-minus-primary quantities are deterministic paired descriptions of two processing outputs from the same libraries and were not assigned paired-test *P* values.
7. DESeq2 models used `~ group`, with combined exposure relative to oxygen control. The primary total-count threshold of at least 10 defines the common 13,926-feature universe for both fits; 13,831 features would independently meet that threshold in GSE289098. Wald confidence intervals equal log2 fold change ± 1.96 standard errors. Benjamini–Hochberg adjusted *P* values were calculated separately within each matrix after DESeq2 independent filtering. Because independent filtering is matrix-specific, adjusted *P* values are reported for completeness but are not used as an equivalence or concordance criterion.
8. Full-cohort and all six one-animal-omission panel estimates remained positive under both count matrices. The estimates changed little in this single same-cohort processed-count comparison, which does not establish biological replication or generalizability across cohorts.
