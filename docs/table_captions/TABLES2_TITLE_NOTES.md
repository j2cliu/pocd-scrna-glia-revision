# Table S2 title and notes

## Title

**Table S2. Reference-partition seed sensitivity and conditional inherited-PC graph-reconstruction analyses for GSE267933 microglia.**

## Parts

- **Part A, fixed stored-graph seed results:** partition count and adjusted Rand index relative to the reference partition for six Leiden seeds at resolution 0.4 on the stored neighborhood graph.
- **Part B, fixed stored-graph pairwise adjusted Rand indices:** all 15 unique pairwise adjusted Rand indices among the same six fixed-graph runs.
- **Part C, reference-partition subsampling Jaccard summaries:** per-partition best-match Jaccard results after 80% cell subsampling and graph reconstruction at resolution 0.4, with variable- and fixed-seed analyses reported separately.
- **Part D, subsampling partition-count distribution:** numbers and percentages of runs yielding five through nine partitions in the 100 variable-seed and 50 fixed-seed graph-reconstruction analyses.
- **Part E, resolution-sweep summary:** full-cohort and modal subsample partition counts and across-partition best-match Jaccard summaries at eight Leiden resolutions.
- **Part F, resolution-sweep per-partition results:** size and mean best-match Jaccard index for each reconstructed full-cohort partition at each resolution.
- **Data dictionary:** row unit and purpose of each part.

## General notes

1. Parts A–B and Parts C–F evaluate different computational objects. Parts A–B hold fixed the neighborhood graph stored in `adata_microglia_subtyped.h5ad`. Parts C–F rebuild 15-neighbor graphs from the first 20 principal components inherited from the whole-cell analysis.
2. The reference population contains 7,461 cells and seven numeric partitions at Leiden resolution 0.4. The names Inflammatory, Transitional-A/B, Homeostatic-A/B/C, and Rare are legacy labels rather than validated biological states. Partition 6 (Rare) contains 94 cells, is doublet-enriched, and is not interpreted as a biological state.
3. An adjusted Rand index of 1 denotes identical partitions. Parts A–B quantify algorithmic agreement among partitions obtained from the same cells and stored graph; they contain no biological-replicate test or *P* value.
4. Each conditional subsampling run draws 5,969 of the 7,461 cells (80%) without replacement and does not resample animals. Cell selection, features, normalization, highly variable genes, and principal components are held fixed, so the procedure does not represent an end-to-end pipeline bootstrap.
5. In the variable-seed analysis, the neighborhood-graph and Leiden seeds vary together from 0–99 across 100 runs. In the fixed-seed analysis, both seeds equal 0 across 50 runs. For each reference partition, the best-match Jaccard index is the largest intersection-over-union with any reconstructed partition among the sampled cells. Standard deviations and recovery proportions describe the empirical distribution across runs; they are not confidence intervals or independent replication.
6. A mean Jaccard index of at least 0.75 is shown as a descriptive reference and is not a biological validation threshold.
7. For the resolution sweep, the full-cohort reference graph and Leiden solution use seed 0. Thirty 80% subsample graphs are reconstructed at each resolution with graph and Leiden seeds varying together from 0–29. Reconstructed partition identifiers are local to each resolution and do not map to the reference partition numbers.
8. The partitions meeting the descriptive 0.75 reference at resolutions 0.5 and 0.6 included a 993-cell and a 990-cell partition, respectively. Higher Jaccard values were therefore not restricted to the small doublet-enriched partitions. These conditional analyses nevertheless do not establish stable biological states or a transcriptional continuum.
9. All analyses use cells as computational subsamples. They add no animals and do not constitute biological replication or orthogonal validation.
