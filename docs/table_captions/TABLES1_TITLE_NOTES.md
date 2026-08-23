# Table S1 title and notes

## Title

**Table S1. Reference-partition summary and one-versus-rest marker statistics for GSE267933 microglia.**

## Part A — Reference-partition summary

Part A reports the numeric reference partition, legacy label, total cell count, numbers contributed by oxygen-control and combined-exposure libraries, Scrublet-predicted doublet count and rate, group-specific doublet rates, and interpretive scope. Partition 6 (legacy label, Rare) contained 82 predicted doublets among 94 cells and was excluded from biological interpretation.

## Part B — One-versus-rest marker statistics

Part B reports all positive one-versus-rest markers with a Benjamini–Hochberg adjusted *P* value below 0.05 and a positive Scanpy log2 fold-change estimate after normalization to 10,000 counts per cell and log1p transformation. Reported fields include rank, gene symbol, Wilcoxon score, Scanpy log2 fold-change estimate, nominal and adjusted *P* values, and detection percentages within the partition and among all remaining cells. `top_50_marker` identifies the first 50 retained markers for inspection and is not an additional selection criterion.

## General notes

- Numeric partitions 0–6 define the reference solution. Legacy names identify the original annotation and are not treated as validated biological states.
- Marker statistics use cells as observations and serve as descriptive annotation measures rather than animal-level differential-expression tests.
- These results do not establish partition stability, differential abundance, or biological-state validity.
- BAM/CAM denotes border- or central-nervous-system-associated macrophage.
