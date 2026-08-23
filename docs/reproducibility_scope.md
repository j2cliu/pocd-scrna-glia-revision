# Reproducibility scope and claim boundary

## Included

- Final Figures 1–4 and Figures S1–S3 as PDF, 300-dpi TIFF, PNG preview, and
  serialized R plot object.
- Main Tables 1–2 and machine-readable Tables S1–S6.
- Frozen panel-ready inputs, selected pseudobulk exports, hashes, software
  sessions, execution metadata, and final plotting code.
- Primary animal-level pseudobulk, systematic one-animal-omission,
  background-adjusted scoring, partition-stability, processed-count
  sensitivity, and external-context code used in the revision.
- Canonical revision Methods, Supplementary Material, legends, and table notes
  as frozen documentation snapshots.

## Deliberately excluded

- Raw GEO matrices and large processed AnnData objects.
- The isolated non-IFN discovery protocol/execution branch.
- GSE222430/GSE303920 and Suo 2022 feasibility branches not integrated into
  the manuscript.
- The mean-UMI depth-thinning exploration, which is not reported in the final
  manuscript.
- Retired immediate-early-gene display/table payload.
- Obsolete IRF7-regulatory, 173-gene, TNF–PTPRS, graphical-abstract, and other
  mechanistic analyses inconsistent with the revised claim boundary.
- Credentials, passphrases, encrypted crosswalks, and machine-local paths.

The machine-readable exclusion list is
`provenance/EXCLUSION_POLICY.csv`.

## Interpretation

This repository makes the analysis transparent; it does not increase the
number of biological replicates. Repeated scoring methods, cell-inclusion
rules, processed-count matrices from the same libraries, and animal-omission
refits are sensitivity or influence analyses. They are not independent
biological validation.

