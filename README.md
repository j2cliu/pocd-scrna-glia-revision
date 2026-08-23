# GLIA revision reproducibility package

This repository contains the frozen analysis code, derived data, final figures,
machine-readable tables, execution metadata, and claim-provenance materials for
the revised manuscript:

> **Animal-level heterogeneity in interferon-responsive transcription after
> sevoflurane anesthesia and laparotomy in aged hippocampal microglia**

The package is specific to the GLIA revision. It does not contain analyses from
retired mechanistic narratives, unrelated follow-up projects, or the isolated
non-IFN discovery branch.

## Scientific scope

The primary dataset is GSE267933: hippocampal single-cell RNA sequencing from
three 18-month-old male mice receiving oxygen control and three receiving
sevoflurane anesthesia plus laparotomy. The mouse and corresponding library,
not the individual cell, is the biological unit for between-group analyses.

- **GSE267933** is the six-animal primary cohort.
- **GSE289098** contains an alternative processed-count matrix derived from the
  same six libraries and cells. It is a same-cohort processing sensitivity
  analysis, not replication or orthogonal validation.
- **GSE283401** supplies external cross-experiment context from bulk RNA
  sequencing of isolated hippocampal microglia. Its 6-h and 48-h experiments
  use different animals and separately deposited matrices. It is not a
  longitudinal time course, replication, or validation of GSE267933.

The revised work is exploratory. Direct expression, background-adjusted
`score_genes` values, and Hallmark enrichment statistics are distinct
estimands. Systematic one-animal omissions are influence diagnostics rather
than independent tests, and no animal is removed from the full-cohort result.

## Repository map

| Location | Contents |
|---|---|
| `analysis/figure1` | Study design, submitted-partition audit, and animal-level composition |
| `analysis/tables1_2` | Main Tables 1–2 |
| `analysis/figure2` | Primary pseudobulk results, selected transcripts, and Table S3 |
| `analysis/figure3` | Estimator/animal-influence displays and clean Table S4 |
| `analysis/figure4` | GSE283401 external context and Table S5 |
| `analysis/figureS1` | Partition annotation/doublet audit and Table S1 |
| `analysis/figureS2` | Partition-stability diagnostics and Table S2 |
| `analysis/figureS3` | GSE289098 same-cohort count sensitivity and Table S6 |
| `analysis/source_audit` | Read-only audit scripts and frozen derived inputs supporting current analyses |
| `docs` | Methods, supplementary material, legends, table notes, and scope documentation |
| `provenance` | Inclusion/exclusion records, evidence ledger, checksums, and prepublication audit |
| `validation` | Repository and source-workspace validation scripts |

A display-to-code-to-data index is available in
[`docs/analysis_map.md`](docs/analysis_map.md).

## Reproducibility levels

### 1. Self-contained display reproduction

The repository includes the frozen panel-ready data and final plotting code for
Figures 1–4 and Figures S1–S3. These plots can be regenerated without
downloading the raw GEO matrices. Run the commands from the repository root:

```bash
Rscript analysis/figure1/scripts/79_plot_glia_R1_figure1.R
Rscript analysis/figure2/scripts/84_plot_glia_R1_figure2_rebalanced.R
Rscript analysis/figure3/scripts/85_plot_glia_R1_figure3.R
Rscript analysis/figure4/scripts/87_plot_glia_R1_figure4.R
Rscript analysis/figureS1/scripts/92_plot_glia_R1_figureS1.R
Rscript analysis/figureS2/scripts/92_plot_glia_R1_figureS2.R
Rscript analysis/figureS3/scripts/90_plot_figureS3.R
```

The final vector PDFs were produced with the macOS Quartz device and Arial;
TIFF and PNG files were rendered with `ragg`. On another operating system,
numeric results should be unchanged, but fonts, PDF metadata, and raster hashes
may differ. Session information for each frozen execution is retained beside
the corresponding analysis.

### 2. Analysis recomputation from public data

Data-preparation and audit scripts are included where they support the revised
analyses. Raw GEO matrices and large AnnData objects are not duplicated in this
repository. Their expected locations and SHA-256 identities are described in
[`docs/data_access.md`](docs/data_access.md) and the analysis manifests.

The exact submitted reference cell population is preserved through frozen
object-derived exports. As stated in the Supplementary Methods, the available
executable record does not contain every intermediate required to regenerate
the submitted whole-cell object and reference partition exactly from raw
matrices. The package therefore makes a narrower, testable claim: the final
statistics, tables, and displays are traceable to frozen derived inputs, while
raw-data scripts and source hashes document the upstream dependencies.

## Public datasets

- [GSE267933](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267933)
- [GSE289098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289098)
- [GSE283401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283401)

## Integrity and audit

Every curated source copy is recorded in
[`provenance/INCLUSION_MANIFEST.csv`](provenance/INCLUSION_MANIFEST.csv), with
source and repository SHA-256 values. Explicit exclusions are listed in
[`provenance/EXCLUSION_POLICY.csv`](provenance/EXCLUSION_POLICY.csv).
The 16 active Results units are linked to their final repository data, code,
display, numeric anchors, and claim ceilings in
[`provenance/RESULTS_EVIDENCE_MAP.csv`](provenance/RESULTS_EVIDENCE_MAP.csv).

Run the bounded public-package audit with:

```bash
Rscript validation/00_validate_public_repository.R
```

The audit verifies required displays/tables, parses all R scripts, compiles all
Python scripts, scans text for credentials and private absolute paths, checks
file sizes, and writes the repository file manifest and audit report.

## Citation and license status

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). A public
license has deliberately not been assigned during local staging; see
[`LICENSE_PENDING.md`](LICENSE_PENDING.md). No remote repository has been
created and no file has been published from this staging package.
