# Data access and expected inputs

## Public accessions

| Accession | Repository role | Raw data included here? |
|---|---|---|
| [GSE267933](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267933) | Primary six-animal single-cell cohort | No |
| [GSE289098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289098) | Alternative processed-count matrix from the same six libraries/cells | No |
| [GSE283401](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283401) | External cross-experiment bulk-microglial context | No |

Raw matrices are not duplicated because they are publicly deposited and, for
GSE267933/GSE289098, substantially larger than the derived analysis payload.
The repository includes the hashes, dimensions, sample mappings, pseudobulk
counts, panel-ready data, and final machine-readable results required to audit
the revised claims.

## Expected project-root layout for upstream scripts

Scripts accepting `--project-root` expect the following structure after data
are downloaded and prepared:

```text
PROJECT_ROOT/
└── data/
    ├── raw/
    │   ├── GSE267933_family.soft.gz
    │   ├── GSE289098/
    │   │   ├── matrix.mtx.gz
    │   │   ├── barcodes.tsv.gz
    │   │   ├── features.tsv.gz
    │   │   └── protocols.txt.gz
    │   └── GSE283401/
    │       ├── GSE283401_family.soft
    │       ├── GSE283401_hp_6h_counts.csv.gz
    │       └── GSE283401_hp_48h_counts.csv.gz
    └── processed/
        ├── adata_raw.h5ad
        ├── adata_microglia_pseudotime.h5ad
        └── adata_microglia_scrublet_annotated.h5ad
```

Exact expected SHA-256 values for GSE289098/GSE267933 inputs are stored in
`analysis/figureS3/data/preflight/gse289098_input_manifest.csv`. GSE283401
source hashes are stored in
`analysis/figure4/manifests/fig04_source_manifest.csv`.

The processed AnnData objects preserve the exact submitted reference-cell
population and are not redistributed in this repository. Figure and table
reproduction from the included panel-ready and pseudobulk-derived files does
not require those objects.

`adata_microglia_pseudotime.h5ad` is a legacy upstream filename. The retained
Figure S1 preparation script reads only its stored cell identifiers, sample and
group fields, numeric reference-partition labels, and traceability labels. It
does not read a pseudotime value or perform diffusion-pseudotime or directional
trajectory analysis.

## Data-role boundary

GSE289098 must not be counted as a second cohort: its barcodes and ordered
features were matched to the GSE267933 cells and libraries. GSE283401 differs
in assay, anesthetic, exposure duration, laboratory, sampling design, and
sample time. It is therefore used only for external context.
