# Software environments

Frozen R session information is stored in each analysis `manifests/` directory.
The final figure executions used R 4.5.2 and recorded package versions for all
loaded plotting and analysis libraries. Python package versions used by the
matrix and object audits are recorded in the corresponding execution manifests
or scripts.

Key R/Bioconductor packages include DESeq2, clusterProfiler, fgsea, msigdbr,
org.Mm.eg.db, AnnotationDbi, ggplot2, patchwork, readr, dplyr, tidyr, ragg,
systemfonts, digest, and jsonlite. Key Python packages include anndata, scanpy,
numpy, pandas, scipy, and Scrublet.

No lockfile is asserted to be an exact cross-platform environment. The frozen
session files and source/output hashes are the authoritative execution record.

The read-only scripts in `analysis/source_audit/scripts` use the environment
variable `POCD_SCRNA_PROJECT_ROOT` for the local project containing the public
raw data and frozen processed objects:

```bash
export POCD_SCRNA_PROJECT_ROOT=/path/to/pocd_scrna
```
