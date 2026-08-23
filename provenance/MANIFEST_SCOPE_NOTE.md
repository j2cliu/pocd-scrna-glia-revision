# Scope of frozen execution manifests

Several original execution manifests record both CSV and TSV versions emitted
by upstream scripts, and the main-table manifest also records a Word export.
To avoid duplicate large files and submission-format artifacts, this public
package retains the CSV or CSV.GZ representation used by the revised
supplementary inventory and the machine-readable main-table CSV files.

Rows for omitted duplicate TSV files or the Word export remain in the frozen
execution manifests solely to preserve the original execution record. Their
absence from the repository is intentional and is not a failed copy. The
authoritative list of files actually curated into this package is
`INCLUSION_MANIFEST.csv`; the authoritative current file hashes are in
`REPOSITORY_FILE_MANIFEST.csv`.

The repository-facing `RESULTS_EVIDENCE_MAP.csv` uses only paths that exist in
the clean package. The internal evidence ledger is retained for historical
claim provenance, including its explicitly retired R4.P1 row and pre-cleanup
Table S4 numbering; it is not the repository navigation index.

