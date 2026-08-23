#!/usr/bin/env python3
"""Build a Scrublet-negative GSE267933 animal-level pseudobulk matrix.

This bounded audit reads the frozen project inputs but writes only beside this
script. It removes cells called ``predicted_doublet`` in the annotated
microglial object, recovers integer UMI counts for the retained barcodes from
``adata_raw.h5ad``, aggregates by animal/library, and applies a full-cohort
gene-total filter of at least 10 UMI. The filtered gene universe is frozen for
all downstream leave-one-animal-out fits.
"""

from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT = Path(os.environ["POCD_SCRNA_PROJECT_ROOT"]).expanduser().resolve()
PROCESSED = PROJECT / "data" / "processed"
SCRUBLET_PATH = PROCESSED / "adata_microglia_scrublet_annotated.h5ad"
RAW_PATH = PROCESSED / "adata_raw.h5ad"
REFERENCE_ANNOTATION_PATH = (
    PROJECT
    / "data"
    / "results"
    / "GSE267933_pseudobulk"
    / "deseq2_all_microglia_annotated.csv"
)
OUT_DIR = Path(__file__).resolve().parent

COUNTS_OUT = OUT_DIR / "scrublet_negative_pseudobulk_counts.csv.gz"
ANNOTATION_OUT = OUT_DIR / "scrublet_negative_gene_annotation.csv"
META_OUT = OUT_DIR / "scrublet_negative_sample_meta.csv"
AUDIT_OUT = OUT_DIR / "scrublet_negative_count_audit.csv"
MANIFEST_OUT = OUT_DIR / "scrublet_negative_input_manifest.json"

EXPECTED_MICROGLIA = 7461
EXPECTED_DOUBLETS = 90
MIN_TOTAL_COUNT = 10


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_csr(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def integer_residual(matrix: sparse.csr_matrix) -> float:
    if matrix.nnz == 0:
        return 0.0
    return float(np.max(np.abs(matrix.data - np.rint(matrix.data))))


def main() -> None:
    scrub = ad.read_h5ad(SCRUBLET_PATH)
    raw = ad.read_h5ad(RAW_PATH)

    if scrub.n_obs != EXPECTED_MICROGLIA:
        raise RuntimeError(
            f"Expected {EXPECTED_MICROGLIA} microglia, found {scrub.n_obs}"
        )
    required = {"sample", "group", "predicted_doublet"}
    missing = required.difference(scrub.obs.columns)
    if missing:
        raise RuntimeError(f"Missing Scrublet metadata columns: {sorted(missing)}")
    if scrub.obs["predicted_doublet"].isna().any():
        raise RuntimeError("predicted_doublet contains missing values")
    if not pd.api.types.is_bool_dtype(scrub.obs["predicted_doublet"]):
        raise RuntimeError(
            "predicted_doublet is not Boolean; refusing an implicit conversion"
        )

    doublet = scrub.obs["predicted_doublet"].to_numpy(dtype=bool)
    n_doublet = int(doublet.sum())
    if n_doublet != EXPECTED_DOUBLETS:
        raise RuntimeError(
            f"Expected {EXPECTED_DOUBLETS} predicted doublets, found {n_doublet}"
        )
    if not scrub.obs_names.is_unique or not raw.obs_names.is_unique:
        raise RuntimeError("Cell barcodes are not unique")
    missing_barcodes = scrub.obs_names.difference(raw.obs_names)
    if len(missing_barcodes):
        raise RuntimeError(
            f"{len(missing_barcodes)} annotated microglial barcodes are absent "
            "from adata_raw.h5ad"
        )

    # Align raw counts exactly to the annotated microglial barcode order.
    mg = raw[scrub.obs_names].copy()
    x = as_csr(mg.X)
    residual = integer_residual(x)
    if residual > 1e-6:
        raise RuntimeError(
            "adata_raw.h5ad does not contain integer-valued UMI counts "
            f"(maximum residual {residual})"
        )

    sample = scrub.obs["sample"].astype(str)
    group = scrub.obs["group"].astype(str)
    raw_sample = mg.obs["sample"].astype(str)
    raw_group = mg.obs["group"].astype(str)
    if not np.array_equal(sample.to_numpy(), raw_sample.to_numpy()):
        raise RuntimeError("Sample labels disagree between Scrublet and raw objects")
    if not np.array_equal(group.to_numpy(), raw_group.to_numpy()):
        raise RuntimeError("Group labels disagree between Scrublet and raw objects")

    sample_ids = sorted(sample.unique())
    expected_samples = ["C1", "C2", "C3", "S1", "S2", "S3"]
    if sample_ids != expected_samples:
        raise RuntimeError(
            f"Expected samples {expected_samples}, found {sample_ids}"
        )
    group_map = (
        pd.DataFrame({"sample": sample, "group": group})
        .drop_duplicates()
        .set_index("sample")["group"]
        .to_dict()
    )
    if any(
        sid.startswith("S") != (group_map[sid] == "Surgery")
        for sid in sample_ids
    ):
        raise RuntimeError(f"Sample prefix/group mismatch: {group_map}")

    cell_totals = np.asarray(x.sum(axis=1)).ravel()
    pseudobulk_columns: dict[str, np.ndarray] = {}
    audit_rows: list[dict] = []
    for sid in sample_ids:
        sample_mask = sample.to_numpy() == sid
        retained_mask = sample_mask & ~doublet
        removed_mask = sample_mask & doublet
        retained_counts = np.asarray(x[retained_mask].sum(axis=0)).ravel()
        pseudobulk_columns[sid] = np.rint(retained_counts).astype(np.int64)
        audit_rows.append(
            {
                "sample": sid,
                "group": group_map[sid],
                "cells_before": int(sample_mask.sum()),
                "predicted_doublets_removed": int(removed_mask.sum()),
                "cells_retained": int(retained_mask.sum()),
                "percent_cells_removed": float(
                    100 * removed_mask.sum() / sample_mask.sum()
                ),
                "umi_before": int(round(cell_totals[sample_mask].sum())),
                "umi_removed": int(round(cell_totals[removed_mask].sum())),
                "umi_retained": int(round(cell_totals[retained_mask].sum())),
                "percent_umi_removed": float(
                    100
                    * cell_totals[removed_mask].sum()
                    / cell_totals[sample_mask].sum()
                ),
                "median_umi_retained_cell": float(
                    np.median(cell_totals[retained_mask])
                ),
            }
        )

    pseudobulk = pd.DataFrame(pseudobulk_columns, index=mg.var_names)
    full_gene_totals = pseudobulk.sum(axis=1)
    keep_gene = full_gene_totals >= MIN_TOTAL_COUNT
    filtered = pseudobulk.loc[keep_gene].copy()
    filtered.index.name = "gene"
    filtered.to_csv(COUNTS_OUT, compression="gzip")

    # Reuse the exact Ensembl-to-symbol mapping used by the all-cell primary
    # analysis. Using the feature-file symbols here would mix doublet exclusion
    # with a gene-nomenclature update and would not be a one-variable sensitivity.
    reference_annotation = pd.read_csv(
        REFERENCE_ANNOTATION_PATH,
        usecols=["gene", "symbol"],
    )
    if reference_annotation["gene"].duplicated().any():
        raise RuntimeError("Reference Ensembl-to-symbol mapping is not unique")
    annotation = (
        pd.DataFrame({"gene": mg.var_names[keep_gene]})
        .merge(
            reference_annotation,
            on="gene",
            how="left",
            validate="one_to_one",
            indicator=True,
        )
    )
    if not (annotation["_merge"] == "both").all():
        raise RuntimeError(
            "At least one Scrublet-negative feature is absent from the "
            "all-cell reference annotation"
        )
    annotation = annotation.drop(columns="_merge")
    annotation.to_csv(ANNOTATION_OUT, index=False)

    meta = pd.DataFrame(
        {
            "sample": sample_ids,
            "group": [group_map[sid] for sid in sample_ids],
            "n_cells": [
                int(((sample.to_numpy() == sid) & ~doublet).sum())
                for sid in sample_ids
            ],
        }
    )
    meta.to_csv(META_OUT, index=False)
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(AUDIT_OUT, index=False)

    manifest = {
        "analysis": "GSE267933 Scrublet-negative pseudobulk extraction",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project_mutation": False,
        "inputs": {
            "scrublet_annotation": {
                "path": str(SCRUBLET_PATH),
                "size_bytes": SCRUBLET_PATH.stat().st_size,
                "sha256": sha256_file(SCRUBLET_PATH),
                "shape": list(scrub.shape),
                "predicted_doublet_column": "predicted_doublet",
            },
            "integer_counts": {
                "path": str(RAW_PATH),
                "size_bytes": RAW_PATH.stat().st_size,
                "sha256": sha256_file(RAW_PATH),
                "shape": list(raw.shape),
                "integer_max_residual_in_microglia": residual,
            },
            "reference_gene_annotation": {
                "path": str(REFERENCE_ANNOTATION_PATH),
                "size_bytes": REFERENCE_ANNOTATION_PATH.stat().st_size,
                "sha256": sha256_file(REFERENCE_ANNOTATION_PATH),
                "purpose": (
                    "reuse the all-cell Ensembl-to-symbol map so doublet "
                    "exclusion is the only changed analytic input"
                ),
            },
        },
        "design": {
            "biological_unit": "animal/library",
            "contrast": "Surgery minus Control",
            "samples": group_map,
            "cells_before": int(scrub.n_obs),
            "predicted_doublets_removed": n_doublet,
            "cells_retained": int((~doublet).sum()),
            "full_cohort_gene_total_filter": f">={MIN_TOTAL_COUNT} UMI",
            "features_before_filter": int(pseudobulk.shape[0]),
            "features_after_filter": int(filtered.shape[0]),
            "allcell_filtered_features": int(len(reference_annotation)),
            "features_lost_after_doublet_exclusion": int(
                len(reference_annotation) - filtered.shape[0]
            ),
            "loo_gene_universe": "frozen full-cohort filtered feature set",
        },
        "software": {
            "python": platform.python_version(),
            "anndata": version("anndata"),
            "numpy": version("numpy"),
            "pandas": version("pandas"),
            "scipy": version("scipy"),
        },
        "outputs": {
            path.name: {
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in (COUNTS_OUT, ANNOTATION_OUT, META_OUT, AUDIT_OUT)
        },
    }
    MANIFEST_OUT.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    print("Scrublet-negative pseudobulk extraction complete")
    print(audit.to_string(index=False))
    print(
        f"Cells: {scrub.n_obs} - {n_doublet} = {(~doublet).sum()}; "
        f"features: {pseudobulk.shape[0]} -> {filtered.shape[0]} "
        f"(total UMI >= {MIN_TOTAL_COUNT})"
    )


if __name__ == "__main__":
    main()
