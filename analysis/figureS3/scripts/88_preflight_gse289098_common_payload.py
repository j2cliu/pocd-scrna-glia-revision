#!/usr/bin/env python3
"""Bounded GSE289098/GSE267933 same-cohort processed-count preflight.

This script does not cluster cells or perform cell-level inference. It:

1. verifies the six-library provenance and barcode-suffix mapping;
2. maps the GSE289098 integrated matrix one-to-one to GSE267933 cells;
3. freezes the same 7,371 Scrublet-negative microglial cells;
4. demonstrates and quantifies processed-count payload differences; and
5. exports matched animal-level pseudobulk inputs for the R estimand analysis.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.io import mmread


SAMPLES = ["C1", "C2", "C3", "S1", "S2", "S3"]
GROUPS = {
    "C1": "Control",
    "C2": "Control",
    "C3": "Control",
    "S1": "Surgery",
    "S2": "Surgery",
    "S3": "Surgery",
}
LIBRARIES = {
    "C1": {
        "suffix": "6",
        "gsm": "GSM8281758",
        "title": "Hippocampus_Control_rep1",
        "biosample": "SAMN41463843",
        "sra_experiment": "SRX24614046",
    },
    "C2": {
        "suffix": "5",
        "gsm": "GSM8281759",
        "title": "Hippocampus_Control_rep2",
        "biosample": "SAMN41463842",
        "sra_experiment": "SRX24614047",
    },
    "C3": {
        "suffix": "2",
        "gsm": "GSM8281760",
        "title": "Hippocampus_Control_rep3",
        "biosample": "SAMN41463841",
        "sra_experiment": "SRX24614048",
    },
    "S1": {
        "suffix": "1",
        "gsm": "GSM8281761",
        "title": "Hippocampus_Surgery_rep1",
        "biosample": "SAMN41463840",
        "sra_experiment": "SRX24614049",
    },
    "S2": {
        "suffix": "3",
        "gsm": "GSM8281762",
        "title": "Hippocampus_Surgery_rep2",
        "biosample": "SAMN41463839",
        "sra_experiment": "SRX24614050",
    },
    "S3": {
        "suffix": "4",
        "gsm": "GSM8281763",
        "title": "Hippocampus_Surgery_rep3",
        "biosample": "SAMN41463838",
        "sra_experiment": "SRX24614051",
    },
}
SELECTED_GENES = ["Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_true(condition: bool, message: str) -> None:
    if not bool(condition):
        raise RuntimeError(message)


def as_csr(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def integer_residual(matrix: sparse.csr_matrix) -> float:
    if matrix.nnz == 0:
        return 0.0
    return float(np.max(np.abs(matrix.data - np.rint(matrix.data))))


def read_gzip_lines(path: Path) -> list[str]:
    # GEO protocol text may contain legacy Windows-1252 punctuation.  The
    # identifiers used for provenance checks are ASCII, so replacement is
    # preferable to aborting on a non-UTF-8 typographic character.
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\n") for line in handle]


def primary_to_integrated_barcode(primary_barcode: str) -> str:
    sample, tenx_barcode = primary_barcode.split("_", 1)
    assert_true(sample in LIBRARIES, f"Unexpected primary sample prefix: {sample}")
    base = tenx_barcode.rsplit("-", 1)[0]
    return f"{base}-{LIBRARIES[sample]['suffix']}"


def write_csv(data: pd.DataFrame, path: Path) -> None:
    data.to_csv(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--figure-root", required=True)
    args = parser.parse_args()

    project = Path(args.project_root).resolve()
    figure_root = Path(args.figure_root).resolve()
    out_dir = figure_root / "data" / "preflight"
    manifest_dir = figure_root / "manifests"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    processed = project / "data" / "processed"
    raw_dir = project / "data" / "raw"
    gse289 = raw_dir / "GSE289098"
    primary_path = processed / "adata_raw.h5ad"
    scrublet_path = processed / "adata_microglia_scrublet_annotated.h5ad"
    family_soft = raw_dir / "GSE267933_family.soft.gz"
    matrix_path = gse289 / "matrix.mtx.gz"
    barcode_path = gse289 / "barcodes.tsv.gz"
    feature_path = gse289 / "features.tsv.gz"
    protocol_path = gse289 / "protocols.txt.gz"
    series_stub = gse289 / "GSE289098_series_matrix.txt.gz"
    raw_tar_stub = gse289 / "GSE289098_RAW.tar"
    primary_pseudobulk_reference = (
        figure_root.parents[2]
        / "independent_audit"
        / "scrublet_negative_pseudobulk_counts.csv.gz"
    )

    required_paths = [
        primary_path,
        scrublet_path,
        family_soft,
        matrix_path,
        barcode_path,
        feature_path,
        protocol_path,
        series_stub,
        raw_tar_stub,
        primary_pseudobulk_reference,
    ]
    for path in required_paths:
        assert_true(path.exists(), f"Missing input: {path}")

    protocols = "\n".join(read_gzip_lines(protocol_path))
    with gzip.open(family_soft, "rt", errors="replace") as handle:
        soft = handle.read()
    for sample, item in LIBRARIES.items():
        expected_protocol = f"{item['title']} ---> barcode-{item['suffix']}"
        assert_true(
            expected_protocol in protocols,
            f"Protocol suffix mapping not found: {expected_protocol}",
        )
        for token in [item["gsm"], item["title"], item["biosample"], item["sra_experiment"]]:
            assert_true(token in soft, f"GSE267933 SOFT lacks provenance token: {token}")
        relation = f"!Sample_relation = Reanalyzed by: GSE289098"
        assert_true(relation in soft, "GSE267933 SOFT lacks the GSE289098 relation")

    barcodes = read_gzip_lines(barcode_path)
    features_split = [line.split("\t") for line in read_gzip_lines(feature_path)]
    assert_true(all(len(row) >= 2 for row in features_split), "Malformed feature table")
    feature_ids = [row[0] for row in features_split]
    feature_symbols = [row[1] for row in features_split]
    assert_true(len(barcodes) == 20684, "Expected 20,684 integrated barcodes")
    assert_true(len(set(barcodes)) == len(barcodes), "Integrated barcodes are duplicated")
    assert_true(len(feature_ids) == 27998, "Expected 27,998 features")
    assert_true(len(set(feature_ids)) == len(feature_ids), "Integrated feature IDs are duplicated")

    print("Loading GSE289098 integrated Matrix Market payload...", flush=True)
    with gzip.open(matrix_path, "rb") as handle:
        integrated_features_by_cells = mmread(handle)
    integrated = integrated_features_by_cells.T.tocsr()
    del integrated_features_by_cells
    gc.collect()
    assert_true(integrated.shape == (20684, 27998), "Unexpected integrated matrix shape")
    assert_true(integer_residual(integrated) == 0.0, "Integrated counts are non-integer")
    assert_true(integrated.data.min() >= 0, "Integrated counts contain negative values")

    print("Loading GSE267933 primary count object...", flush=True)
    primary = ad.read_h5ad(primary_path)
    primary_counts = as_csr(primary.X)
    assert_true(primary.shape == (20684, 27998), "Unexpected primary object shape")
    assert_true(integer_residual(primary_counts) == 0.0, "Primary counts are non-integer")
    assert_true(primary_counts.data.min() >= 0, "Primary counts contain negative values")
    assert_true(primary.obs_names.is_unique, "Primary cell identifiers are duplicated")
    assert_true(list(primary.var_names) == feature_ids, "Feature ID order differs")
    assert_true(
        primary.var["gene_symbol"].astype(str).tolist() == feature_symbols,
        "Feature-symbol order differs",
    )

    mapped_barcodes = [primary_to_integrated_barcode(x) for x in primary.obs_names]
    integrated_index = pd.Index(barcodes)
    mapping_index = integrated_index.get_indexer(mapped_barcodes)
    assert_true(np.all(mapping_index >= 0), "At least one primary barcode is absent")
    assert_true(len(np.unique(mapping_index)) == len(mapping_index), "Barcode mapping is not one-to-one")
    assert_true(set(mapped_barcodes) == set(barcodes), "Full barcode universes differ")

    scrublet = ad.read_h5ad(scrublet_path, backed="r")
    required_obs = {"sample", "group", "predicted_doublet"}
    assert_true(required_obs.issubset(scrublet.obs.columns), "Scrublet metadata is incomplete")
    assert_true(scrublet.n_obs == 7461, "Expected 7,461 submitted microglia")
    assert_true(scrublet.obs_names.is_unique, "Scrublet cell identifiers are duplicated")
    assert_true(
        pd.api.types.is_bool_dtype(scrublet.obs["predicted_doublet"]),
        "predicted_doublet is not Boolean",
    )
    retained_names = scrublet.obs_names[~scrublet.obs["predicted_doublet"].to_numpy()]
    assert_true(len(retained_names) == 7371, "Expected 7,371 Scrublet-negative cells")
    retained_primary_index = primary.obs_names.get_indexer(retained_names)
    assert_true(np.all(retained_primary_index >= 0), "Retained cells are absent from primary counts")
    retained_integrated_index = mapping_index[retained_primary_index]
    assert_true(len(np.unique(retained_integrated_index)) == 7371, "Retained mapping is not unique")

    primary_cell_totals = np.asarray(primary_counts.sum(axis=1)).ravel().astype(np.int64)
    integrated_cell_totals_native = np.asarray(integrated.sum(axis=1)).ravel().astype(np.int64)
    integrated_cell_totals = integrated_cell_totals_native[mapping_index]
    primary_cell_nnz = np.diff(primary_counts.indptr).astype(np.int64)
    integrated_cell_nnz = np.diff(integrated.indptr).astype(np.int64)[mapping_index]

    print("Comparing full matched count payloads...", flush=True)
    n_different_entries = 0
    n_primary_greater = 0
    n_integrated_greater = 0
    absolute_difference_sum = 0
    signed_primary_minus_integrated_sum = 0
    maximum_absolute_difference = 0
    cells_with_any_difference = np.zeros(primary.n_obs, dtype=bool)
    chunk_size = 256
    for start in range(0, primary.n_obs, chunk_size):
        end = min(start + chunk_size, primary.n_obs)
        local_primary = primary_counts[start:end]
        local_integrated = integrated[mapping_index[start:end]]
        difference = (local_primary - local_integrated).tocsr()
        difference.eliminate_zeros()
        if difference.nnz:
            data = difference.data.astype(np.int64)
            n_different_entries += int(difference.nnz)
            n_primary_greater += int(np.sum(data > 0))
            n_integrated_greater += int(np.sum(data < 0))
            absolute_difference_sum += int(np.abs(data).sum())
            signed_primary_minus_integrated_sum += int(data.sum())
            maximum_absolute_difference = max(
                maximum_absolute_difference, int(np.abs(data).max())
            )
            cells_with_any_difference[start:end] = difference.getnnz(axis=1) > 0

    # Do not use a single float32 sparse reduction for the primary total: its
    # accumulation can round at this scale.  Per-cell integer totals preserve
    # the exact count sum and must reconcile with the signed entrywise audit.
    primary_total_umi_exact = int(primary_cell_totals.sum(dtype=np.int64))
    integrated_total_umi_exact = int(
        integrated_cell_totals.sum(dtype=np.int64)
    )
    assert_true(
        primary_total_umi_exact - integrated_total_umi_exact
        == signed_primary_minus_integrated_sum,
        "Exact payload totals do not reconcile with entrywise differences",
    )
    assert_true(
        n_integrated_greater == 0
        and absolute_difference_sum == signed_primary_minus_integrated_sum,
        "GSE289098 is not a count-nonincreasing payload of the primary matrix",
    )

    sample_values = primary.obs["sample"].astype(str).to_numpy()
    group_values = primary.obs["group"].astype(str).to_numpy()
    assert_true(set(sample_values) == set(SAMPLES), "Unexpected primary sample labels")
    for sample in SAMPLES:
        assert_true(
            set(group_values[sample_values == sample]) == {GROUPS[sample]},
            f"Primary group mapping failed for {sample}",
        )

    library_rows = []
    barcode_rows = []
    for sample in SAMPLES:
        mask = sample_values == sample
        scrub_mask = scrublet.obs["sample"].astype(str).to_numpy() == sample
        retained_mask = (
            scrub_mask & ~scrublet.obs["predicted_doublet"].to_numpy()
        )
        item = LIBRARIES[sample]
        primary_total = int(primary_cell_totals[mask].sum())
        integrated_total = int(integrated_cell_totals[mask].sum())
        library_rows.append(
            {
                "sample": sample,
                "group": GROUPS[sample],
                "gsm": item["gsm"],
                "integrated_suffix": item["suffix"],
                "n_cells": int(mask.sum()),
                "primary_total_umi": primary_total,
                "integrated_total_umi": integrated_total,
                "umi_removed": primary_total - integrated_total,
                "umi_removed_percent": 100 * (primary_total - integrated_total) / primary_total,
                "integrated_to_primary_umi_ratio": integrated_total / primary_total,
                "primary_nnz": int(primary_cell_nnz[mask].sum()),
                "integrated_nnz": int(integrated_cell_nnz[mask].sum()),
                "nnz_removed": int(
                    primary_cell_nnz[mask].sum() - integrated_cell_nnz[mask].sum()
                ),
                "cells_with_any_payload_difference": int(cells_with_any_difference[mask].sum()),
                "bit_identical_library": bool(not cells_with_any_difference[mask].any()),
                "median_cell_umi_ratio": float(
                    np.median(integrated_cell_totals[mask] / primary_cell_totals[mask])
                ),
            }
        )
        barcode_rows.append(
            {
                "sample": sample,
                "group": GROUPS[sample],
                "gsm": item["gsm"],
                "biosample": item["biosample"],
                "sra_experiment": item["sra_experiment"],
                "integrated_suffix": item["suffix"],
                "primary_all_cells": int(mask.sum()),
                "integrated_mapped_cells": int(mask.sum()),
                "submitted_microglia": int(scrub_mask.sum()),
                "scrublet_negative_microglia": int(retained_mask.sum()),
                "missing_primary_to_integrated": 0,
                "duplicate_integrated_matches": 0,
            }
        )

    primary_target = primary_counts[retained_primary_index]
    integrated_target = integrated[retained_integrated_index]
    retained_sample = scrublet.obs.loc[retained_names, "sample"].astype(str).to_numpy()
    retained_group = scrublet.obs.loc[retained_names, "group"].astype(str).to_numpy()
    assert_true(set(retained_sample) == set(SAMPLES), "Retained sample labels are incomplete")

    target_difference = (primary_target - integrated_target).tocsr()
    target_difference.eliminate_zeros()
    target_difference_data = target_difference.data.astype(np.int64)
    assert_true(
        target_difference.nnz > 0 and np.all(target_difference_data > 0),
        "Common-cell payload is not strictly count-nonincreasing",
    )
    target_primary_total_umi = int(
        np.asarray(primary_target.sum(axis=1))
        .ravel()
        .astype(np.int64)
        .sum(dtype=np.int64)
    )
    target_integrated_total_umi = int(
        np.asarray(integrated_target.sum(axis=1))
        .ravel()
        .astype(np.int64)
        .sum(dtype=np.int64)
    )
    target_cells_with_difference = int(
        np.sum(target_difference.getnnz(axis=1) > 0)
    )
    target_primary_nonzeros_removed = int(
        primary_target.nnz - integrated_target.nnz
    )
    target_shared_nonzeros_reduced = int(
        target_difference.nnz - target_primary_nonzeros_removed
    )
    assert_true(
        target_primary_total_umi - target_integrated_total_umi
        == int(target_difference_data.sum()),
        "Common-cell payload totals do not reconcile with entrywise differences",
    )

    primary_pseudobulk = np.zeros((len(feature_ids), len(SAMPLES)), dtype=np.int64)
    integrated_pseudobulk = np.zeros_like(primary_pseudobulk)
    for column, sample in enumerate(SAMPLES):
        mask = retained_sample == sample
        primary_pseudobulk[:, column] = np.asarray(
            primary_target[mask].sum(axis=0)
        ).ravel().astype(np.int64)
        integrated_pseudobulk[:, column] = np.asarray(
            integrated_target[mask].sum(axis=0)
        ).ravel().astype(np.int64)

    primary_gene_total = primary_pseudobulk.sum(axis=1)
    primary_filter = primary_gene_total >= 10
    assert_true(int(primary_filter.sum()) == 13926, "Primary >=10 universe is not 13,926")

    reference = pd.read_csv(primary_pseudobulk_reference)
    assert_true(reference.columns[0] == "gene", "Unexpected primary pseudobulk reference")
    reference = reference.set_index("gene").loc[np.asarray(feature_ids)[primary_filter], SAMPLES]
    assert_true(
        np.array_equal(reference.to_numpy(dtype=np.int64), primary_pseudobulk[primary_filter]),
        "Recomputed primary pseudobulk does not match the frozen audit",
    )

    symbol_to_indices: dict[str, list[int]] = {}
    for index, symbol in enumerate(feature_symbols):
        if symbol in SELECTED_GENES:
            symbol_to_indices.setdefault(symbol, []).append(index)
    assert_true(set(symbol_to_indices) == set(SELECTED_GENES), "Selected-gene coverage failed")

    selected_rows = []
    panel_rows = []
    retained_count_rows = []
    for sample_index, sample in enumerate(SAMPLES):
        sample_mask = retained_sample == sample
        primary_denominator = int(primary_pseudobulk[:, sample_index].sum())
        integrated_denominator = int(integrated_pseudobulk[:, sample_index].sum())
        retained_count_rows.append(
            {
                "sample": sample,
                "group": GROUPS[sample],
                "n_common_scrublet_negative_cells": int(sample_mask.sum()),
                "primary_total_umi": primary_denominator,
                "integrated_total_umi": integrated_denominator,
                "integrated_to_primary_umi_ratio": integrated_denominator / primary_denominator,
            }
        )
        primary_panel_values = []
        integrated_panel_values = []
        for gene in SELECTED_GENES:
            indices = symbol_to_indices[gene]
            primary_gene_umi = int(primary_pseudobulk[indices, sample_index].sum())
            integrated_gene_umi = int(integrated_pseudobulk[indices, sample_index].sum())
            primary_log2cpm = float(
                np.log2((primary_gene_umi + 0.5) / (primary_denominator + 1) * 1e6)
            )
            integrated_log2cpm = float(
                np.log2((integrated_gene_umi + 0.5) / (integrated_denominator + 1) * 1e6)
            )
            primary_panel_values.append(primary_log2cpm)
            integrated_panel_values.append(integrated_log2cpm)
            selected_rows.append(
                {
                    "sample": sample,
                    "group": GROUPS[sample],
                    "gene": gene,
                    "n_features_summed": len(indices),
                    "primary_gene_umi": primary_gene_umi,
                    "integrated_gene_umi": integrated_gene_umi,
                    "primary_total_umi": primary_denominator,
                    "integrated_total_umi": integrated_denominator,
                    "primary_log2_cpm": primary_log2cpm,
                    "integrated_log2_cpm": integrated_log2cpm,
                    "integrated_minus_primary_log2_cpm": integrated_log2cpm - primary_log2cpm,
                }
            )
        panel_rows.append(
            {
                "sample": sample,
                "group": GROUPS[sample],
                "primary_panel_mean_log2_cpm": float(np.mean(primary_panel_values)),
                "integrated_panel_mean_log2_cpm": float(np.mean(integrated_panel_values)),
                "integrated_minus_primary_panel_value": float(
                    np.mean(integrated_panel_values) - np.mean(primary_panel_values)
                ),
            }
        )

    common_manifest = pd.DataFrame(
        {
            "primary_barcode": retained_names,
            "integrated_barcode": [mapped_barcodes[i] for i in retained_primary_index],
            "sample": retained_sample,
            "group": retained_group,
            "scrublet_predicted_doublet": False,
        }
    )
    common_manifest.to_csv(
        out_dir / "gse289098_common_cell_manifest.csv.gz",
        index=False,
        compression="gzip",
    )

    annotation = pd.DataFrame({"gene": feature_ids, "symbol": feature_symbols})
    write_csv(annotation, out_dir / "gse289098_gene_annotation.csv")
    write_csv(pd.DataFrame(barcode_rows), out_dir / "gse289098_barcode_audit.csv")
    write_csv(pd.DataFrame(library_rows), out_dir / "gse289098_payload_library_audit.csv")
    write_csv(pd.DataFrame(retained_count_rows), out_dir / "gse289098_common_cell_count_audit.csv")
    write_csv(pd.DataFrame(selected_rows), out_dir / "gse289098_selected_gene_payload_values.csv")
    write_csv(pd.DataFrame(panel_rows), out_dir / "gse289098_panel_animal_values.csv")

    filtered_ids = np.asarray(feature_ids)[primary_filter]
    primary_filtered = pd.DataFrame(primary_pseudobulk[primary_filter], columns=SAMPLES)
    primary_filtered.insert(0, "gene", filtered_ids)
    integrated_filtered = pd.DataFrame(integrated_pseudobulk[primary_filter], columns=SAMPLES)
    integrated_filtered.insert(0, "gene", filtered_ids)
    primary_filtered.to_csv(
        out_dir / "gse267933_primary_commoncell_pseudobulk_counts.csv.gz",
        index=False,
        compression="gzip",
    )
    integrated_filtered.to_csv(
        out_dir / "gse289098_integrated_commoncell_pseudobulk_counts.csv.gz",
        index=False,
        compression="gzip",
    )

    mapping_rows = []
    for sample in SAMPLES:
        item = LIBRARIES[sample]
        mapping_rows.append(
            {
                "sample": sample,
                "source_group": GROUPS[sample],
                "gse267933_gsm": item["gsm"],
                "gse267933_title": item["title"],
                "biosample": item["biosample"],
                "sra_experiment": item["sra_experiment"],
                "gse289098_barcode_suffix": item["suffix"],
                "gse267933_relation": "Reanalyzed by GSE289098",
                "biological_replication_added": False,
            }
        )
    write_csv(pd.DataFrame(mapping_rows), out_dir / "gse289098_library_mapping.csv")

    global_audit = pd.DataFrame(
        [
            {
                "n_cells_primary": primary.n_obs,
                "n_cells_integrated": integrated.shape[0],
                "n_features_primary": primary.n_vars,
                "n_features_integrated": integrated.shape[1],
                "n_common_cells_all": len(mapping_index),
                "n_common_scrublet_negative_microglia": len(retained_names),
                "primary_nnz": int(primary_counts.nnz),
                "integrated_nnz": int(integrated.nnz),
                "primary_total_umi": primary_total_umi_exact,
                "integrated_total_umi": integrated_total_umi_exact,
                "umi_removed": primary_total_umi_exact - integrated_total_umi_exact,
                "umi_removed_percent": 100
                * (primary_total_umi_exact - integrated_total_umi_exact)
                / primary_total_umi_exact,
                "integrated_to_primary_umi_ratio": integrated_total_umi_exact
                / primary_total_umi_exact,
                "nnz_removed": int(primary_counts.nnz - integrated.nnz),
                "nnz_removed_percent": 100
                * (primary_counts.nnz - integrated.nnz)
                / primary_counts.nnz,
                "n_different_matrix_entries": n_different_entries,
                "n_entries_primary_greater": n_primary_greater,
                "n_entries_integrated_greater": n_integrated_greater,
                "primary_nonzeros_removed": int(primary_counts.nnz - integrated.nnz),
                "shared_nonzeros_reduced": int(
                    n_different_entries - (primary_counts.nnz - integrated.nnz)
                ),
                "shared_nonzeros_unchanged": int(
                    integrated.nnz
                    - (n_different_entries - (primary_counts.nnz - integrated.nnz))
                ),
                "absolute_difference_sum": absolute_difference_sum,
                "signed_primary_minus_integrated_sum": signed_primary_minus_integrated_sum,
                "maximum_absolute_entry_difference": maximum_absolute_difference,
                "cells_with_any_payload_difference": int(cells_with_any_difference.sum()),
                "target_primary_nnz": int(primary_target.nnz),
                "target_integrated_nnz": int(integrated_target.nnz),
                "target_primary_total_umi": target_primary_total_umi,
                "target_integrated_total_umi": target_integrated_total_umi,
                "target_umi_removed": target_primary_total_umi
                - target_integrated_total_umi,
                "target_umi_removed_percent": 100
                * (target_primary_total_umi - target_integrated_total_umi)
                / target_primary_total_umi,
                "target_different_matrix_entries": int(target_difference.nnz),
                "target_primary_nonzeros_removed": target_primary_nonzeros_removed,
                "target_shared_nonzeros_reduced": target_shared_nonzeros_reduced,
                "target_cells_with_any_payload_difference": target_cells_with_difference,
                "target_cells_bit_identical": int(
                    len(retained_names) - target_cells_with_difference
                ),
                "primary_filtered_gene_universe": int(primary_filter.sum()),
                "integrated_genes_ge10_on_common_cells": int(
                    (integrated_pseudobulk.sum(axis=1) >= 10).sum()
                ),
                "malformed_local_series_matrix_stub": series_stub.stat().st_size <= 1000,
                "malformed_local_raw_tar_stub": raw_tar_stub.stat().st_size <= 1000,
            }
        ]
    )
    write_csv(global_audit, out_dir / "gse289098_payload_global_audit.csv")

    input_rows = []
    source_types = {
        primary_path: "GSE267933 QC-retained integer count object",
        scrublet_path: "Submitted-microglia metadata and Scrublet calls",
        family_soft: "GSE267933 GEO SOFT provenance",
        matrix_path: "GSE289098 integrated processed integer count matrix",
        barcode_path: "GSE289098 integrated barcodes",
        feature_path: "GSE289098 feature identifiers and symbols",
        protocol_path: "GSE289098 processing and barcode-suffix protocol",
        series_stub: "Malformed local download stub; not analytic input",
        raw_tar_stub: "Malformed local download stub; not analytic input",
        primary_pseudobulk_reference: "Frozen primary pseudobulk cross-check",
    }
    for path in required_paths:
        input_rows.append(
            {
                "source": source_types[path],
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "analytic_role": (
                    "excluded malformed stub"
                    if path in {series_stub, raw_tar_stub}
                    else "required or cross-check"
                ),
            }
        )
    write_csv(pd.DataFrame(input_rows), out_dir / "gse289098_input_manifest.csv")

    payload_diff_pass = n_different_entries > 0 and signed_primary_minus_integrated_sum != 0
    gates = pd.DataFrame(
        [
            {
                "gate": "Input",
                "status": "PARTIAL",
                "evidence": (
                    "Required matrix/barcode/feature/protocol and GSE267933 provenance files are valid; "
                    "two 992-byte nonanalytic local download stubs are malformed"
                ),
            },
            {
                "gate": "Identity",
                "status": "PASS",
                "evidence": "All 20,684 cells map one-to-one across six named source libraries",
            },
            {
                "gate": "Common cell",
                "status": "PASS",
                "evidence": "All 7,371 frozen Scrublet-negative microglia map exactly once",
            },
            {
                "gate": "Coverage",
                "status": "PASS",
                "evidence": "All seven selected genes use identical feature IDs and are measurable in both payloads",
            },
            {
                "gate": "Replication",
                "status": "PASS",
                "evidence": "Same six animal/library units retained at n=3 per group; no new replication claimed",
            },
            {
                "gate": "Payload difference",
                "status": "PASS" if payload_diff_pass else "FAIL",
                "evidence": f"{n_different_entries:,} matched matrix entries differ",
            },
            {
                "gate": "Estimand concordance",
                "status": "PENDING",
                "evidence": "Requires the frozen R animal-level effect and leave-one-animal-out analysis",
            },
        ]
    )
    write_csv(gates, out_dir / "gse289098_gate_table.csv")

    output_paths = sorted(out_dir.glob("*"))
    output_manifest = pd.DataFrame(
        [
            {
                "relative_path": str(path.relative_to(figure_root)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in output_paths
            if path.is_file()
        ]
    )
    write_csv(output_manifest, manifest_dir / "figS03_preflight_output_manifest.csv")

    execution = {
        "analysis": "GSE289098 same-cohort processed-count bounded preflight",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "claim_role": "same-six-library processed-count sensitivity only",
        "biological_unit": "animal/library",
        "primary_analysis_cells": 7371,
        "plotting": "none; plotting is reserved for the R Figure S3 script",
        "software": {
            "python": platform.python_version(),
            "anndata": version("anndata"),
            "numpy": version("numpy"),
            "pandas": version("pandas"),
            "scipy": scipy.__version__,
        },
        "gates": gates.to_dict(orient="records"),
    }
    with (manifest_dir / "figS03_preflight_execution.json").open("w") as handle:
        json.dump(execution, handle, indent=2)

    with (manifest_dir / "figS03_preflight_sessionInfo.txt").open("w") as handle:
        handle.write(f"Python: {platform.python_version()}\n")
        handle.write(f"Platform: {platform.platform()}\n")
        for package in ["anndata", "numpy", "pandas", "scipy"]:
            handle.write(f"{package}: {version(package)}\n")

    scrublet.file.close()
    print("Preflight complete.")
    print(global_audit.to_string(index=False))
    print(gates.to_string(index=False))


if __name__ == "__main__":
    main()
