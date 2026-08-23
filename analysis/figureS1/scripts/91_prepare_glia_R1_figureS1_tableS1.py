#!/usr/bin/env python3
"""Prepare the retrospective submitted-partition audit for Figure S1/Table S1.

This is a traceability analysis of the seven-way partition used in the submitted
manuscript.  It does not validate those partitions as reproducible biological
states.  One-versus-rest Wilcoxon statistics are descriptive cell-level marker
statistics and are not biological-replicate tests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


PARTITION_MAP = {
    "0": "Inflammatory",
    "1": "Transitional-A",
    "2": "Homeostatic-A",
    "3": "Homeostatic-B",
    "4": "Homeostatic-C",
    "5": "Transitional-B",
    "6": "Rare",
}

SELECTED_GENES = [
    "P2ry12",
    "Tmem119",
    "Cx3cr1",
    "Tnf",
    "Il1b",
    "Irf7",
    "Ifitm3",
    "Nr1d1",
    "Dbp",
]

DETECTION_GENE_GROUPS = {
    "Microglia core": ["P2ry12", "Tmem119", "Cx3cr1", "Hexb", "Csf1r"],
    "BAM/CAM": ["Pf4", "Mrc1", "Lyve1", "Cd163", "F13a1"],
    "Choroid-plexus-associated": ["Ttr"],
    "Oligodendrocyte lineage": ["Plp1", "Mbp", "Mal", "Cldn11", "Mobp"],
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_true(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dense_vector(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def package_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-root", type=Path, default=None)
    parser.add_argument(
        "--project-root",
        type=Path,
        required=True,
        help="Project root containing data/processed and data/raw",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    figure_root = (
        args.figure_root.resolve()
        if args.figure_root is not None
        else script_path.parent.parent.resolve()
    )
    project_root = args.project_root.resolve()
    panel_dir = figure_root / "data" / "panel_ready"
    table_dir = figure_root / "outputs" / "tableS1"
    manifest_dir = figure_root / "manifests"
    for directory in (panel_dir, table_dir, manifest_dir):
        directory.mkdir(parents=True, exist_ok=True)

    mg_path = project_root / "data" / "processed" / "adata_microglia_pseudotime.h5ad"
    raw_path = project_root / "data" / "processed" / "adata_raw.h5ad"
    doublet_path = (
        project_root
        / "data"
        / "results"
        / "scrublet"
        / "microglia_doublet_summary.csv"
    )
    input_paths = [mg_path, raw_path, doublet_path]
    assert_true(all(path.exists() for path in input_paths), "A required input is missing.")

    mg = ad.read_h5ad(mg_path)
    required_obs = {"sample", "group", "mg_leiden", "mg_subtype2"}
    assert_true(required_obs.issubset(mg.obs.columns), "Submitted partition metadata are missing.")
    assert_true(mg.n_obs == 7461, "Expected 7,461 submitted microglial cells.")

    obs = mg.obs.loc[:, ["sample", "group", "mg_leiden", "mg_subtype2"]].copy()
    obs["partition"] = obs["mg_leiden"].astype(str)
    obs["traceability_name"] = obs["partition"].map(PARTITION_MAP)
    assert_true(obs["traceability_name"].notna().all(), "Unexpected submitted partition label.")
    observed_mapping = (
        obs.loc[:, ["partition", "mg_subtype2"]]
        .drop_duplicates()
        .set_index("partition")["mg_subtype2"]
        .astype(str)
        .to_dict()
    )
    assert_true(observed_mapping == PARTITION_MAP, "Numeric-to-traceability mapping changed.")

    raw = ad.read_h5ad(raw_path)
    assert_true(set(obs.index).issubset(raw.obs_names), "Microglial cell IDs do not map to raw counts.")
    counts = raw[obs.index].copy()
    symbols = counts.var["gene_symbol"].astype(str).to_numpy()
    valid_symbol = pd.Series(symbols).notna().to_numpy() & (symbols != "") & (symbols != "nan")
    first_symbol = ~pd.Series(symbols).duplicated(keep="first").to_numpy()
    counts = counts[:, valid_symbol & first_symbol].copy()
    counts.var_names = symbols[valid_symbol & first_symbol]
    detected = np.asarray((counts.X > 0).sum(axis=0)).ravel() > 0
    counts = counts[:, detected].copy()
    counts.obs = obs.copy()
    counts.obs["partition"] = pd.Categorical(
        counts.obs["partition"], categories=list(PARTITION_MAP), ordered=True
    )
    # All seven submitted partitions retain 17,878 detected symbols.  The
    # corresponding non-Rare set contains 17,799; the latter is not the
    # estimand for this all-partition traceability audit.
    assert_true(counts.n_vars == 17878, "Detected symbol-deduplicated universe changed.")

    raw_detected = counts.X.copy()
    raw_totals = np.asarray(raw_detected.sum(axis=1)).ravel()
    assert_true(np.all(raw_totals > 0), "A submitted microglial cell has zero total counts.")

    # Panel A: author-selected descriptive genes on log1p(total-count/10,000).
    missing_selected = sorted(set(SELECTED_GENES) - set(counts.var_names))
    assert_true(not missing_selected, f"Selected genes missing: {missing_selected}")
    normalized = counts.copy()
    sc.pp.normalize_total(normalized, target_sum=1e4)
    sc.pp.log1p(normalized)

    selected_rows: list[dict] = []
    for partition, trace_name in PARTITION_MAP.items():
        mask = np.asarray(normalized.obs["partition"].astype(str) == partition)
        for gene_order, gene in enumerate(SELECTED_GENES, start=1):
            raw_values = dense_vector(counts[mask, gene].X)
            log_values = dense_vector(normalized[mask, gene].X)
            selected_rows.append(
                {
                    "partition": int(partition),
                    "traceability_name": trace_name,
                    "gene": gene,
                    "gene_order": gene_order,
                    "n_cells": int(mask.sum()),
                    "mean_log1p_cp10k": float(log_values.mean()),
                    "pct_detected": float(100 * np.mean(raw_values > 0)),
                    "selection_status": "author-selected descriptive panel",
                }
            )
    selected_df = pd.DataFrame(selected_rows)

    # Panel B and Table S1A: submitted partition composition and Scrublet calls.
    doublets = pd.read_csv(doublet_path)
    summary = (
        obs.groupby(["partition", "traceability_name"], observed=True)
        .agg(
            n_cells=("sample", "size"),
            n_oxygen_control=("group", lambda values: int(np.sum(values == "Control"))),
            n_combined_exposure=("group", lambda values: int(np.sum(values == "Surgery"))),
        )
        .reset_index()
    )
    summary["partition"] = summary["partition"].astype(int)
    doublets = doublets.rename(columns={"subtype": "traceability_name"})
    summary = summary.merge(
        doublets[
            [
                "traceability_name",
                "n_cells",
                "n_doublets",
                "doublet_rate",
                "control_doublet_rate",
                "surgery_doublet_rate",
            ]
        ],
        on=["traceability_name", "n_cells"],
        how="left",
        validate="one_to_one",
    ).sort_values("partition")
    assert_true(summary["n_doublets"].sum() == 90, "Expected 90 Scrublet-predicted doublets.")
    assert_true(
        int(summary.loc[summary["partition"] == 6, "n_doublets"].iloc[0]) == 82,
        "Expected 82/94 Scrublet calls in partition 6.",
    )
    summary["doublet_rate_percent"] = 100 * summary["doublet_rate"]
    summary["interpretation_ceiling"] = np.where(
        summary["partition"] == 6,
        "Doublet-enriched; excluded from biological interpretation",
        "Submitted traceability partition; not a validated biological state",
    )

    # Panel C: raw-count detection of identity, contamination and ambient markers.
    detection_rows: list[dict] = []
    for partition, trace_name in PARTITION_MAP.items():
        mask = np.asarray(counts.obs["partition"].astype(str) == partition)
        for group_order, (gene_group, genes) in enumerate(
            DETECTION_GENE_GROUPS.items(), start=1
        ):
            for within_order, gene in enumerate(genes, start=1):
                assert_true(gene in counts.var_names, f"Detection gene missing: {gene}")
                values = dense_vector(counts[mask, gene].X)
                detection_rows.append(
                    {
                        "partition": int(partition),
                        "traceability_name": trace_name,
                        "gene_group": gene_group,
                        "gene_group_order": group_order,
                        "gene": gene,
                        "within_group_order": within_order,
                        "n_cells": int(mask.sum()),
                        "pct_detected": float(100 * np.mean(values > 0)),
                    }
                )
    detection_df = pd.DataFrame(detection_rows)

    # Table S1B: complete positive one-versus-rest descriptive marker audit.
    sc.tl.rank_genes_groups(
        normalized,
        groupby="partition",
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        use_raw=False,
        n_genes=normalized.n_vars,
        pts=True,
        tie_correct=True,
    )
    marker_frames: list[pd.DataFrame] = []
    for partition, trace_name in PARTITION_MAP.items():
        frame = sc.get.rank_genes_groups_df(normalized, group=partition).copy()
        frame = frame.rename(
            columns={
                "names": "gene",
                "scores": "wilcoxon_score",
                "logfoldchanges": "log2_fold_change_scanpy",
                "pvals": "p_value",
                "pvals_adj": "adjusted_p_value_bh",
                "pct_nz_group": "fraction_detected_partition",
                "pct_nz_reference": "fraction_detected_rest",
            }
        )
        expected_columns = {
            "gene",
            "wilcoxon_score",
            "log2_fold_change_scanpy",
            "p_value",
            "adjusted_p_value_bh",
            "fraction_detected_partition",
            "fraction_detected_rest",
        }
        assert_true(expected_columns.issubset(frame.columns), "Marker output schema changed.")
        frame.insert(0, "rank_all", np.arange(1, len(frame) + 1))
        frame = frame.loc[
            (frame["adjusted_p_value_bh"] < 0.05)
            & (frame["log2_fold_change_scanpy"] > 0)
        ].copy()
        frame.insert(0, "positive_marker_rank", np.arange(1, len(frame) + 1))
        frame.insert(0, "traceability_name", trace_name)
        frame.insert(0, "partition", int(partition))
        frame["pct_detected_partition"] = 100 * frame["fraction_detected_partition"]
        frame["pct_detected_rest"] = 100 * frame["fraction_detected_rest"]
        frame["top_50_marker"] = frame["positive_marker_rank"] <= 50
        frame["statistical_unit"] = "cell-level descriptive ranking"
        marker_frames.append(frame)
    markers = pd.concat(marker_frames, ignore_index=True)
    assert_true(markers["partition"].nunique() == 7, "Not all partitions have marker output.")

    metadata_df = pd.DataFrame(
        [
            {
                "n_submitted_cells": int(counts.n_obs),
                "n_submitted_partitions": int(counts.obs["partition"].nunique()),
                "n_detected_symbol_deduplicated_genes": int(counts.n_vars),
                "normalization": "total-count normalization to 10,000 followed by log1p",
                "marker_method": "one-versus-rest Wilcoxon rank-sum with tie correction; BH within partition",
                "marker_retention": "adjusted P < 0.05 and Scanpy log2 fold change > 0",
                "claim_ceiling": "retrospective traceability audit; no biological-state validation or animal-level inference",
                "published_signature_comparison": "withdrawn",
                "published_signature_reason": "submitted partition is seed-sensitive and annotation rule was not reproducible",
            }
        ]
    )

    output_paths = {
        "metadata": panel_dir / "figS01_metadata.csv",
        "selected_gene_dotplot": panel_dir / "figS01_selected_gene_dotplot.csv",
        "doublet_summary": panel_dir / "figS01_doublet_summary.csv",
        "detection_heatmap": panel_dir / "figS01_detection_heatmap.csv",
        "tableS1_partition_summary_csv": table_dir / "TableS1_partA_partition_summary.csv",
        "tableS1_partition_summary_tsv": table_dir / "TableS1_partA_partition_summary.tsv",
        "tableS1_marker_audit_csv_gz": table_dir / "TableS1_partB_marker_audit.csv.gz",
        "tableS1_marker_audit_tsv_gz": table_dir / "TableS1_partB_marker_audit.tsv.gz",
    }
    metadata_df.to_csv(output_paths["metadata"], index=False)
    selected_df.to_csv(output_paths["selected_gene_dotplot"], index=False)
    summary.to_csv(output_paths["doublet_summary"], index=False)
    detection_df.to_csv(output_paths["detection_heatmap"], index=False)
    summary.to_csv(output_paths["tableS1_partition_summary_csv"], index=False)
    summary.to_csv(output_paths["tableS1_partition_summary_tsv"], sep="\t", index=False)
    markers.to_csv(output_paths["tableS1_marker_audit_csv_gz"], index=False, compression="gzip")
    markers.to_csv(
        output_paths["tableS1_marker_audit_tsv_gz"],
        sep="\t",
        index=False,
        compression="gzip",
    )

    source_manifest = pd.DataFrame(
        [
            {
                "role": "source",
                "path": str(path),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in input_paths
        ]
        + [
            {
                "role": "script",
                "path": str(script_path),
                "sha256": sha256_file(script_path),
                "size_bytes": script_path.stat().st_size,
            }
        ]
    )
    source_manifest.to_csv(manifest_dir / "figS01_source_manifest.csv", index=False)

    output_manifest = pd.DataFrame(
        [
            {
                "role": key,
                "path": str(path),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for key, path in output_paths.items()
        ]
    )
    output_manifest.to_csv(manifest_dir / "figS01_data_output_manifest.csv", index=False)

    execution = {
        "analysis": "Figure S1 and Table S1 retrospective submitted-partition audit",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(script_path),
        "script_sha256": sha256_file(script_path),
        "python": platform.python_version(),
        "packages": {
            name: package_version(name)
            for name in ["anndata", "numpy", "pandas", "scanpy", "scipy"]
        },
        "anchors": {
            "n_cells": int(counts.n_obs),
            "n_partitions": int(counts.obs["partition"].nunique()),
            "n_detected_symbol_deduplicated_genes": int(counts.n_vars),
            "n_scrublet_calls": int(summary["n_doublets"].sum()),
            "partition_6_scrublet_calls": int(
                summary.loc[summary["partition"] == 6, "n_doublets"].iloc[0]
            ),
            "partition_6_cells": int(
                summary.loc[summary["partition"] == 6, "n_cells"].iloc[0]
            ),
            "n_retained_positive_markers": int(len(markers)),
        },
        "claim_ceiling": metadata_df.loc[0, "claim_ceiling"],
    }
    (manifest_dir / "figS01_data_preparation_execution.json").write_text(
        json.dumps(execution, indent=2) + "\n"
    )

    print(json.dumps(execution["anchors"], indent=2))


if __name__ == "__main__":
    main()
