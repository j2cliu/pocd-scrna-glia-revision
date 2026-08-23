#!/usr/bin/env python3
"""Independent sample-level audit of the revised POCD IFN/ISG claim.

This script intentionally does not use Scanpy's score_genes. It reconstructs
animal-level summaries directly from the raw count matrix and reports:

1. pseudobulk log2-CPM for each prespecified gene;
2. mean per-cell log1p(CP10K) expression for each gene;
3. unweighted panel summaries at the animal level;
4. leave-one-animal-out mean differences and Cohen's d;
5. basic sample-level QC quantities.

It reads the project files but writes only to this independent audit workspace.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT = Path(os.environ["POCD_SCRNA_PROJECT_ROOT"]).expanduser().resolve()
RAW_PATH = PROJECT / "data/processed/adata_raw.h5ad"
FINAL_PATH = PROJECT / "data/processed/adata_microglia_final.h5ad"
OUT_DIR = Path(__file__).resolve().parent

PANELS = {
    "core_isg7": ["Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3"],
    "broad_isg15": [
        "Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3",
        "Oasl2", "Oas1a", "Oas2", "Rsad2", "Bst2", "Stat1", "Stat2", "Irf9",
    ],
    "holtman_isg3": ["Irf7", "Ifitm3", "Isg15"],
    "holtman_nonisg4": ["Il1b", "Lgals3", "Cybb", "Fcer1g"],
    "ieg13": [
        "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Ier2",
        "Dusp1", "Zfp36", "Hspa1a", "Hspa1b", "Socs3",
    ],
}


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    pooled = np.sqrt(
        ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
        / (len(a) + len(b) - 2)
    )
    return float((a.mean() - b.mean()) / pooled) if pooled > 0 else np.nan


def exact_perm_p(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    pooled = np.concatenate([a, b])
    observed = abs(a.mean() - b.mean())
    hits = 0
    total = 0
    for idx in combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        candidate = abs(pooled[mask].mean() - pooled[~mask].mean())
        hits += candidate >= observed - 1e-12
        total += 1
    return hits / total


def to_1d(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def summarize_panel(
    values: pd.DataFrame, panel: str, metric: str, group_map: dict[str, str]
) -> tuple[dict, list[dict]]:
    """Summarize an animal x gene table using an equal-weight gene mean."""
    genes = [g for g in PANELS[panel] if g in values.columns]
    score = values[genes].mean(axis=1)
    surgery_ids = [s for s in score.index if group_map[s] == "Surgery"]
    control_ids = [s for s in score.index if group_map[s] == "Control"]
    surgery = score.loc[surgery_ids].to_numpy()
    control = score.loc[control_ids].to_numpy()
    full_diff = float(surgery.mean() - control.mean())
    full_d = cohens_d(surgery, control)
    row = {
        "metric": metric,
        "panel": panel,
        "n_genes": len(genes),
        "genes": ";".join(genes),
        "surgery_mean": float(surgery.mean()),
        "control_mean": float(control.mean()),
        "mean_difference": full_diff,
        "cohens_d": full_d,
        "exact_perm_p": exact_perm_p(surgery, control),
        "surgery_values": ";".join(f"{s}:{score.loc[s]:.8g}" for s in surgery_ids),
        "control_values": ";".join(f"{s}:{score.loc[s]:.8g}" for s in control_ids),
    }
    loo_rows = []
    for dropped in score.index:
        kept_s = [s for s in surgery_ids if s != dropped]
        kept_c = [s for s in control_ids if s != dropped]
        a = score.loc[kept_s].to_numpy()
        b = score.loc[kept_c].to_numpy()
        loo_rows.append(
            {
                "metric": metric,
                "panel": panel,
                "dropped_animal": dropped,
                "n_surgery": len(a),
                "n_control": len(b),
                "mean_difference": float(a.mean() - b.mean()),
                "cohens_d": cohens_d(a, b),
                "sign_matches_full": bool(
                    np.sign(a.mean() - b.mean()) == np.sign(full_diff)
                ),
            }
        )
    return row, loo_rows


def main() -> None:
    raw = ad.read_h5ad(RAW_PATH)
    final = ad.read_h5ad(FINAL_PATH)
    if not final.obs_names.is_unique:
        raise RuntimeError("Microglia cell identifiers are not unique")
    missing_cells = final.obs_names.difference(raw.obs_names)
    if len(missing_cells):
        raise RuntimeError(f"{len(missing_cells)} microglia cells absent from raw object")

    mg = raw[final.obs_names].copy()
    mg.obs["sample"] = final.obs["sample"].astype(str).to_numpy()
    mg.obs["group"] = final.obs["group"].astype(str).to_numpy()
    mg.obs["subtype"] = final.obs["mg_subtype2"].astype(str).to_numpy()

    gene_symbols = mg.var["gene_symbol"].astype(str).to_numpy()
    if sparse.issparse(mg.X):
        counts = mg.X.tocsr()
    else:
        counts = sparse.csr_matrix(np.asarray(mg.X))

    integer_residual = float(np.max(np.abs(counts.data - np.rint(counts.data))))
    if integer_residual > 1e-6:
        raise RuntimeError(
            f"adata_raw.X does not appear to contain integer counts "
            f"(max residual={integer_residual})"
        )

    sample_ids = sorted(mg.obs["sample"].unique())
    group_map = (
        mg.obs[["sample", "group"]]
        .drop_duplicates()
        .set_index("sample")["group"]
        .to_dict()
    )
    if any(s.startswith("S") != (group_map[s] == "Surgery") for s in sample_ids):
        raise RuntimeError(f"Sample prefix/group mismatch: {group_map}")

    # Compare the analysis with and without the 94-cell Rare label. The primary
    # audit excludes Rare to match the revised scripts.
    masks = {
        "rare_excluded": mg.obs["subtype"].to_numpy() != "Rare",
        "all_microglia_labels": np.ones(mg.n_obs, dtype=bool),
    }

    all_genes = sorted({g for genes in PANELS.values() for g in genes})
    gene_to_cols: dict[str, list[int]] = {}
    for idx, symbol in enumerate(gene_symbols):
        if symbol in all_genes:
            gene_to_cols.setdefault(symbol, []).append(idx)

    manifest = {
        "raw_path": str(RAW_PATH),
        "final_path": str(FINAL_PATH),
        "raw_shape": list(raw.shape),
        "microglia_shape": list(mg.shape),
        "samples": group_map,
        "integer_count_max_residual": integer_residual,
        "genes_with_multiple_features": {
            g: len(cols) for g, cols in gene_to_cols.items() if len(cols) > 1
        },
        "missing_panel_genes": [g for g in all_genes if g not in gene_to_cols],
    }

    qc_rows: list[dict] = []
    gene_rows: list[dict] = []
    panel_rows: list[dict] = []
    loo_rows: list[dict] = []
    dominance_rows: list[dict] = []
    gene_effect_rows: list[dict] = []

    mt_cols = np.array(
        [i for i, symbol in enumerate(gene_symbols) if symbol.lower().startswith("mt-")],
        dtype=int,
    )

    for mask_name, analysis_mask in masks.items():
        obs = mg.obs.iloc[np.flatnonzero(analysis_mask)]
        x = counts[analysis_mask]
        cell_totals = to_1d(x.sum(axis=1))
        cell_detected = to_1d((x > 0).sum(axis=1))
        mt_counts = (
            to_1d(x[:, mt_cols].sum(axis=1))
            if len(mt_cols)
            else np.zeros(x.shape[0], dtype=float)
        )

        per_gene_logcpm: dict[str, dict[str, float]] = {
            g: {} for g in gene_to_cols
        }
        per_gene_cellmean: dict[str, dict[str, float]] = {
            g: {} for g in gene_to_cols
        }

        for sample in sample_ids:
            local_mask = obs["sample"].to_numpy() == sample
            sx = x[local_mask]
            stotal = float(sx.sum())
            qc_rows.append(
                {
                    "mask": mask_name,
                    "sample": sample,
                    "group": group_map[sample],
                    "n_cells": int(local_mask.sum()),
                    "pseudobulk_total_umi": stotal,
                    "median_umi_per_cell": float(np.median(cell_totals[local_mask])),
                    "median_genes_per_cell": float(np.median(cell_detected[local_mask])),
                    "mean_mt_fraction": float(
                        np.mean(
                            np.divide(
                                mt_counts[local_mask],
                                cell_totals[local_mask],
                                out=np.zeros_like(mt_counts[local_mask], dtype=float),
                                where=cell_totals[local_mask] > 0,
                            )
                        )
                    ),
                }
            )

            for gene, cols in gene_to_cols.items():
                gene_counts = to_1d(sx[:, cols].sum(axis=1))
                pseudobulk_count = float(gene_counts.sum())
                logcpm = float(
                    np.log2((pseudobulk_count + 0.5) / (stotal + 1.0) * 1e6)
                )
                norm_cell = np.log1p(
                    gene_counts
                    * np.divide(
                        1e4,
                        cell_totals[local_mask],
                        out=np.zeros_like(cell_totals[local_mask], dtype=float),
                        where=cell_totals[local_mask] > 0,
                    )
                )
                cellmean = float(norm_cell.mean())
                per_gene_logcpm[gene][sample] = logcpm
                per_gene_cellmean[gene][sample] = cellmean
                gene_rows.append(
                    {
                        "mask": mask_name,
                        "sample": sample,
                        "group": group_map[sample],
                        "gene": gene,
                        "n_features_collapsed": len(cols),
                        "pseudobulk_count": pseudobulk_count,
                        "pseudobulk_log2cpm": logcpm,
                        "mean_cell_log1p_cp10k": cellmean,
                        "pct_cells_detected": float((gene_counts > 0).mean() * 100),
                    }
                )

        for metric, source in (
            ("pseudobulk_log2cpm", per_gene_logcpm),
            ("mean_cell_log1p_cp10k", per_gene_cellmean),
        ):
            table = pd.DataFrame(source).loc[sample_ids]
            surgery_ids = [s for s in sample_ids if group_map[s] == "Surgery"]
            control_ids = [s for s in sample_ids if group_map[s] == "Control"]
            for gene in table.columns:
                surgery = table.loc[surgery_ids, gene].to_numpy()
                control = table.loc[control_ids, gene].to_numpy()
                surgery_without_s3 = table.loc[
                    [s for s in surgery_ids if s != "S3"], gene
                ].to_numpy()
                gene_effect_rows.append(
                    {
                        "mask": mask_name,
                        "metric": metric,
                        "gene": gene,
                        "mean_difference": float(surgery.mean() - control.mean()),
                        "cohens_d": cohens_d(surgery, control),
                        "exact_perm_p": exact_perm_p(surgery, control),
                        "perfect_separation": bool(
                            surgery.min() > control.max()
                            or surgery.max() < control.min()
                        ),
                        "mean_difference_without_s3": float(
                            surgery_without_s3.mean() - control.mean()
                        ),
                        "n_surgery_above_control_max": int(
                            (surgery > control.max()).sum()
                        ),
                        "n_surgery_above_control_mean": int(
                            (surgery > control.mean()).sum()
                        ),
                    }
                )
            for panel in PANELS:
                row, rows_loo = summarize_panel(table, panel, metric, group_map)
                row["mask"] = mask_name
                panel_rows.append(row)
                for r in rows_loo:
                    r["mask"] = mask_name
                    loo_rows.append(r)

            if metric == "pseudobulk_log2cpm":
                for panel, genes in PANELS.items():
                    used = [g for g in genes if g in table.columns]
                    control_mean = table.loc[
                        [s for s in sample_ids if group_map[s] == "Control"], used
                    ].mean()
                    for sample in sample_ids:
                        deltas = table.loc[sample, used] - control_mean
                        dominance_rows.append(
                            {
                                "mask": mask_name,
                                "panel": panel,
                                "sample": sample,
                                "group": group_map[sample],
                                "n_genes_above_control_mean": int((deltas > 0).sum()),
                                "n_genes_below_control_mean": int((deltas < 0).sum()),
                                "mean_gene_delta_vs_control_mean": float(deltas.mean()),
                            }
                        )

    qc = pd.DataFrame(qc_rows)
    genes = pd.DataFrame(gene_rows)
    panels = pd.DataFrame(panel_rows)
    loo = pd.DataFrame(loo_rows)
    dominance = pd.DataFrame(dominance_rows)
    gene_effects = pd.DataFrame(gene_effect_rows)

    qc.to_csv(OUT_DIR / "sample_qc.csv", index=False)
    genes.to_csv(OUT_DIR / "per_sample_gene_values.csv", index=False)
    panels.to_csv(OUT_DIR / "panel_effects.csv", index=False)
    loo.to_csv(OUT_DIR / "panel_leave_one_animal_out.csv", index=False)
    dominance.to_csv(OUT_DIR / "panel_gene_direction_by_animal.csv", index=False)
    gene_effects.to_csv(OUT_DIR / "gene_effects.csv", index=False)
    (OUT_DIR / "input_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print("INPUT")
    print(json.dumps(manifest, indent=2))
    print("\nPRIMARY PANEL EFFECTS (Rare excluded)")
    print(
        panels[
            (panels["mask"] == "rare_excluded")
            & (panels["panel"].isin(["core_isg7", "broad_isg15", "ieg13"]))
        ][
            [
                "metric",
                "panel",
                "mean_difference",
                "cohens_d",
                "exact_perm_p",
                "surgery_values",
                "control_values",
            ]
        ].to_string(index=False)
    )
    print("\nCORE ISG7 GENE DIRECTIONS BY ANIMAL (pseudobulk log2-CPM)")
    print(
        dominance[
            (dominance["mask"] == "rare_excluded")
            & (dominance["panel"] == "core_isg7")
        ].to_string(index=False)
    )
    print("\nCORE ISG7 LOO (pseudobulk log2-CPM)")
    print(
        loo[
            (loo["mask"] == "rare_excluded")
            & (loo["metric"] == "pseudobulk_log2cpm")
            & (loo["panel"] == "core_isg7")
        ].to_string(index=False)
    )
    print("\nCORE ISG7 PER-GENE EFFECTS (pseudobulk log2-CPM)")
    print(
        gene_effects[
            (gene_effects["mask"] == "rare_excluded")
            & (gene_effects["metric"] == "pseudobulk_log2cpm")
            & (gene_effects["gene"].isin(PANELS["core_isg7"]))
        ].sort_values("gene").to_string(index=False)
    )


if __name__ == "__main__":
    main()
