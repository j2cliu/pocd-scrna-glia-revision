#!/usr/bin/env python3
"""Sensitivity of Scanpy score_genes to its matched-control specification.

The revised analysis treated a sign reversal after dropping S3 as evidence that
the IFN signal was a single-animal artefact. This script checks whether that
conclusion is specific to Scanpy's randomly sampled control-gene subtraction.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


PROJECT = Path(os.environ["POCD_SCRNA_PROJECT_ROOT"]).expanduser().resolve()
OUT_DIR = Path(__file__).resolve().parent
CORE = ["Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3"]
BROAD = [
    "Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3",
    "Oasl2", "Oas1a", "Oas2", "Rsad2", "Bst2", "Stat1", "Stat2", "Irf9",
]


def cohens_d(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    pooled = np.sqrt(
        ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
        / (len(a) + len(b) - 2)
    )
    return float((a.mean() - b.mean()) / pooled) if pooled > 0 else np.nan


def exact_perm_p(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    pooled = np.concatenate([a, b])
    observed = abs(a.mean() - b.mean())
    hits = total = 0
    for idx in combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), bool)
        mask[list(idx)] = True
        hits += abs(pooled[mask].mean() - pooled[~mask].mean()) >= observed - 1e-12
        total += 1
    return hits / total


def main():
    raw = ad.read_h5ad(PROJECT / "data/processed/adata_raw.h5ad")
    final = ad.read_h5ad(PROJECT / "data/processed/adata_microglia_final.h5ad")
    a = raw[final.obs_names].copy()
    symbols = a.var["gene_symbol"].astype(str).to_numpy()
    keep = ~pd.Series(symbols).duplicated(keep="first").to_numpy()
    a = a[:, keep].copy()
    a.var_names = symbols[keep]
    a.var_names_make_unique()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.obs["sample"] = final.obs["sample"].astype(str).to_numpy()
    a.obs["group"] = final.obs["group"].astype(str).to_numpy()
    a.obs["subtype"] = final.obs["mg_subtype2"].astype(str).to_numpy()
    a = a[a.obs["subtype"] != "Rare"].copy()

    surgery_ids = ["S1", "S2", "S3"]
    control_ids = ["C1", "C2", "C3"]
    configs = []
    for seed in range(20):
        configs.append(("seed_sweep", seed, 50, 25))
    for ctrl_size in (25, 50, 100, 200):
        configs.append(("control_size_sweep", 42, ctrl_size, 25))
    for n_bins in (10, 25, 50):
        configs.append(("bin_sweep", 42, 50, n_bins))

    rows = []
    for panel_name, genes in (("core_isg7", CORE), ("broad_isg15", BROAD)):
        for sweep, seed, ctrl_size, n_bins in configs:
            key = f"score_{panel_name}_{sweep}_{seed}_{ctrl_size}_{n_bins}"
            sc.tl.score_genes(
                a,
                gene_list=genes,
                score_name=key,
                ctrl_size=ctrl_size,
                n_bins=n_bins,
                random_state=seed,
                use_raw=False,
            )
            per_animal = a.obs.groupby("sample", observed=True)[key].mean()
            surgery = per_animal.loc[surgery_ids].to_numpy()
            control = per_animal.loc[control_ids].to_numpy()
            no_s3 = per_animal.loc[["S1", "S2"]].to_numpy()
            full_diff = float(surgery.mean() - control.mean())
            no_s3_diff = float(no_s3.mean() - control.mean())
            rows.append(
                {
                    "panel": panel_name,
                    "sweep": sweep,
                    "seed": seed,
                    "ctrl_size": ctrl_size,
                    "n_bins": n_bins,
                    "mean_difference": full_diff,
                    "cohens_d": cohens_d(surgery, control),
                    "exact_perm_p": exact_perm_p(surgery, control),
                    "mean_difference_without_s3": no_s3_diff,
                    "cohens_d_without_s3": cohens_d(no_s3, control),
                    "sign_preserved_without_s3": bool(
                        np.sign(no_s3_diff) == np.sign(full_diff)
                    ),
                    "surgery_values": ";".join(
                        f"{s}:{per_animal.loc[s]:.8g}" for s in surgery_ids
                    ),
                    "control_values": ";".join(
                        f"{s}:{per_animal.loc[s]:.8g}" for s in control_ids
                    ),
                }
            )
            del a.obs[key]

    result = pd.DataFrame(rows)
    result.to_csv(OUT_DIR / "score_genes_sensitivity.csv", index=False)

    summary = (
        result.groupby(["panel", "sweep"], observed=True)
        .agg(
            n_runs=("cohens_d", "size"),
            d_min=("cohens_d", "min"),
            d_median=("cohens_d", "median"),
            d_max=("cohens_d", "max"),
            diff_no_s3_min=("mean_difference_without_s3", "min"),
            diff_no_s3_max=("mean_difference_without_s3", "max"),
            n_sign_preserved_without_s3=("sign_preserved_without_s3", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "score_genes_sensitivity_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
