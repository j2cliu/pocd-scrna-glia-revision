#!/usr/bin/env python3
"""Audit the graph inputs and claims used in the Claude cluster-stability pass."""

from __future__ import annotations

import itertools
from pathlib import Path
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score


PROJECT = Path(os.environ["POCD_SCRNA_PROJECT_ROOT"]).expanduser().resolve()
OUT_DIR = Path(__file__).resolve().parent
SEEDS = [0, 1, 7, 42, 123, 2024]


def graph_summary(name, obj):
    graph = obj.obsp["connectivities"].tocsr()
    return {
        "object": name,
        "n_obs": obj.n_obs,
        "nnz": graph.nnz,
        "sum_weights": float(graph.sum()),
        "neighbor_params": repr(obj.uns.get("neighbors", {}).get("params", {})),
    }


def main():
    scored = ad.read_h5ad(PROJECT / "data/processed/adata_microglia_scored.h5ad")
    subtyped = ad.read_h5ad(PROJECT / "data/processed/adata_microglia_subtyped.h5ad")
    final = ad.read_h5ad(PROJECT / "data/processed/adata_microglia_final.h5ad")
    if not np.array_equal(scored.obs_names, subtyped.obs_names):
        raise RuntimeError("Cell order differs")
    if not np.array_equal(scored.obs_names, final.obs_names):
        raise RuntimeError("Final object cell order differs")

    summaries = pd.DataFrame(
        [graph_summary("scored_pre_subtyping", scored), graph_summary("subtyped_output", subtyped)]
    )
    summaries.to_csv(OUT_DIR / "cluster_graph_inputs.csv", index=False)

    g_scored = scored.obsp["connectivities"].tocsr()
    g_subtyped = subtyped.obsp["connectivities"].tocsr()
    graph_delta = g_scored - g_subtyped
    graph_comparison = {
        "same_shape": g_scored.shape == g_subtyped.shape,
        "same_sparse_values": bool(
            graph_delta.nnz == 0
            or np.allclose(graph_delta.data, 0, rtol=0, atol=1e-12)
        ),
        "n_differing_entries": int(graph_delta.nnz),
        "absolute_weight_difference": float(abs(graph_delta).sum()),
    }
    pd.DataFrame([graph_comparison]).to_csv(
        OUT_DIR / "cluster_graph_comparison.csv", index=False
    )

    reference = subtyped.obs["mg_leiden"].astype(str)
    labels = {}
    seed_rows = []
    for seed in SEEDS:
        key = f"audit_leiden_{seed}"
        sc.tl.leiden(
            subtyped,
            resolution=0.4,
            random_state=seed,
            key_added=key,
        )
        labels[seed] = subtyped.obs[key].astype(str).copy()
        seed_rows.append(
            {
                "seed": seed,
                "n_clusters": labels[seed].nunique(),
                "ari_vs_original": adjusted_rand_score(reference, labels[seed]),
            }
        )
    seed_df = pd.DataFrame(seed_rows)
    pair_rows = []
    for a, b in itertools.combinations(SEEDS, 2):
        pair_rows.append(
            {"seed_a": a, "seed_b": b, "ari": adjusted_rand_score(labels[a], labels[b])}
        )
    pair_df = pd.DataFrame(pair_rows)
    seed_df.to_csv(OUT_DIR / "subtyped_stored_graph_seed_stability.csv", index=False)
    pair_df.to_csv(OUT_DIR / "subtyped_stored_graph_pairwise_ari.csv", index=False)

    # Reconstruct the current-version graph used by script 74 and map the stable
    # resolution-0.5/0.6 clusters back to the submitted labels and animals.
    base = scored.copy()
    base.obsm["X_pca"] = base.obsm["X_pca"][:, :20].copy()
    sc.pp.neighbors(base, n_neighbors=15, n_pcs=20, random_state=0)
    base.obs["submitted_subtype"] = final.obs["mg_subtype2"].astype(str).to_numpy()
    base.obs["sample"] = final.obs["sample"].astype(str).to_numpy()
    mapping_rows = []
    for resolution in (0.5, 0.6):
        key = f"res_{resolution}"
        sc.tl.leiden(base, resolution=resolution, random_state=0, key_added=key)
        stability = pd.read_csv(
            PROJECT / "data/results/cluster_stability/resolution_sweep_percluster.csv"
        )
        stable_ids = stability[
            (stability["resolution"] == resolution)
            & (stability["jaccard_mean"] >= 0.75)
        ]["cluster"].astype(str)
        for cluster_id in stable_ids:
            mask = base.obs[key].astype(str) == cluster_id
            submitted = (
                base.obs.loc[mask, "submitted_subtype"]
                .value_counts(normalize=True)
                .to_dict()
            )
            animals = base.obs.loc[mask, "sample"].value_counts().to_dict()
            mapping_rows.append(
                {
                    "resolution": resolution,
                    "cluster": cluster_id,
                    "n_cells": int(mask.sum()),
                    "submitted_subtype_mix": ";".join(
                        f"{k}:{v:.3f}" for k, v in submitted.items()
                    ),
                    "animal_counts": ";".join(f"{k}:{v}" for k, v in sorted(animals.items())),
                }
            )
    mapping = pd.DataFrame(mapping_rows)
    mapping.to_csv(OUT_DIR / "stable_cluster_mapping.csv", index=False)

    print("GRAPH INPUTS")
    print(summaries.to_string(index=False))
    print("\nGRAPH COMPARISON")
    print(pd.DataFrame([graph_comparison]).to_string(index=False))
    print("\nCORRECT STORED GRAPH: LEIDEN SEED STABILITY")
    print(seed_df.to_string(index=False))
    print(
        f"pairwise ARI mean={pair_df['ari'].mean():.3f}, "
        f"min={pair_df['ari'].min():.3f}, max={pair_df['ari'].max():.3f}"
    )
    print("\nSTABLE CLUSTERS REPORTED BY SCRIPT 74, MAPPED TO SUBMITTED LABELS")
    print(mapping.to_string(index=False))


if __name__ == "__main__":
    main()
