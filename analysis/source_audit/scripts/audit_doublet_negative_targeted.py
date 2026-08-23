#!/usr/bin/env python3
"""Doublet-negative animal-level audit of the ISG7 and IEG13 panels.

This bounded analysis reads, but does not modify, the project AnnData files. It:

1. takes the 7,461 frozen microglial barcodes and ``predicted_doublet`` calls
   from ``adata_microglia_scrublet_annotated.h5ad``;
2. removes all Scrublet-positive cells;
3. retrieves integer UMI counts for the retained barcodes from
   ``adata_raw.h5ad``;
4. computes per-animal gene-level pseudobulk log2-CPM values and equal-weight
   means for the author-selected ISG7 and exploratory IEG13 panels; and
5. reports full 3-vs-3 effects and systematic leave-one-animal-out influence
   diagnostics.

Full-cohort uncertainty is a Welch--Satterthwaite 95% confidence interval for
the exposed-minus-control mean difference. Cohen's d is reported as a point
estimate only. Exact two-sided permutation p values enumerate all C(6, 3) = 20
label assignments. Leave-one-animal-out rows are influence diagnostics and are
not assigned permutation p values.
"""

from __future__ import annotations

from itertools import combinations
from importlib.metadata import version
from pathlib import Path
import os
import platform

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.stats import t


PROJECT = Path(os.environ["POCD_SCRNA_PROJECT_ROOT"]).expanduser().resolve()
SCRUBLET_PATH = (
    PROJECT / "data/processed/adata_microglia_scrublet_annotated.h5ad"
)
RAW_PATH = PROJECT / "data/processed/adata_raw.h5ad"
OUT_DIR = Path(__file__).resolve().parent

ISG7_GENES = ["Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3"]
IEG13_GENES = [
    "Fos",
    "Fosb",
    "Jun",
    "Junb",
    "Jund",
    "Egr1",
    "Atf3",
    "Ier2",
    "Dusp1",
    "Zfp36",
    "Hspa1a",
    "Hspa1b",
    "Socs3",
]
PANELS = {
    "author_selected_isg7": ISG7_GENES,
    "exploratory_ieg13": IEG13_GENES,
}
TARGET_GENES = list(dict.fromkeys(ISG7_GENES + IEG13_GENES))
SOURCE_EXPOSED_LABEL = "Surgery"
SOURCE_CONTROL_LABEL = "Control"


def to_1d(x) -> np.ndarray:
    """Return a dense one-dimensional float array from a matrix-like object."""
    if sparse.issparse(x):
        return np.asarray(x.toarray(), dtype=float).ravel()
    return np.asarray(x, dtype=float).ravel()


def cohens_d(exposed: np.ndarray, control: np.ndarray) -> float:
    """Pooled-SD Cohen's d, oriented exposed minus control."""
    exposed = np.asarray(exposed, dtype=float)
    control = np.asarray(control, dtype=float)
    pooled_variance = (
        (len(exposed) - 1) * exposed.var(ddof=1)
        + (len(control) - 1) * control.var(ddof=1)
    ) / (len(exposed) + len(control) - 2)
    if pooled_variance <= 0:
        return np.nan
    return float((exposed.mean() - control.mean()) / np.sqrt(pooled_variance))


def welch_mean_difference_ci(
    exposed: np.ndarray, control: np.ndarray, confidence: float = 0.95
) -> dict[str, float]:
    """Welch--Satterthwaite CI for an unpaired mean difference."""
    exposed = np.asarray(exposed, dtype=float)
    control = np.asarray(control, dtype=float)
    exposed_component = exposed.var(ddof=1) / len(exposed)
    control_component = control.var(ddof=1) / len(control)
    variance = exposed_component + control_component
    difference = float(exposed.mean() - control.mean())
    if variance <= 0:
        return {
            "mean_difference": difference,
            "welch_se": 0.0,
            "welch_df": np.nan,
            "ci_level": confidence,
            "mean_difference_ci_low": difference,
            "mean_difference_ci_high": difference,
        }
    standard_error = float(np.sqrt(variance))
    denominator = (
        exposed_component**2 / (len(exposed) - 1)
        + control_component**2 / (len(control) - 1)
    )
    degrees_freedom = float(variance**2 / denominator)
    critical_value = float(t.ppf((1.0 + confidence) / 2.0, degrees_freedom))
    return {
        "mean_difference": difference,
        "welch_se": standard_error,
        "welch_df": degrees_freedom,
        "ci_level": confidence,
        "mean_difference_ci_low": difference - critical_value * standard_error,
        "mean_difference_ci_high": difference + critical_value * standard_error,
    }


def exact_two_sided_permutation(
    exposed: np.ndarray, control: np.ndarray
) -> dict[str, float | int]:
    """Enumerate all allocations with the observed exposed-group size."""
    exposed = np.asarray(exposed, dtype=float)
    control = np.asarray(control, dtype=float)
    pooled = np.concatenate([exposed, control])
    observed = abs(float(exposed.mean() - control.mean()))
    hits = 0
    total = 0
    for exposed_idx in combinations(range(len(pooled)), len(exposed)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(exposed_idx)] = True
        candidate = abs(float(pooled[mask].mean() - pooled[~mask].mean()))
        hits += int(candidate >= observed - 1e-12)
        total += 1
    return {
        "exact_perm_hits": hits,
        "exact_perm_allocations": total,
        "exact_perm_p": hits / total,
    }


def effect_summary(
    values: pd.Series,
    group_map: dict[str, str],
    outcome: str,
    outcome_type: str,
    genes: list[str],
) -> dict:
    """Full 3-vs-3 effect summary for one gene or the panel."""
    exposed_ids = [
        sample for sample in values.index
        if group_map[sample] == SOURCE_EXPOSED_LABEL
    ]
    control_ids = [
        sample for sample in values.index
        if group_map[sample] == SOURCE_CONTROL_LABEL
    ]
    exposed = values.loc[exposed_ids].to_numpy(dtype=float)
    control = values.loc[control_ids].to_numpy(dtype=float)
    row = {
        "outcome_type": outcome_type,
        "outcome": outcome,
        "metric": "pseudobulk_log2cpm",
        "n_genes": len(genes),
        "genes": ";".join(genes),
        "n_exposed": len(exposed),
        "n_control": len(control),
        "exposed_mean": float(exposed.mean()),
        "control_mean": float(control.mean()),
        "cohens_d_point_estimate": cohens_d(exposed, control),
        "perfect_group_separation": bool(
            exposed.min() > control.max() or exposed.max() < control.min()
        ),
        "exposed_values": ";".join(
            f"{sample}:{values.loc[sample]:.10g}" for sample in exposed_ids
        ),
        "control_values": ";".join(
            f"{sample}:{values.loc[sample]:.10g}" for sample in control_ids
        ),
    }
    row.update(welch_mean_difference_ci(exposed, control))
    row.update(exact_two_sided_permutation(exposed, control))
    return row


def loo_summaries(
    values: pd.Series,
    group_map: dict[str, str],
    outcome: str,
    outcome_type: str,
    full_difference: float,
) -> list[dict]:
    """Drop each animal exactly once and recompute descriptive effects."""
    rows: list[dict] = []
    for dropped_animal in values.index:
        retained = values.drop(index=dropped_animal)
        exposed_ids = [
            sample for sample in retained.index
            if group_map[sample] == SOURCE_EXPOSED_LABEL
        ]
        control_ids = [
            sample for sample in retained.index
            if group_map[sample] == SOURCE_CONTROL_LABEL
        ]
        exposed = retained.loc[exposed_ids].to_numpy(dtype=float)
        control = retained.loc[control_ids].to_numpy(dtype=float)
        row = {
            "outcome_type": outcome_type,
            "outcome": outcome,
            "metric": "pseudobulk_log2cpm",
            "dropped_animal": dropped_animal,
            "dropped_source_group": group_map[dropped_animal],
            "n_exposed": len(exposed),
            "n_control": len(control),
            "exposed_mean": float(exposed.mean()),
            "control_mean": float(control.mean()),
            "cohens_d_point_estimate": cohens_d(exposed, control),
        }
        row.update(welch_mean_difference_ci(exposed, control))
        row["sign_matches_full"] = bool(
            np.sign(row["mean_difference"]) == np.sign(full_difference)
        )
        rows.append(row)
    return rows


def main() -> None:
    scrublet = ad.read_h5ad(SCRUBLET_PATH)
    raw = ad.read_h5ad(RAW_PATH)

    required_obs = {"sample", "group", "predicted_doublet", "mg_subtype2"}
    missing_obs = required_obs.difference(scrublet.obs.columns)
    if missing_obs:
        raise RuntimeError(f"Missing Scrublet metadata fields: {sorted(missing_obs)}")
    if not scrublet.obs_names.is_unique:
        raise RuntimeError("Scrublet microglial cell identifiers are not unique")
    if scrublet.obs["predicted_doublet"].isna().any():
        raise RuntimeError("predicted_doublet contains missing values")
    if scrublet.obs["predicted_doublet"].dtype != bool:
        raise RuntimeError(
            "predicted_doublet is not Boolean; refusing implicit coercion"
        )

    missing_in_raw = scrublet.obs_names.difference(raw.obs_names)
    if len(missing_in_raw):
        raise RuntimeError(
            f"{len(missing_in_raw)} frozen microglial barcodes are absent from raw"
        )

    original_obs = scrublet.obs[
        ["sample", "group", "predicted_doublet", "mg_subtype2"]
    ].copy()
    original_obs["sample"] = original_obs["sample"].astype(str)
    original_obs["group"] = original_obs["group"].astype(str)
    original_obs["mg_subtype2"] = original_obs["mg_subtype2"].astype(str)
    group_counts = (
        original_obs[["sample", "group"]].drop_duplicates().groupby("sample").size()
    )
    if not (group_counts == 1).all():
        raise RuntimeError("A sample maps to more than one source group")
    group_map = (
        original_obs[["sample", "group"]]
        .drop_duplicates()
        .set_index("sample")["group"]
        .to_dict()
    )
    if set(group_map.values()) != {SOURCE_CONTROL_LABEL, SOURCE_EXPOSED_LABEL}:
        raise RuntimeError(f"Unexpected source groups: {group_map}")
    if list(sorted(group_map)) != ["C1", "C2", "C3", "S1", "S2", "S3"]:
        raise RuntimeError(f"Unexpected sample identifiers: {sorted(group_map)}")

    retained_mask = ~original_obs["predicted_doublet"].to_numpy()
    retained_names = scrublet.obs_names[retained_mask]
    retained_obs = original_obs.loc[retained_names].copy()
    removed_n = int((~retained_mask).sum())
    retained_n = int(retained_mask.sum())
    if scrublet.n_obs != 7461 or removed_n != 90 or retained_n != 7371:
        raise RuntimeError(
            "Unexpected cell counts: "
            f"input={scrublet.n_obs}, removed={removed_n}, retained={retained_n}"
        )
    rare_mask = original_obs["mg_subtype2"].to_numpy() == "Rare"
    rare_total = int(rare_mask.sum())
    rare_doublets = int(
        original_obs.loc[rare_mask, "predicted_doublet"].sum()
    )
    retained_rare = int((rare_mask & retained_mask).sum())
    nonrare_doublets = int((~rare_mask & ~retained_mask).sum())

    microglia = raw[retained_names].copy()
    if sparse.issparse(microglia.X):
        counts = microglia.X.tocsr()
    else:
        counts = sparse.csr_matrix(np.asarray(microglia.X))
    if counts.data.size == 0:
        raise RuntimeError("Retained raw count matrix is empty")
    integer_residual = float(
        np.max(np.abs(counts.data - np.rint(counts.data)))
    )
    if integer_residual > 1e-6 or counts.data.min() < 0:
        raise RuntimeError(
            "adata_raw.X failed nonnegative-integer UMI validation: "
            f"max residual={integer_residual}, min={counts.data.min()}"
        )

    if "gene_symbol" not in microglia.var.columns:
        raise RuntimeError("adata_raw.var lacks gene_symbol")
    symbols = microglia.var["gene_symbol"].astype(str).to_numpy()
    gene_columns: dict[str, list[int]] = {}
    for column, symbol in enumerate(symbols):
        if symbol in TARGET_GENES:
            gene_columns.setdefault(symbol, []).append(column)
    missing_genes = [gene for gene in TARGET_GENES if gene not in gene_columns]
    if missing_genes:
        raise RuntimeError(f"Panel genes absent from adata_raw: {missing_genes}")

    sample_ids = ["C1", "C2", "C3", "S1", "S2", "S3"]
    sample_rows: list[dict] = []
    gene_rows: list[dict] = []

    for sample in sample_ids:
        sample_original = original_obs["sample"].to_numpy() == sample
        sample_retained = retained_obs["sample"].to_numpy() == sample
        sample_counts = counts[sample_retained]
        total_umi = float(sample_counts.sum())
        n_original = int(sample_original.sum())
        n_doublets = int(
            original_obs.loc[sample_original, "predicted_doublet"].sum()
        )
        n_retained = int(sample_retained.sum())
        if n_retained != n_original - n_doublets:
            raise RuntimeError(f"Cell accounting failed for {sample}")
        sample_rows.append(
            {
                "sample": sample,
                "source_group": group_map[sample],
                "analysis_group": (
                    "combined_exposure"
                    if group_map[sample] == SOURCE_EXPOSED_LABEL
                    else "oxygen_control"
                ),
                "n_original_microglia": n_original,
                "n_scrublet_doublets_removed": n_doublets,
                "n_retained_microglia": n_retained,
                "retained_fraction": n_retained / n_original,
                "retained_pseudobulk_total_umi": total_umi,
            }
        )
        for gene in TARGET_GENES:
            columns = gene_columns[gene]
            gene_umi = float(sample_counts[:, columns].sum())
            log2cpm = float(
                np.log2((gene_umi + 0.5) / (total_umi + 1.0) * 1e6)
            )
            gene_rows.append(
                {
                    "sample": sample,
                    "source_group": group_map[sample],
                    "analysis_group": (
                        "combined_exposure"
                        if group_map[sample] == SOURCE_EXPOSED_LABEL
                        else "oxygen_control"
                    ),
                    "gene": gene,
                    "n_features_summed": len(columns),
                    "gene_umi": gene_umi,
                    "pseudobulk_total_umi": total_umi,
                    "pseudobulk_log2cpm": log2cpm,
                }
            )

    samples = pd.DataFrame(sample_rows)
    gene_values = pd.DataFrame(gene_rows)
    animal_by_gene = gene_values.pivot(
        index="sample", columns="gene", values="pseudobulk_log2cpm"
    ).loc[sample_ids, TARGET_GENES]
    panel_values = {
        panel: animal_by_gene[genes].mean(axis=1)
        for panel, genes in PANELS.items()
    }

    effects: list[dict] = []
    loo: list[dict] = []
    for gene in TARGET_GENES:
        values = animal_by_gene[gene]
        effect = effect_summary(values, group_map, gene, "gene", [gene])
        effects.append(effect)
        loo.extend(
            loo_summaries(
                values,
                group_map,
                gene,
                "gene",
                effect["mean_difference"],
            )
        )
    for panel, genes in PANELS.items():
        panel_effect = effect_summary(
            panel_values[panel], group_map, panel, "panel", genes
        )
        effects.append(panel_effect)
        loo.extend(
            loo_summaries(
                panel_values[panel],
                group_map,
                panel,
                "panel",
                panel_effect["mean_difference"],
            )
        )

    panel_animal_rows: list[dict] = []
    for panel, genes in PANELS.items():
        for sample in sample_ids:
            panel_animal_rows.append(
                {
                    "sample": sample,
                    "source_group": group_map[sample],
                    "analysis_group": (
                        "combined_exposure"
                        if group_map[sample] == SOURCE_EXPOSED_LABEL
                        else "oxygen_control"
                    ),
                    "panel": panel,
                    "n_genes": len(genes),
                    "genes": ";".join(genes),
                    "panel_pseudobulk_log2cpm_mean": panel_values[panel].loc[
                        sample
                    ],
                }
            )
    panel_animal_values = pd.DataFrame(panel_animal_rows)
    effects_df = pd.DataFrame(effects)
    loo_df = pd.DataFrame(loo)
    loo_summary = (
        loo_df.groupby(["outcome_type", "outcome"], sort=False)
        .agg(
            n_drops=("dropped_animal", "size"),
            mean_difference_min=("mean_difference", "min"),
            mean_difference_max=("mean_difference", "max"),
            cohens_d_min=("cohens_d_point_estimate", "min"),
            cohens_d_max=("cohens_d_point_estimate", "max"),
            n_sign_matches_full=("sign_matches_full", "sum"),
        )
        .reset_index()
    )

    input_audit = pd.DataFrame(
        [
            {
                "scrublet_path": str(SCRUBLET_PATH),
                "raw_path": str(RAW_PATH),
                "scrublet_n_obs": scrublet.n_obs,
                "scrublet_n_vars": scrublet.n_vars,
                "raw_n_obs": raw.n_obs,
                "raw_n_vars": raw.n_vars,
                "microglia_barcodes_missing_in_raw": len(missing_in_raw),
                "predicted_doublets_removed": removed_n,
                "retained_microglia": retained_n,
                "submitted_rare_label_total": rare_total,
                "rare_labeled_doublets_removed": rare_doublets,
                "rare_labeled_singlets_retained": retained_rare,
                "nonrare_doublets_removed": nonrare_doublets,
                "integer_count_max_residual": integer_residual,
                "minimum_raw_count": float(counts.data.min()),
                "panels": ";".join(PANELS),
                "panel_gene_sets": "|".join(
                    f"{panel}={';'.join(genes)}"
                    for panel, genes in PANELS.items()
                ),
                "duplicate_feature_counts": ";".join(
                    f"{gene}:{len(gene_columns[gene])}" for gene in TARGET_GENES
                ),
                "python_version": platform.python_version(),
                "anndata_version": version("anndata"),
                "numpy_version": np.__version__,
                "pandas_version": pd.__version__,
                "scipy_version": scipy.__version__,
            }
        ]
    )

    input_audit.to_csv(
        OUT_DIR / "doublet_negative_input_audit.csv", index=False
    )
    samples.to_csv(OUT_DIR / "doublet_negative_sample_qc.csv", index=False)
    gene_values.to_csv(
        OUT_DIR / "doublet_negative_gene_log2cpm.csv", index=False
    )
    panel_animal_values.to_csv(
        OUT_DIR / "doublet_negative_panel_animal_values.csv", index=False
    )
    effects_df.to_csv(
        OUT_DIR / "doublet_negative_gene_and_panel_effects.csv", index=False
    )
    loo_df.to_csv(
        OUT_DIR / "doublet_negative_gene_and_panel_loo.csv", index=False
    )
    loo_summary.to_csv(
        OUT_DIR / "doublet_negative_gene_and_panel_loo_summary.csv",
        index=False,
    )

    print(input_audit.to_string(index=False))
    print("\nSAMPLE ACCOUNTING")
    print(samples.to_string(index=False))
    print("\nFULL 3-vs-3 EFFECTS")
    print(
        effects_df[
            [
                "outcome_type",
                "outcome",
                "mean_difference",
                "mean_difference_ci_low",
                "mean_difference_ci_high",
                "welch_df",
                "cohens_d_point_estimate",
                "exact_perm_hits",
                "exact_perm_allocations",
                "exact_perm_p",
            ]
        ].to_string(index=False)
    )
    print("\nSYSTEMATIC LEAVE-ONE-ANIMAL-OUT SUMMARY")
    print(loo_summary.to_string(index=False))


if __name__ == "__main__":
    main()
