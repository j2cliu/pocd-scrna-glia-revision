#!/usr/bin/env python3
"""Scrublet-negative score sensitivity for the ISG7 and IEG13 panels.

This reproduces the parameter sweeps in ``audit_score_genes_sensitivity.py``
while changing only the cell-inclusion rule:

* prior audit: exclude cells whose historical ``mg_subtype2`` label is Rare;
* this audit: retain every cell with ``predicted_doublet == False``.

Expression is rebuilt from the integer-count ``adata_raw.h5ad`` object.  The
Scrublet-annotated object supplies only cell identities and metadata. The
script also runs fixed-parameter ISG7 leave-one-gene-out diagnostics and
canonical ISG7--IEG13 cell- and animal-level correlations.
"""

from __future__ import annotations

from itertools import combinations
import os
from pathlib import Path
import tempfile

# Some sandboxed conda launches cannot use Numba's package-local cache locator.
# Give Numba an explicit disposable cache directory before importing Scanpy.
NUMBA_CACHE_DIR = Path(tempfile.gettempdir()) / "pocd_scrna_numba_cache"
NUMBA_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(NUMBA_CACHE_DIR))
MPL_CACHE_DIR = Path(tempfile.gettempdir()) / "pocd_scrna_matplotlib_cache"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))
XDG_CACHE_DIR = Path(tempfile.gettempdir()) / "pocd_scrna_xdg_cache"
XDG_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(XDG_CACHE_DIR))

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import pearsonr, spearmanr, t


PROJECT = Path(os.environ["POCD_SCRNA_PROJECT_ROOT"]).expanduser().resolve()
OUT_DIR = Path(__file__).resolve().parent
CORE = ["Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3"]
IEG = [
    "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Ier2",
    "Dusp1", "Zfp36", "Hspa1a", "Hspa1b", "Socs3",
]
PANELS = {"core_isg7": CORE, "ieg13": IEG}
SURGERY_IDS = ["S1", "S2", "S3"]
CONTROL_IDS = ["C1", "C2", "C3"]
SOURCE_EXPOSED_LABEL = "Surgery"
SOURCE_CONTROL_LABEL = "Control"


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    pooled = np.sqrt(
        ((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
        / (len(a) + len(b) - 2)
    )
    return float((a.mean() - b.mean()) / pooled) if pooled > 0 else np.nan


def exact_perm_p(a: np.ndarray, b: np.ndarray) -> float:
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


def welch_mean_difference_ci(
    a: np.ndarray, b: np.ndarray, confidence: float = 0.95
) -> tuple[float, float, float]:
    """Welch--Satterthwaite CI for the exposed-minus-control mean."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    component_a = a.var(ddof=1) / len(a)
    component_b = b.var(ddof=1) / len(b)
    variance = component_a + component_b
    difference = float(a.mean() - b.mean())
    if variance <= 0:
        return difference, difference, np.nan
    denominator = (
        component_a**2 / (len(a) - 1)
        + component_b**2 / (len(b) - 1)
    )
    degrees_freedom = float(variance**2 / denominator)
    critical = float(t.ppf((1.0 + confidence) / 2.0, degrees_freedom))
    margin = critical * np.sqrt(variance)
    return difference - margin, difference + margin, degrees_freedom


def main() -> None:
    raw = ad.read_h5ad(PROJECT / "data/processed/adata_raw.h5ad")
    scrub = ad.read_h5ad(
        PROJECT / "data/processed/adata_microglia_scrublet_annotated.h5ad"
    )
    if not scrub.obs_names.is_unique:
        raise ValueError("Scrublet-annotated cell identifiers are not unique")
    if not scrub.obs_names.isin(raw.obs_names).all():
        raise ValueError("Some Scrublet-annotated cells are absent from adata_raw")
    if scrub.obs["predicted_doublet"].isna().any():
        raise ValueError("predicted_doublet contains missing values")
    if scrub.obs["predicted_doublet"].dtype != bool:
        raise ValueError("predicted_doublet must be Boolean")

    a = raw[scrub.obs_names].copy()
    raw_values = a.X.data if sparse.issparse(a.X) else np.asarray(a.X).ravel()
    integer_residual = float(np.max(np.abs(raw_values - np.rint(raw_values))))
    if integer_residual > 1e-6 or raw_values.min() < 0:
        raise ValueError(
            "adata_raw.X failed nonnegative-integer UMI validation: "
            f"max residual={integer_residual}, min={raw_values.min()}"
        )
    symbols = a.var["gene_symbol"].astype(str).to_numpy()
    keep = ~pd.Series(symbols).duplicated(keep="first").to_numpy()
    a = a[:, keep].copy()
    a.var_names = symbols[keep]
    a.var_names_make_unique()

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.obs["sample"] = scrub.obs["sample"].astype(str).to_numpy()
    a.obs["group"] = scrub.obs["group"].astype(str).to_numpy()
    a.obs["subtype"] = scrub.obs["mg_subtype2"].astype(str).to_numpy()
    a.obs["predicted_doublet"] = scrub.obs["predicted_doublet"].to_numpy()
    a = a[~a.obs["predicted_doublet"]].copy()
    if a.n_obs != 7371:
        raise ValueError(f"Expected 7,371 Scrublet-negative cells, found {a.n_obs}")

    missing = sorted(set(CORE + IEG) - set(a.var_names))
    if missing:
        raise ValueError(f"Panel genes missing after symbol mapping: {missing}")

    retained = (
        a.obs.groupby(["group", "sample"], observed=True)
        .size()
        .rename("n_retained")
        .reset_index()
    )
    retained["selection"] = "predicted_doublet_false"
    retained.to_csv(
        OUT_DIR / "score_genes_doublet_negative_retained_cells.csv", index=False
    )

    configs: list[tuple[str, int, int, int]] = []
    for seed in range(20):
        configs.append(("seed_sweep", seed, 50, 25))
    for ctrl_size in (25, 50, 100, 200):
        configs.append(("control_size_sweep", 42, ctrl_size, 25))
    for n_bins in (10, 25, 50):
        configs.append(("bin_sweep", 42, 50, n_bins))

    rows: list[dict[str, object]] = []
    for panel_name, panel_genes in PANELS.items():
        for sweep, seed, ctrl_size, n_bins in configs:
            key = f"score_{panel_name}_{sweep}_{seed}_{ctrl_size}_{n_bins}"
            sc.tl.score_genes(
                a,
                gene_list=panel_genes,
                score_name=key,
                ctrl_size=ctrl_size,
                n_bins=n_bins,
                random_state=seed,
                use_raw=False,
                ctrl_as_ref=True,
            )
            per_animal = a.obs.groupby("sample", observed=True)[key].mean()
            surgery = per_animal.loc[SURGERY_IDS].to_numpy()
            control = per_animal.loc[CONTROL_IDS].to_numpy()
            no_s3 = per_animal.loc[["S1", "S2"]].to_numpy()
            full_diff = float(surgery.mean() - control.mean())
            no_s3_diff = float(no_s3.mean() - control.mean())
            ci_low, ci_high, welch_df = welch_mean_difference_ci(
                surgery, control
            )
            rows.append(
                {
                    "selection": "predicted_doublet_false",
                    "n_cells": a.n_obs,
                    "n_genes": a.n_vars,
                    "panel": panel_name,
                    "panel_genes": ";".join(panel_genes),
                    "sweep": sweep,
                    "seed": seed,
                    "ctrl_size": ctrl_size,
                    "n_bins": n_bins,
                    "mean_difference": full_diff,
                    "mean_difference_ci_low": ci_low,
                    "mean_difference_ci_high": ci_high,
                    "welch_df": welch_df,
                    "cohens_d": cohens_d(surgery, control),
                    "exact_perm_p": exact_perm_p(surgery, control),
                    "mean_difference_without_s3": no_s3_diff,
                    "cohens_d_without_s3": cohens_d(no_s3, control),
                    "sign_preserved_without_s3": bool(
                        np.sign(no_s3_diff) == np.sign(full_diff)
                    ),
                    "surgery_values": ";".join(
                        f"{sample}:{per_animal.loc[sample]:.8g}"
                        for sample in SURGERY_IDS
                    ),
                    "control_values": ";".join(
                        f"{sample}:{per_animal.loc[sample]:.8g}"
                        for sample in CONTROL_IDS
                    ),
                }
            )
            del a.obs[key]

    result = pd.DataFrame(rows)
    result.to_csv(
        OUT_DIR / "score_genes_sensitivity_doublet_negative.csv", index=False
    )
    summary = (
        result.groupby(["selection", "panel", "sweep"], observed=True)
        .agg(
            n_runs=("cohens_d", "size"),
            n_cells=("n_cells", "first"),
            n_genes=("n_genes", "first"),
            mean_difference_min=("mean_difference", "min"),
            mean_difference_median=("mean_difference", "median"),
            mean_difference_max=("mean_difference", "max"),
            d_min=("cohens_d", "min"),
            d_median=("cohens_d", "median"),
            d_max=("cohens_d", "max"),
            exact_perm_p_min=("exact_perm_p", "min"),
            exact_perm_p_max=("exact_perm_p", "max"),
            diff_no_s3_min=("mean_difference_without_s3", "min"),
            diff_no_s3_max=("mean_difference_without_s3", "max"),
            d_no_s3_min=("cohens_d_without_s3", "min"),
            d_no_s3_max=("cohens_d_without_s3", "max"),
            n_sign_preserved_without_s3=("sign_preserved_without_s3", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(
        OUT_DIR / "score_genes_sensitivity_doublet_negative_summary.csv",
        index=False,
    )

    # Fixed-parameter leave-one-gene-out analysis for the author-selected ISG7
    # panel. This tests dependence on panel membership, not animal robustness.
    canonical = result.loc[
        (result["panel"] == "core_isg7")
        & (result["seed"] == 42)
        & (result["ctrl_size"] == 50)
        & (result["n_bins"] == 25)
    ].iloc[0]
    loo_gene_rows: list[dict[str, object]] = []
    for dropped_gene in CORE:
        retained_genes = [gene for gene in CORE if gene != dropped_gene]
        key = f"score_core_isg7_drop_{dropped_gene}"
        sc.tl.score_genes(
            a,
            gene_list=retained_genes,
            score_name=key,
            ctrl_size=50,
            n_bins=25,
            random_state=42,
            use_raw=False,
            ctrl_as_ref=True,
        )
        per_animal = a.obs.groupby("sample", observed=True)[key].mean()
        surgery = per_animal.loc[SURGERY_IDS].to_numpy()
        control = per_animal.loc[CONTROL_IDS].to_numpy()
        no_s3 = per_animal.loc[["S1", "S2"]].to_numpy()
        difference = float(surgery.mean() - control.mean())
        difference_no_s3 = float(no_s3.mean() - control.mean())
        ci_low, ci_high, welch_df = welch_mean_difference_ci(
            surgery, control
        )
        loo_gene_rows.append(
            {
                "selection": "predicted_doublet_false",
                "panel": "core_isg7",
                "dropped_gene": dropped_gene,
                "n_remaining": len(retained_genes),
                "genes_remaining": ";".join(retained_genes),
                "seed": 42,
                "ctrl_size": 50,
                "n_bins": 25,
                "mean_difference": difference,
                "mean_difference_ci_low": ci_low,
                "mean_difference_ci_high": ci_high,
                "welch_df": welch_df,
                "cohens_d": cohens_d(surgery, control),
                "exact_perm_p": exact_perm_p(surgery, control),
                "mean_difference_without_s3": difference_no_s3,
                "cohens_d_without_s3": cohens_d(no_s3, control),
                "full_isg7_mean_difference": canonical["mean_difference"],
                "full_isg7_cohens_d": canonical["cohens_d"],
                "surgery_values": ";".join(
                    f"{sample}:{per_animal.loc[sample]:.8g}"
                    for sample in SURGERY_IDS
                ),
                "control_values": ";".join(
                    f"{sample}:{per_animal.loc[sample]:.8g}"
                    for sample in CONTROL_IDS
                ),
            }
        )
        del a.obs[key]
    loo_gene = pd.DataFrame(loo_gene_rows)
    loo_gene.to_csv(
        OUT_DIR / "score_genes_doublet_negative_leave_one_gene_out.csv",
        index=False,
    )

    # Canonical-score co-variation for the nuisance assessment described in
    # Methods. Cell-level correlations are descriptive; the six animal means
    # remain the only biological-unit correlation.
    canonical_score_columns: dict[str, str] = {}
    for panel_name, panel_genes in PANELS.items():
        key = f"canonical_{panel_name}"
        sc.tl.score_genes(
            a,
            gene_list=panel_genes,
            score_name=key,
            ctrl_size=50,
            n_bins=25,
            random_state=42,
            use_raw=False,
            ctrl_as_ref=True,
        )
        canonical_score_columns[panel_name] = key
    canonical_animal = (
        a.obs.groupby("sample", observed=True)[
            list(canonical_score_columns.values())
        ]
        .mean()
        .rename(
            columns={
                canonical_score_columns["core_isg7"]: "core_isg7_score",
                canonical_score_columns["ieg13"]: "ieg13_score",
            }
        )
        .loc[CONTROL_IDS + SURGERY_IDS]
    )
    canonical_animal.insert(
        0,
        "source_group",
        [SOURCE_CONTROL_LABEL] * 3 + [SOURCE_EXPOSED_LABEL] * 3,
    )
    canonical_animal.insert(0, "sample", canonical_animal.index)
    canonical_animal.reset_index(drop=True, inplace=True)
    canonical_animal.to_csv(
        OUT_DIR / "score_genes_doublet_negative_canonical_animal_values.csv",
        index=False,
    )

    cell_isg = a.obs[canonical_score_columns["core_isg7"]].to_numpy()
    cell_ieg = a.obs[canonical_score_columns["ieg13"]].to_numpy()
    animal_isg = canonical_animal["core_isg7_score"].to_numpy()
    animal_ieg = canonical_animal["ieg13_score"].to_numpy()
    cell_pearson = pearsonr(cell_isg, cell_ieg)
    cell_spearman = spearmanr(cell_isg, cell_ieg)
    animal_pearson = pearsonr(animal_isg, animal_ieg)
    animal_spearman = spearmanr(animal_isg, animal_ieg)
    correlations = pd.DataFrame(
        [
            {
                "analysis_level": "cell_descriptive",
                "n_units": a.n_obs,
                "pearson_r": cell_pearson.statistic,
                "pearson_p": cell_pearson.pvalue,
                "spearman_rho": cell_spearman.statistic,
                "spearman_p": cell_spearman.pvalue,
            },
            {
                "analysis_level": "animal",
                "n_units": len(canonical_animal),
                "pearson_r": animal_pearson.statistic,
                "pearson_p": animal_pearson.pvalue,
                "spearman_rho": animal_spearman.statistic,
                "spearman_p": animal_spearman.pvalue,
            },
        ]
    )
    correlations.insert(0, "selection", "predicted_doublet_false")
    correlations.insert(1, "seed", 42)
    correlations.insert(2, "ctrl_size", 50)
    correlations.insert(3, "n_bins", 25)
    correlations.to_csv(
        OUT_DIR / "score_genes_doublet_negative_isg_ieg_correlation.csv",
        index=False,
    )
    for key in canonical_score_columns.values():
        del a.obs[key]

    # Machine-readable comparison with the previously generated Rare-excluded
    # sensitivity audit. That file is reproducible from
    # audit_score_genes_sensitivity.py in the same directory.
    prior_path = OUT_DIR / "score_genes_sensitivity.csv"
    if not prior_path.exists():
        raise FileNotFoundError(
            "Run audit_score_genes_sensitivity.py before this script so the "
            "selection-rule comparison can be generated"
        )
    prior = pd.read_csv(prior_path)
    prior = prior.loc[prior["panel"] == "core_isg7"].copy()
    current_core = result.loc[result["panel"] == "core_isg7"].copy()
    keys = ["panel", "sweep", "seed", "ctrl_size", "n_bins"]
    comparison = prior.merge(
        current_core,
        on=keys,
        how="inner",
        validate="one_to_one",
        suffixes=("_rare_excluded", "_doublet_negative"),
    )
    if len(comparison) != len(current_core):
        raise ValueError(
            f"Expected {len(current_core)} matched settings, "
            f"found {len(comparison)}"
        )
    comparison.insert(0, "old_selection", "submitted_rare_label_excluded")
    comparison.insert(1, "old_n_cells", 7367)
    comparison.insert(2, "new_selection", "predicted_doublet_false")
    comparison.insert(3, "new_n_cells", a.n_obs)
    comparison["change_mean_difference"] = (
        comparison["mean_difference_doublet_negative"]
        - comparison["mean_difference_rare_excluded"]
    )
    comparison["change_cohens_d"] = (
        comparison["cohens_d_doublet_negative"]
        - comparison["cohens_d_rare_excluded"]
    )
    comparison["change_mean_difference_without_s3"] = (
        comparison["mean_difference_without_s3_doublet_negative"]
        - comparison["mean_difference_without_s3_rare_excluded"]
    )
    comparison["change_cohens_d_without_s3"] = (
        comparison["cohens_d_without_s3_doublet_negative"]
        - comparison["cohens_d_without_s3_rare_excluded"]
    )
    comparison["full_sign_changed"] = (
        np.sign(comparison["mean_difference_doublet_negative"])
        != np.sign(comparison["mean_difference_rare_excluded"])
    )
    comparison["without_s3_sign_changed"] = (
        np.sign(comparison["mean_difference_without_s3_doublet_negative"])
        != np.sign(comparison["mean_difference_without_s3_rare_excluded"])
    )
    comparison.to_csv(
        OUT_DIR / "score_genes_selection_rule_comparison.csv", index=False
    )
    comparison_summary = (
        comparison.groupby("sweep", observed=True)
        .agg(
            n_matched_settings=("seed", "size"),
            max_abs_change_mean_difference=(
                "change_mean_difference",
                lambda x: float(np.abs(x).max()),
            ),
            max_abs_change_cohens_d=(
                "change_cohens_d",
                lambda x: float(np.abs(x).max()),
            ),
            max_abs_change_mean_difference_without_s3=(
                "change_mean_difference_without_s3",
                lambda x: float(np.abs(x).max()),
            ),
            max_abs_change_cohens_d_without_s3=(
                "change_cohens_d_without_s3",
                lambda x: float(np.abs(x).max()),
            ),
            n_full_sign_changed=("full_sign_changed", "sum"),
            n_without_s3_sign_changed=("without_s3_sign_changed", "sum"),
        )
        .reset_index()
    )
    comparison_summary.to_csv(
        OUT_DIR / "score_genes_selection_rule_comparison_summary.csv",
        index=False,
    )
    print(retained.to_string(index=False))
    print()
    print(summary.to_string(index=False))
    print()
    print(comparison_summary.to_string(index=False))
    print()
    print(loo_gene.to_string(index=False))
    print()
    print(correlations.to_string(index=False))


if __name__ == "__main__":
    main()
