#!/usr/bin/env python3
"""Freeze and decompose the fixed Scanpy ``score_genes`` configuration.

This bounded audit is specific to the author-selected seven-transcript panel
and the primary 7,371-cell Scrublet-negative analysis set.  It reconstructs
the exact matched-control genes used by Scanpy 1.11.5 at seed 42,
``ctrl_size=50``, ``n_bins=25``, ``use_raw=False``, and
``ctrl_as_ref=True``.  It then separates the per-cell score into:

    selected-gene component - matched-control component = final score

Components are averaged within animal before condition contrasts are formed.
The matched-control component is an algorithmic reference, not a biologically
neutral background program.  Full-cohort and animal-omission rows are
descriptive decompositions of the same six animals, not additional tests or
biological replication.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import inspect
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import sys
import tempfile

# Scanpy/Numba need writable cache locations in the sandboxed environment.
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
from scanpy.tools._score_genes import (  # type: ignore[attr-defined]
    _check_score_genes_args,
    _nan_means,
    _score_genes_bins,
)
from scipy import sparse


DEFAULT_PROJECT = Path(os.environ["POCD_SCRNA_PROJECT_ROOT"]).expanduser().resolve()
CORE = ["Irf7", "Ifitm3", "Isg15", "Mx1", "Ifit1", "Ifit2", "Ifit3"]
SAMPLES = ["C1", "C2", "C3", "S1", "S2", "S3"]
SOURCE_GROUP = {
    "C1": "Control",
    "C2": "Control",
    "C3": "Control",
    "S1": "Surgery",
    "S2": "Surgery",
    "S3": "Surgery",
}
DISPLAY_GROUP = {
    "Control": "Oxygen control",
    "Surgery": "Combined exposure",
}
SEED = 42
CTRL_SIZE = 50
N_BINS = 25
NORMALIZATION_TARGET = 10_000.0
EXPECTED_N_CELLS = 7_371
EXPECTED_N_FEATURES = 27_933
EXPECTED_N_CONTROL_GENES = 249
EXPECTED_SCANPY_VERSION = "1.11.5"
EXPECTED_ANNDATA_VERSION = "0.12.10"
EXPECTED_RAW_SHA256 = (
    "72bc45907d1d939a76f609c11087897f8bd4636f48205e67bbbcc28f6d981b80"
)
EXPECTED_SCRUB_SHA256 = (
    "fc6c595b1589caafd7d84ab71bd688182ac926694479ba734c5b071c9faecaeb"
)
EXPECTED_CANONICAL_SCORE_SHA256 = (
    "b814dc3905b26f25326b7b6267500f615b9b2def9d5871d5b718e0126dd42359"
)
EXPECTED_FEATURE_UNIVERSE_SHA256 = (
    "d45a6462afe691448a866b2787a567fdf3fd14e8b533e6cc816dadfa04c60e39"
)
EXPECTED_PANEL_SHA256 = (
    "618740475974430a371ca4926e660cfd372d25b3b8f4a2129bfb2764300626ea"
)
EXPECTED_CONTROL_SET_SHA256 = (
    "38af968fe25becbf6dd7fa95d05fb508819a9a51ac90f1cd0c3a5d35fedb6a1f"
)
EXPECTED_SCANPY_SCORE_GENES_SOURCE_SHA256 = (
    "1974378ba371cbcd317bda1afb50a3c4ab91db2289cd28c0efa0d4b6c22bab1f"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=DEFAULT_PROJECT,
        help="POCD scRNA-seq project root containing data/processed.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for machine-readable audit outputs.",
    )
    parser.add_argument(
        "--canonical-score-path",
        type=Path,
        default=(
            Path(__file__).resolve().parent
            / "score_genes_doublet_negative_canonical_animal_values.csv"
        ),
        help=(
            "Frozen canonical animal-score CSV used only as a cross-check; "
            "defaults to the audit-script directory."
        ),
    )
    return parser.parse_args()


def assert_true(condition: bool, message: str) -> None:
    if not bool(condition):
        raise AssertionError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: list[str]) -> str:
    payload = "".join(f"{value}\n" for value in values).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def contrast_rows(animal: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for scenario_order, omitted_animal in enumerate(["(none)", *SAMPLES], 1):
        local = animal.copy()
        if omitted_animal != "(none)":
            local = local.loc[local["sample"] != omitted_animal].copy()
        exposed = local.loc[local["source_group"] == "Surgery"]
        control = local.loc[local["source_group"] == "Control"]

        row: dict[str, object] = {
            "scenario_order": scenario_order,
            "scenario": (
                "Full cohort"
                if omitted_animal == "(none)"
                else f"Omit {omitted_animal}"
            ),
            "omitted_animal": omitted_animal,
            "omitted_group": (
                pd.NA
                if omitted_animal == "(none)"
                else DISPLAY_GROUP[SOURCE_GROUP[omitted_animal]]
            ),
            "n_combined_exposure": len(exposed),
            "n_oxygen_control": len(control),
        }
        for component in (
            "selected_component",
            "matched_control_component",
            "score",
        ):
            exposed_mean = float(exposed[component].mean())
            control_mean = float(control[component].mean())
            row[f"{component}_combined_exposure_mean"] = exposed_mean
            row[f"{component}_oxygen_control_mean"] = control_mean
            row[f"{component}_mean_difference"] = exposed_mean - control_mean

        row["contrast_identity_residual"] = float(
            row["score_mean_difference"]
            - (
                row["selected_component_mean_difference"]
                - row["matched_control_component_mean_difference"]
            )
        )
        row["inference_status"] = (
            "Fixed-score descriptive component decomposition; same six animals"
            if omitted_animal == "(none)"
            else (
                "Fixed-score descriptive component decomposition after one "
                "animal omission; influence diagnostic, not replication"
            )
        )
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    project_root = args.project_root.expanduser().resolve(strict=True)
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    script_path = Path(__file__).resolve()
    raw_path = project_root / "data/processed/adata_raw.h5ad"
    scrub_path = (
        project_root
        / "data/processed/adata_microglia_scrublet_annotated.h5ad"
    )
    expected_score_path = args.canonical_score_path.expanduser().resolve(strict=True)
    score_genes_source_path = Path(
        inspect.getsourcefile(_score_genes_bins) or ""
    ).resolve(strict=True)
    for path in (raw_path, scrub_path, expected_score_path):
        assert_true(path.exists(), f"Required input is missing: {path}")

    frozen_input_hashes = {
        raw_path: EXPECTED_RAW_SHA256,
        scrub_path: EXPECTED_SCRUB_SHA256,
        expected_score_path: EXPECTED_CANONICAL_SCORE_SHA256,
        score_genes_source_path: EXPECTED_SCANPY_SCORE_GENES_SOURCE_SHA256,
    }
    for path, expected_hash in frozen_input_hashes.items():
        observed_hash = sha256_file(path)
        assert_true(
            observed_hash == expected_hash,
            f"Frozen input/source SHA-256 changed for {path}: {observed_hash}",
        )

    assert_true(
        version("scanpy") == EXPECTED_SCANPY_VERSION,
        f"Expected Scanpy {EXPECTED_SCANPY_VERSION}, found {version('scanpy')}",
    )
    assert_true(
        version("anndata") == EXPECTED_ANNDATA_VERSION,
        f"Expected anndata {EXPECTED_ANNDATA_VERSION}, found {version('anndata')}",
    )

    raw = ad.read_h5ad(raw_path)
    scrub = ad.read_h5ad(scrub_path)
    assert_true(scrub.obs_names.is_unique, "Scrublet cell identifiers are not unique")
    assert_true(
        scrub.obs_names.isin(raw.obs_names).all(),
        "Some Scrublet-annotated cells are absent from adata_raw",
    )
    assert_true(
        scrub.obs["predicted_doublet"].dtype == bool,
        "predicted_doublet must be Boolean",
    )
    assert_true(
        not scrub.obs["predicted_doublet"].isna().any(),
        "predicted_doublet contains missing values",
    )

    a = raw[scrub.obs_names].copy()
    raw_values = a.X.data if sparse.issparse(a.X) else np.asarray(a.X).ravel()
    integer_residual = float(np.max(np.abs(raw_values - np.rint(raw_values))))
    assert_true(
        integer_residual <= 1e-6 and float(raw_values.min()) >= 0,
        "adata_raw.X failed nonnegative-integer UMI validation",
    )

    symbols = a.var["gene_symbol"].astype(str).to_numpy()
    keep = ~pd.Series(symbols).duplicated(keep="first").to_numpy()
    a = a[:, keep].copy()
    a.var_names = symbols[keep]
    a.var_names_make_unique()

    sc.pp.normalize_total(a, target_sum=NORMALIZATION_TARGET)
    sc.pp.log1p(a)
    a.obs["sample"] = scrub.obs["sample"].astype(str).to_numpy()
    a.obs["source_group"] = scrub.obs["group"].astype(str).to_numpy()
    a.obs["predicted_doublet"] = scrub.obs["predicted_doublet"].to_numpy()
    a = a[~a.obs["predicted_doublet"]].copy()

    assert_true(a.n_obs == EXPECTED_N_CELLS, f"Expected {EXPECTED_N_CELLS} cells")
    assert_true(
        a.n_vars == EXPECTED_N_FEATURES,
        f"Expected {EXPECTED_N_FEATURES} symbol-deduplicated features",
    )
    assert_true(
        set(a.obs["sample"].astype(str)) == set(SAMPLES),
        "Unexpected sample identities after filtering",
    )
    observed_group = (
        a.obs[["sample", "source_group"]]
        .drop_duplicates()
        .set_index("sample")["source_group"]
        .astype(str)
        .to_dict()
    )
    assert_true(observed_group == SOURCE_GROUP, "Sample/group mapping changed")
    missing_panel = sorted(set(CORE) - set(a.var_names))
    assert_true(not missing_panel, f"Panel genes are missing: {missing_panel}")
    feature_universe_sha256 = sha256_lines(a.var_names.astype(str).tolist())
    panel_sha256 = sha256_lines(CORE)
    assert_true(
        feature_universe_sha256 == EXPECTED_FEATURE_UNIVERSE_SHA256,
        "Symbol-deduplicated feature-universe SHA-256 changed",
    )
    assert_true(panel_sha256 == EXPECTED_PANEL_SHA256, "Panel SHA-256 changed")

    # Reconstruct the exact private Scanpy control-selection path.  Version
    # assertions above intentionally fail if this private implementation drifts.
    gene_list, gene_pool, get_subset = _check_score_genes_args(
        a,
        CORE,
        None,
        use_raw=False,
        layer=None,
    )
    np.random.seed(SEED)
    control_genes = pd.Index([], dtype="string")
    for genes in _score_genes_bins(
        gene_list,
        gene_pool,
        ctrl_as_ref=True,
        ctrl_size=CTRL_SIZE,
        n_bins=N_BINS,
        get_subset=get_subset,
    ):
        control_genes = control_genes.union(genes)
    assert_true(
        len(control_genes) == EXPECTED_N_CONTROL_GENES,
        f"Expected {EXPECTED_N_CONTROL_GENES} control genes, found {len(control_genes)}",
    )
    control_set_sha256 = sha256_lines(control_genes.astype(str).tolist())
    assert_true(
        control_set_sha256 == EXPECTED_CONTROL_SET_SHA256,
        "Fixed matched-control set SHA-256 changed",
    )

    selected_component = np.asarray(
        _nan_means(get_subset(gene_list), axis=1, dtype="float64")
    ).ravel()
    matched_control_component = np.asarray(
        _nan_means(get_subset(control_genes), axis=1, dtype="float64")
    ).ravel()
    reconstructed_score = selected_component - matched_control_component

    score_name = "fixed_author_selected_isg7"
    sc.tl.score_genes(
        a,
        gene_list=CORE,
        score_name=score_name,
        ctrl_size=CTRL_SIZE,
        gene_pool=None,
        n_bins=N_BINS,
        random_state=SEED,
        use_raw=False,
        ctrl_as_ref=True,
    )
    scanpy_score = a.obs[score_name].to_numpy(dtype=float)
    per_cell_max_abs_residual = float(
        np.max(np.abs(scanpy_score - reconstructed_score))
    )
    assert_true(
        per_cell_max_abs_residual == 0.0,
        "Reconstructed selected-minus-control score does not exactly match Scanpy",
    )

    # Recreate Scanpy's expression bins for an auditable control-gene manifest.
    obs_avg = pd.Series(
        np.asarray(_nan_means(get_subset(gene_pool), axis=0)).ravel(),
        index=gene_pool,
        dtype=float,
    )
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (N_BINS - 1)))
    obs_cut = obs_avg.rank(method="min") // n_items
    relevant_bins = sorted(int(value) for value in np.unique(obs_cut.loc[gene_list]))

    control_manifest = pd.DataFrame(
        {
            "control_gene_order": np.arange(1, len(control_genes) + 1),
            "control_gene": control_genes.astype(str),
        }
    )
    control_manifest["matrix_feature_index_1based"] = [
        int(a.var_names.get_loc(gene)) + 1
        for gene in control_manifest["control_gene"]
    ]
    control_manifest["expression_bin"] = [
        int(obs_cut.loc[gene]) for gene in control_manifest["control_gene"]
    ]
    control_manifest["mean_log1p_normalized_expression"] = [
        float(obs_avg.loc[gene]) for gene in control_manifest["control_gene"]
    ]
    control_manifest["fixed_configuration"] = (
        "Scanpy 1.11.5; seed=42; ctrl_size=50; n_bins=25; "
        "use_raw=False; ctrl_as_ref=True; full gene pool"
    )
    control_manifest["interpretation_limit"] = (
        "Algorithmic matched-control gene; not a validated neutral biological background"
    )

    cell_components = pd.DataFrame(
        {
            "sample": a.obs["sample"].astype(str).to_numpy(),
            "selected_component": selected_component,
            "matched_control_component": matched_control_component,
            "score": scanpy_score,
        },
        index=a.obs_names,
    )
    animal_components = (
        cell_components.groupby("sample", observed=True)
        .agg(
            n_cells=("score", "size"),
            selected_component=("selected_component", "mean"),
            matched_control_component=("matched_control_component", "mean"),
            score=("score", "mean"),
        )
        .reindex(SAMPLES)
        .reset_index()
    )
    animal_components.insert(1, "sample_order", np.arange(1, len(SAMPLES) + 1))
    animal_components.insert(
        2,
        "source_group",
        animal_components["sample"].map(SOURCE_GROUP),
    )
    animal_components.insert(
        3,
        "display_group",
        animal_components["source_group"].map(DISPLAY_GROUP),
    )
    animal_components["component_identity_residual"] = (
        animal_components["score"]
        - (
            animal_components["selected_component"]
            - animal_components["matched_control_component"]
        )
    )
    assert_true(
        int(animal_components["n_cells"].sum()) == EXPECTED_N_CELLS,
        "Animal component cell counts do not sum to 7,371",
    )
    assert_true(
        float(animal_components["component_identity_residual"].abs().max())
        < 1e-15,
        "Animal component identity residual exceeds tolerance",
    )

    expected_scores = pd.read_csv(expected_score_path).set_index("sample")
    expected_scores = expected_scores.reindex(SAMPLES)
    expected_score_max_abs_residual = float(
        np.max(
            np.abs(
                animal_components.set_index("sample")["score"].to_numpy()
                - expected_scores["core_isg7_score"].to_numpy(dtype=float)
            )
        )
    )
    assert_true(
        expected_score_max_abs_residual < 1e-15,
        "Animal score components differ from the frozen canonical animal scores",
    )

    effects = contrast_rows(animal_components)
    assert_true(
        float(effects["contrast_identity_residual"].abs().max()) < 1e-15,
        "Contrast component identity residual exceeds tolerance",
    )

    full = effects.loc[effects["omitted_animal"] == "(none)"].iloc[0]
    omit_s3 = effects.loc[effects["omitted_animal"] == "S3"].iloc[0]
    expected_anchors = {
        "full_selected": 0.0574906460808467,
        "full_control": 0.010995934902713278,
        "full_score": 0.04649471117813342,
        "omit_s3_selected": 0.005293989180154848,
        "omit_s3_control": 0.011383651690556565,
        "omit_s3_score": -0.00608966251040171,
    }
    observed_anchors = {
        "full_selected": full["selected_component_mean_difference"],
        "full_control": full["matched_control_component_mean_difference"],
        "full_score": full["score_mean_difference"],
        "omit_s3_selected": omit_s3["selected_component_mean_difference"],
        "omit_s3_control": omit_s3["matched_control_component_mean_difference"],
        "omit_s3_score": omit_s3["score_mean_difference"],
    }
    assert_true(
        max(
            abs(float(observed_anchors[key]) - value)
            for key, value in expected_anchors.items()
        )
        < 1e-12,
        "Component contrasts differ from independently checked anchors",
    )

    control_path = out_dir / "score_genes_fixed_component_control_genes.csv"
    animal_path = out_dir / "score_genes_fixed_component_animal_values.csv"
    effects_path = out_dir / "score_genes_fixed_component_effects_full_and_loo.csv"
    config_path = out_dir / "score_genes_fixed_component_configuration.csv"
    input_manifest_path = out_dir / "score_genes_fixed_component_input_manifest.csv"
    execution_path = out_dir / "score_genes_fixed_component_execution.json"
    session_path = out_dir / "score_genes_fixed_component_session.txt"
    output_manifest_path = out_dir / "score_genes_fixed_component_output_manifest.csv"

    control_manifest.to_csv(control_path, index=False)
    animal_components.to_csv(animal_path, index=False)
    effects.to_csv(effects_path, index=False)

    configuration = pd.DataFrame(
        [
            {
                "analysis_set": "Primary: 7,371 Scrublet-negative microglia",
                "biological_unit": "Animal/library",
                "panel": "author_selected_isg7",
                "panel_genes": ";".join(CORE),
                "panel_size": len(CORE),
                "selection_status": "Author-selected, post hoc, exploratory",
                "normalization": "normalize_total target_sum=10000, then log1p",
                "symbol_deduplication": "Retain first Ensembl feature per gene symbol",
                "n_expression_features": a.n_vars,
                "gene_pool": "Full symbol-deduplicated matrix",
                "seed": SEED,
                "ctrl_size": CTRL_SIZE,
                "n_bins": N_BINS,
                "use_raw": False,
                "ctrl_as_ref": True,
                "n_control_genes": len(control_genes),
                "relevant_expression_bins": ";".join(map(str, relevant_bins)),
                "control_set_sha256": control_set_sha256,
                "feature_universe_sha256": feature_universe_sha256,
                "panel_sha256": panel_sha256,
                "scanpy_version": version("scanpy"),
                "scanpy_score_genes_source_sha256": sha256_file(
                    score_genes_source_path
                ),
                "anndata_version": version("anndata"),
                "numpy_distribution_version": version("numpy"),
                "numpy_imported_runtime_version": np.__version__,
                "numpy_version_metadata_runtime_discordant": (
                    version("numpy") != np.__version__
                ),
                "pandas_version": version("pandas"),
                "scipy_version": version("scipy"),
                "per_cell_score_identity_max_abs_residual": per_cell_max_abs_residual,
                "canonical_animal_score_max_abs_residual": expected_score_max_abs_residual,
                "interpretation_limit": (
                    "The matched-control component is an algorithmic reference; "
                    "the decomposition is not a biological pathway or mechanism"
                ),
            }
        ]
    )
    configuration.to_csv(config_path, index=False)

    input_manifest = pd.DataFrame(
        [
            {
                "input_key": "raw_integer_counts",
                "path": str(raw_path),
                "bytes": raw_path.stat().st_size,
                "sha256": sha256_file(raw_path),
                "role": "Integer UMI matrix and gene symbols",
            },
            {
                "input_key": "scrublet_annotations",
                "path": str(scrub_path),
                "bytes": scrub_path.stat().st_size,
                "sha256": sha256_file(scrub_path),
                "role": "Frozen cell identities, samples, groups, and doublet calls",
            },
            {
                "input_key": "canonical_animal_score_crosscheck",
                "path": str(expected_score_path),
                "bytes": expected_score_path.stat().st_size,
                "sha256": sha256_file(expected_score_path),
                "role": "Cross-check only; not used to derive components",
            },
            {
                "input_key": "scanpy_score_genes_implementation",
                "path": str(score_genes_source_path),
                "bytes": score_genes_source_path.stat().st_size,
                "sha256": sha256_file(score_genes_source_path),
                "role": (
                    "Frozen Scanpy control-selection and score implementation; "
                    "private helper semantics are hash-locked"
                ),
            },
        ]
    )
    input_manifest.to_csv(input_manifest_path, index=False)

    session_lines = [
        f"python={platform.python_version()}",
        f"platform={platform.platform()}",
        f"scanpy={version('scanpy')}",
        f"anndata={version('anndata')}",
        f"numpy_distribution_metadata={version('numpy')}",
        f"numpy_imported_runtime={np.__version__}",
        "numpy_metadata_runtime_discordant="
        f"{version('numpy') != np.__version__}",
        f"pandas={version('pandas')}",
        f"scipy={version('scipy')}",
    ]
    session_path.write_text("\n".join(session_lines) + "\n", encoding="utf-8")

    execution = {
        "audit": "fixed score_genes component decomposition",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(script_path),
        "script_sha256": sha256_file(script_path),
        "analysis_set_cells": a.n_obs,
        "symbol_deduplicated_features": a.n_vars,
        "selected_genes": CORE,
        "control_gene_count": len(control_genes),
        "control_set_sha256": control_set_sha256,
        "feature_universe_sha256": feature_universe_sha256,
        "panel_sha256": panel_sha256,
        "scanpy_score_genes_source": str(score_genes_source_path),
        "scanpy_score_genes_source_sha256": sha256_file(
            score_genes_source_path
        ),
        "configuration": {
            "scanpy": version("scanpy"),
            "numpy_distribution_metadata": version("numpy"),
            "numpy_imported_runtime": np.__version__,
            "numpy_metadata_runtime_discordant": version("numpy") != np.__version__,
            "seed": SEED,
            "ctrl_size": CTRL_SIZE,
            "n_bins": N_BINS,
            "use_raw": False,
            "ctrl_as_ref": True,
            "gene_pool": "full symbol-deduplicated matrix",
        },
        "identity_checks": {
            "per_cell_max_abs_residual": per_cell_max_abs_residual,
            "canonical_animal_score_max_abs_residual": expected_score_max_abs_residual,
            "animal_component_max_abs_residual": float(
                animal_components["component_identity_residual"].abs().max()
            ),
            "contrast_component_max_abs_residual": float(
                effects["contrast_identity_residual"].abs().max()
            ),
        },
        "claim_ceiling": (
            "Fixed-score mathematical decomposition only; matched controls are "
            "not a validated neutral biological background, and omission rows "
            "are not independent replication"
        ),
        "environment_note": (
            "The installed NumPy distribution metadata and imported runtime "
            "version are recorded separately because this environment reports "
            "different values; reproducibility is therefore anchored to the "
            "frozen control-gene manifest/hash and exact score-identity checks."
        ),
    }
    execution_path.write_text(
        json.dumps(execution, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    tracked_outputs = [
        (control_path, "scientific_data", True),
        (animal_path, "scientific_data", True),
        (effects_path, "scientific_data", True),
        (config_path, "scientific_configuration", True),
        (input_manifest_path, "provenance", True),
        (session_path, "environment_provenance", True),
        (execution_path, "timestamped_run_metadata", False),
    ]
    output_manifest = pd.DataFrame(
        [
            {
                "output_file": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "output_class": output_class,
                "expected_byte_stable_on_identical_inputs": expected_stable,
            }
            for path, output_class, expected_stable in tracked_outputs
        ]
    )
    output_manifest.to_csv(output_manifest_path, index=False)

    print(
        "Fixed score_genes component audit complete: "
        f"{a.n_obs:,} cells, {a.n_vars:,} features, "
        f"{len(control_genes)} matched-control genes."
    )
    print(
        "Full contrasts: selected="
        f"{float(full['selected_component_mean_difference']):+.9f}; "
        "matched-control="
        f"{float(full['matched_control_component_mean_difference']):+.9f}; "
        f"score={float(full['score_mean_difference']):+.9f}."
    )
    print(
        "Omit-S3 contrasts: selected="
        f"{float(omit_s3['selected_component_mean_difference']):+.9f}; "
        "matched-control="
        f"{float(omit_s3['matched_control_component_mean_difference']):+.9f}; "
        f"score={float(omit_s3['score_mean_difference']):+.9f}."
    )


if __name__ == "__main__":
    main()
