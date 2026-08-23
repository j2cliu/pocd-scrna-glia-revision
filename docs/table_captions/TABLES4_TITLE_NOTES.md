# Table S4 title and notes

## Title

**Table S4. Complete Hallmark enrichment, background-adjusted scoring sensitivity, animal-influence summaries, and fixed-score component analyses for GSE267933 microglia.**

## General note

The primary molecular analysis used 7,371 fixed reference microglia after exclusion of 90 Scrublet-predicted doublets. C1–C3 denote oxygen-control animals; S1–S3 denote animals receiving combined anesthesia-plus-laparotomy exposure. The animal/library is the biological unit. Direct-expression differences in log2 counts per million, background-adjusted `score_genes` differences, and normalized enrichment scores (NES) are distinct estimands and are not treated as mutually validating measurements.

Animal omissions assess individual-animal influence. Scoring-parameter sweeps assess estimator sensitivity, and transcript omissions assess dependence on panel definition. These analyses use the same six-animal cohort, provide no additional biological replication or orthogonal validation, and did not remove any animal from the full-cohort estimate.

## Part A — Complete mouse MSigDB Hallmark analysis

Part A reports all 50 mouse MSigDB 2026.1.Mm Hallmark sets for the full animal-level pseudobulk fit and six systematic one-animal-omission refits. Genes were ranked by the signed DESeq2 Wald statistic. Positive NES values denote ranking toward combined exposure and negative values denote ranking toward oxygen control; NES does not measure protein abundance or pathway activity. Benjamini–Hochberg adjustment was performed across all 50 sets separately within each fit. Animal-omission rows are influence summaries rather than independent hypothesis tests.

## Part B — Selected-panel `score_genes` parameter sensitivity

Part B reports all one-parameter-at-a-time scoring runs for the seven selected interferon-responsive transcripts. Seeds 0–19 were evaluated with `ctrl_size = 50` and `n_bins = 25`; `ctrl_size` values of 25, 50, 100, and 200 were evaluated with seed 42 and `n_bins = 25`; and `n_bins` values of 10, 25, and 50 were evaluated with seed 42 and `ctrl_size = 50`. The seed-42, `ctrl_size = 50`, `n_bins = 25` configuration appears in both the control-size and bin-number sweeps by design and is not an additional analysis. `n_expression_features` describes the symbol-deduplicated expression universe; `panel_size` and `panel_genes` describe panel membership.

Each configuration reports the full-cohort effect and the contrast after omitting S3. The latter assesses whether the sign change identified by the systematic fixed-score animal-omission analysis persisted across scoring configurations; it is not an independent experiment.

## Part C — Selected-panel parameter-sweep summaries

Part C summarizes Part B by parameter family. It reports the number of runs and ranges or medians of full-cohort mean differences, Cohen’s *d*, exact permutation *P* values, and omit-S3 contrasts. `n_sign_preserved_omit_s3` counts runs in which the full-cohort direction was retained after S3 was omitted. These ranges describe scoring choices within the same six animals and are not confidence intervals, multiplicity-adjusted inference, or biological replication.

## Part D — Leave-one-transcript-out panel-definition sensitivity

Part D reports seven fixed-configuration runs in which one transcript was omitted from the seven-transcript panel. `omitted_transcript` identifies the excluded transcript and is distinct from the animal-omission fields in Parts A, E, F, H, and K. Each six-transcript analysis reports the full-cohort score and the contrast after omitting S3. These rows assess dependence on panel membership rather than biological-replicate robustness.

## Part E — Selected-panel predicted-doublet-inclusion sensitivity

Part E compares direct selected-panel mean differences from the primary 7,371-cell Scrublet-negative population with estimates from all 7,461 fixed reference microglia, including the 90 Scrublet-predicted doublets. It reports the full cohort and all six animal omissions. Fields prefixed `all_7461_microglia_` refer to all fixed reference microglia rather than all brain-cell types. The same six animals contribute under both cell-inclusion rules, so similarity is a within-cohort inclusion-sensitivity result rather than biological replication or validation.

## Part F — Focused-Hallmark predicted-doublet-inclusion sensitivity

Part F reports paired NES values for the four focused Hallmark sets under the 7,371-cell and 7,461-cell inclusion rules for the full cohort and six animal omissions. `scrublet_negative_adjusted_p_value` and `all_7461_microglia_adjusted_p_value` were each calculated across all 50 Hallmarks within the corresponding fit, although only the four focused sets are included in this part. Agreement in NES direction or false-discovery-rate classification does not establish pathway activity, biological replication, or robustness across independent cohorts.

## Part G — Fixed selected-panel score animal values

Part G reports the six animal-level values underlying Figure 3A. Values are animal means from `scanpy.tl.score_genes` with seed 42, `ctrl_size = 50`, and `n_bins = 25`. The score represents selected-transcript expression relative to algorithmically matched control-gene expression on the log-normalized matrix. It is distinct from direct expression in log2 counts per million and does not measure pathway activity.

## Part H — Fixed selected-panel score effects and animal omissions

Part H reports the full-cohort effect and six systematic one-animal-omission contrasts calculated from the fixed animal values in Part G. The full three-versus-three row includes the combined-exposure-minus-oxygen-control mean difference, Welch–Satterthwaite confidence interval, Cohen’s *d*, and descriptive exact label-permutation *P* value. Exact permutation results are not assigned to unequal-group omission rows. Omission rows are influence summaries calculated from fixed scores rather than independent hypothesis tests or score refits. S3 remains in the full-cohort estimate.

## Part I — Fixed selected-panel matched-control gene manifest

Part I reports the exact 249 algorithm-selected matched-control genes used for the fixed seven-transcript score, with matrix indices, expression bins, mean log1p-normalized expression, scoring configuration, and control-set SHA-256. The set was selected once from the 7,371-cell, 27,933-feature symbol-deduplicated matrix with seed 42, `ctrl_size = 50`, `n_bins = 25`, `use_raw = False`, and `ctrl_as_ref = True`; it was held fixed for the component and animal-omission calculations in Parts J–K. These genes are algorithmic expression-matched references rather than a prespecified negative-control panel or biological background program.

## Part J — Fixed selected-panel score components by animal

Part J reports the selected-gene mean, fixed matched-control-gene mean, and their difference for each animal. The difference reconstructs the stored Scanpy score within numerical tolerance. Components were calculated at the cell level and averaged within animal; the cell count is reported for each animal. The decomposition is on the log-normalized score scale and is not interchangeable with the direct animal-level expression measure in Figure 2.

## Part K — Fixed selected-panel component contrasts and animal omissions

Part K reports combined-exposure and oxygen-control means and their differences separately for the selected, fixed matched-control, and final score components in the full cohort and six animal-omission scenarios. In each row, the selected-component difference minus the matched-control-component difference reconstructs the score difference within numerical tolerance. Omission rows reuse the fixed component values and are influence summaries rather than score refits or independent tests. The arithmetic decomposition does not assign biological meaning to the matched-control genes and does not negate the separate direct-expression result.
