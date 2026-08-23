# GLIA major revision: frozen Results–display map v2.1

Date: 2026-08-22  
Status: **frozen for display construction and subsequent Results rewriting**  
Scope: display hierarchy only. This file does not change an analysis population,
estimand, statistical model, numerical result, or claim boundary. The earlier
`RESULTS_FIGURE_MAPPING.md` v1.7 remains a historical build record; where its
display hierarchy conflicts with this file, this v2.0 map governs.

## 1. Evidence hierarchy

The main displays follow the analysis order in the rewritten Main Methods:

1. study design, analysis populations, and partition limitations;
2. animal-level transcriptome-wide differential expression;
3. secondary descriptive expression of seven selected interferon-responsive
   transcripts;
4. complementary Hallmark, background-adjusted, animal-influence, and
   processing-sensitivity analyses;
5. external cross-experiment context.

The seven-transcript panel is a secondary descriptive summary. It is not the
sole or first display of the animal-level molecular results. Main-figure space
must therefore show the complete 13,926-feature analysis before displaying the
selected transcripts.

## 2. Frozen Results sections and display assignments

| Results section | Primary display | Supporting displays | Required interpretive boundary |
|---|---|---|---|
| 1. The reference partition was computationally sensitive and did not show consistent animal-level expansion | Figure 1 | Tables 1–2; Figures S1–S2; Tables S1–S2 | The reference partition is computationally sensitive; composition results are conditional and do not establish reproducible state expansion. |
| 2. Transcriptome-wide differences were mixed, while selected interferon-responsive transcripts varied among animals | Figure 2 | Table S3 | Transcriptome-wide results are mixed. The seven-transcript panel is a secondary descriptive summary of all six animals. |
| 3. Broader interferon estimates depended on the estimator and animal inclusion | Figure 3A–B | Table S4 | Direct expression, background-adjusted scores, and Hallmark enrichment are distinct estimands. One-animal omissions are influence diagnostics, not replication or independent tests. |
| 4. Same-cohort sensitivity analyses produced similar selected estimates | Figure 3C | Figure S3; Tables S4 and S6 | Cell-inclusion and processed-count sensitivities reuse the same six libraries and provide no additional biological replication. |
| 5. External 6-h and 48-h experiments showed different cross-sectional exposure contrasts | Figure 4 | Table S5 | GSE283401 is external context, not replication, validation, a longitudinal trajectory, or an isolated effect of time. |

## 3. Frozen main-figure architecture

### Figure 1 — Study design and limits of the reference partition

- **A:** source design and the three analysis populations: 7,461 cells for the
  partition audit, 7,367 non-Rare cells for conditional composition, and 7,371
  Scrublet-negative cells for molecular analyses.
- **B:** UMAP with numeric reference partitions 0–6.
- **C:** fixed-graph Leiden random-seed sensitivity.
- **D:** six-animal conditional partition percentages and group contrasts.
- **Action:** retain data and layout. Re-render exact *P* values with leading
  zeros. Do not add a new analysis.

### Figure 2 — Transcriptome-wide results and selected-transcript expression

- **A:** MA plot for all 13,926 DESeq2-tested features. Distinguish the 40
  features meeting BH FDR < 0.05, comprising 11 with positive and 29 with
  negative log2 fold changes. Identify the seven selected transcripts without
  implying gene-level multiplicity-adjusted significance. The five features
  meeting BH FDR < 0.05 with the same direction in all six omission refits may
  be marked as influence-stable threshold hits; *Ttr* requires its existing
  process/ambient-expression caveat in the legend.
- **B:** all 42 direct-expression values for the seven selected transcripts
  across C1–C3 and S1–S3, with every animal identifiable.
- **C:** all six equal-weight panel values and group means.
- **D:** full-cohort and six one-animal-omission panel mean differences.
- **Disposition of former Figure 2B:** the seven transcript-specific DESeq2
  estimates remain in Table S3 and are removed from the main artwork.
- **Action:** rebuild from existing validated Table S3 and panel-ready files.
  No differential-expression model is to be refitted for this redesign.

### Figure 3 — Estimator, animal-influence, and cell-inclusion sensitivity

- **A:** background-adjusted `score_genes` values and full/omission contrasts.
- **B:** focused Hallmark results for the full cohort and all animal omissions.
- **C:** 7,371-versus-7,461 cell-inclusion sensitivity.
- **Action:** retain numerical panels. Synchronize terminology and legend with
  the rewritten Methods. Direct expression, module scores, and NES remain on
  separate scales.

### Figure 4 — External cross-experiment context

- **A:** differences between GSE267933 and GSE283401, including independent
  animals at the two GSE283401 endpoints.
- **B:** separate old-animal 6-h and 48-h gene-level treatment contrasts.
- **C:** formal 48-h-versus-6-h Hallmark interaction contrasts.
- **Action:** retain numerical panels. Use `carrier-gas control` for the
  external cohort where space permits and preserve the no-replication,
  no-validation, and no-trajectory boundary.

## 4. Main and supplementary tables

- **Tables 1–2:** retain values and structure.
- **Table S3:** remains the complete source for the 13,926-feature DESeq2
  result, the 40 full-cohort FDR hits, the seven selected-transcript estimates,
  animal values, and all one-animal-omission summaries.
- **Table S4:** retains complete all-50-Hallmark and estimator-sensitivity
  results.
- **Table S5:** retains the GSE283401 gene and Hallmark interaction results.
- **Table S6/Figure S3:** retain GSE289098 only as same-cohort processed-count
  sensitivity.

## 5. Global display-language rules

- Use American English: `neighborhood`, `neighbor`, and `gray`.
- Use leading zeros for decimal values: `0.05`, `0.20`, and `P = 0.70`.
- Use uppercase panel labels A–D.
- Define `log2 counts per million` before using `log2 CPM` in a legend or axis.
- Describe the seven-transcript display neutrally as a secondary descriptive
  analysis of transcripts selected during exploratory inspection. Detailed
  reviewer-facing provenance belongs in the response letter.
- Never describe GSE289098 or GSE283401 as replication or validation.
- Never convert animal-omission rows into independent hypothesis tests.
- Do not add the mean-UMI equalization audit to the manuscript displays. Its
  formal result was indeterminate and it remains an internal sensitivity
  record.

## 6. Construction and writing order

1. Rebuild and visually verify Figure 2.
2. Apply the global language and formatting rules to Figures 1, 3, 4 and the
   supplementary displays.
3. Export vector PDF masters and 300-dpi LZW TIFF files from R; retain the R
   scripts, plot objects, source/output hashes, warnings, and session records.
4. Rewrite Results from the frozen displays. Results text must not be used to
   justify a later display change unless this map is versioned explicitly.
