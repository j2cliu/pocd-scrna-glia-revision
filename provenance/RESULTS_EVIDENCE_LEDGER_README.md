# GLIA revised Results evidence ledger

This ledger is the paragraph-level source of truth for drafting the rebuilt Results. It locks each paragraph to a claim ceiling, exact numeric anchors, machine-readable source rows, canonical scripts, display cross-references, required caveats, and prohibited language.

## Drafting rules

1. Draft only the paragraph identified by `ledger_id`; do not import prose from the submitted Results.
2. Recheck every number against `source_data` using `row_selector` before it enters prose.
3. Keep direct pseudobulk log2-CPM, DESeq2 log2 fold change, `score_genes` difference, and Hallmark NES as separate estimands.
4. Retain S3 in all primary full-cohort estimates. Omission rows are influence diagnostics, not grounds for exclusion.
5. Use `combined sevoflurane-plus-laparotomy exposure` or `combined exposure` versus `oxygen control` for GSE267933. Source labels may remain only in code or machine-readable provenance fields.
6. Treat submitted numeric partitions as traceability partitions rather than validated biological states.
7. Treat GSE289098 as a same-cohort processed-count sensitivity analysis and GSE283401 as bounded external biological context.
8. If a proposed sentence exceeds `permitted_claim`, contradicts `required_caveat`, or uses any `prohibited_wording`, revise or omit it.

## Hallmark sign-count guardrail

Across all 50 Hallmarks, 30 retained the same NES sign across the full fit and all six omission refits: 27 were consistently negative and three consistently positive. Therefore, “30 retained the same sign” and “27 remained negative” are different statements and must not be interchanged.

## Scaffold map

- Section 1: R1.P1–R1.P4
- Section 2: R2.P0–R2.P3
- Section 3: R3.P1–R3.P3
- Section 4: R4.P1–R4.P3
- Section 5: R5.P1–R5.P3

R4.P1 is retained in the ledger as a historical record of the immediate-early-response nuisance analysis, but it was retired from the manuscript when the corresponding Table S4 material was removed. It is therefore excluded from the active scaffold-ID check.

The matching empty manuscript scaffold is `revision_work/manuscript_staging/results_glia_R1_rebuild.md`.
