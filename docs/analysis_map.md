# Analysis and display map

All final plots were produced in R. Python is used only for selected data
extraction, matrix-identity auditing, and independent checks.

| Manuscript item | Biological/data role | Frozen input | Primary code | Final output / table |
|---|---|---|---|---|
| Figure 1 | GSE267933 design, submitted partition, seed sensitivity, composition | `analysis/figure1/data/panel_ready` | `analysis/figure1/scripts/79_plot_glia_R1_figure1.R` | `analysis/figure1/outputs/Figure1_study_design_partition_audit.*` |
| Table 1 | Per-animal cell/UMI audit | Figure 1 panel data | `analysis/tables1_2/scripts/80_prepare_glia_R1_tables1_2.R` | `analysis/tables1_2/outputs/Table1_*` |
| Table 2 | Per-animal partition composition | Figure 1 panel data | `analysis/tables1_2/scripts/80_prepare_glia_R1_tables1_2.R` | `analysis/tables1_2/outputs/Table2_*` |
| Figure 2 | Full pseudobulk ranking and selected-transcript summaries | `analysis/figure2/data/panel_ready`; Table S3D | `analysis/figure2/scripts/84_plot_glia_R1_figure2_rebalanced.R` | `analysis/figure2/outputs/Figure2_transcriptome_and_selected_transcripts.*` |
| Table S3 | Complete primary pseudobulk and selected-transcript results | Frozen source-audit outputs | `analysis/figure2/scripts/82_prepare_glia_R1_figure2_data.R`; `95_build_glia_R1_tableS3_primary_transcriptome.R` | `analysis/figure2/outputs/tableS3` |
| Figure 3 | Estimator dependence, animal influence, doublet-inclusion sensitivity | `analysis/figure3/data/panel_ready` | `analysis/figure3/scripts/85_plot_glia_R1_figure3.R` | `analysis/figure3/outputs/Figure3_animal_influence_estimator_sensitivity.*` |
| Table S4 | Complete Hallmarks, score sensitivity, and component decomposition | Frozen source-audit outputs | `analysis/source_audit/scripts`; final table assembly is frozen in Parts A–K | `analysis/figure3/outputs/tableS4` |
| Figure 4 | GSE283401 external cross-experiment context | `analysis/figure4/data/panel_ready` | `analysis/figure4/scripts/86_prepare_glia_R1_figure4_data.R`; `87_plot_glia_R1_figure4.R` | `analysis/figure4/outputs/Figure4_independent_temporal_context.*` |
| Table S5 | GSE283401 sample, gene, and Hallmark results | Public GSE283401 source matrices | `analysis/figure4/scripts/86_prepare_glia_R1_figure4_data.R` | `analysis/figure4/outputs/tableS5` |
| Figure S1 | Partition annotation, detection, and predicted-doublet audit | `analysis/figureS1/data/panel_ready` | `analysis/figureS1/scripts/91_prepare_glia_R1_figureS1_tableS1.py`; `92_plot_glia_R1_figureS1.R` | `analysis/figureS1/outputs/FigureS1_*` |
| Table S1 | Partition/marker audit | Frozen reference-object exports | Figure S1 preparation script | `analysis/figureS1/outputs/tableS1` |
| Figure S2 | Partition stability diagnostics | `analysis/figureS2/data/panel_ready` | `analysis/figureS2/scripts/91_prepare_glia_R1_figureS2_tableS2.R`; `92_plot_glia_R1_figureS2.R` | `analysis/figureS2/outputs/FigureS2_*` |
| Table S2 | Complete seed, subsampling, and resolution diagnostics | Frozen stability exports | Figure S2 preparation script | `analysis/figureS2/outputs/tableS2` |
| Figure S3 | GSE289098 same-cell processed-count sensitivity | `analysis/figureS3/data` | `analysis/figureS3/scripts/88_preflight_gse289098_common_payload.py`; `89_prepare_figureS3_tableS6.R`; `90_plot_figureS3.R` | `analysis/figureS3/outputs/FigureS3_*` |
| Table S6 | GSE289098 identity, mapping, and estimate comparisons | Same-cell preflight outputs | Figure S3 preparation script | `analysis/figureS3/outputs/tableS6` |

The historical workspace used an internal Figure S6 label for the GSE289098
sensitivity analysis. Repository filenames, metadata, and scripts have been
renumbered to its final manuscript designation, Figure S3. Table S6 retains its
final manuscript number.

