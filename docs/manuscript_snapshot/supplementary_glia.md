# Supplementary Material

## Supplementary Methods

These procedures provide the implementation detail supporting the analyses described in the Main Methods. Dataset roles, biological units, and inferential boundaries are stated here where needed to make each procedure self-contained.

### SM1. Dataset provenance and cohort roles

This study was a secondary, exploratory analysis of publicly available transcriptomic data and involved no new animal experiment. The primary dataset, GSE267933 (Suo et al., 2025), comprises hippocampal single-cell RNA sequencing from 18-month-old male C57BL/6 mice assigned to sevoflurane anesthesia plus exploratory laparotomy or an oxygen-control condition. The source report states that mice were randomly assigned to the two groups but does not describe the allocation procedure. Exposed animals received 2.5% sevoflurane in 50% oxygen for 30 min plus exploratory laparotomy (n = 3). Control animals received the same 50% oxygen exposure without anesthesia or laparotomy (n = 3). Hippocampi were collected 24 h after the procedure and processed with the 10x Genomics Chromium Single Cell 3′ v2 platform. Repository metadata use the source labels `Surgery` and `Control`; these labels were retained in machine-readable fields, whereas text and displays refer to combined sevoflurane-plus-laparotomy exposure and oxygen control. Because the experiment lacked an anesthesia-only arm, the contrast cannot distinguish effects of anesthesia, laparotomy, or their interaction.

GSE289098 (Zheng et al., 2025) contains a Cell Ranger v3.0.2-aggregated processed-count matrix derived from the same six source libraries as GSE267933. Sample-aware matching of 10x barcodes and source-library identifiers established one-to-one correspondence across all 20,684 Cell Ranger-filtered barcodes. The 27,998 Ensembl feature identifiers and corresponding gene-symbol rows also occurred in the same order. GSE289098 therefore contributed no additional animals or libraries and was used only to examine an alternative processed-count construction. It was not treated as biological replication, a technical replicate, or orthogonal validation.

GSE283401 (Yin et al., 2024) provided external context from FACS-isolated hippocampal microglial bulk RNA sequencing. The accession contained 48 hippocampal and 48 hemisphere libraries from male C57BL/6 mice; only hippocampal libraries were considered. Formal models were restricted to 32 old mice aged 20–22 months: eight source-labeled controls and seven isoflurane-plus-laparotomy animals at 6 h, and eight controls and nine exposed animals at 48 h. Different animals contributed at the two sampling times. Exposed mice received 1.2% isoflurane in 30% oxygen balanced with nitrogen for 2 h plus laparotomy. Controls received the carrier gas for 2 h and the same bupivacaine and meloxicam regimen, but no anesthesia or laparotomy. The 6- and 48-h endpoints were measured from the start of anesthesia or carrier-gas administration. The 16 young-animal hippocampal libraries aged 3–5 months, with four animals in each time-by-treatment cell, were excluded from the formal models. GSE283401 differs from GSE267933 in age, anesthetic, exposure duration, assay, laboratory, and sampling design and was used as external cross-experiment context rather than replication, validation, or a longitudinal time course.

### SM2. Primary single-cell preprocessing and object provenance

The six GSE267933 feature-barcode matrices were concatenated in Scanpy (Wolf et al., 2018), with source-library and exposure labels retained for every barcode. Cells expressing fewer than 200 or more than 5,000 genes, or with a mitochondrial UMI fraction greater than or equal to 20%, were excluded. Genes detected in fewer than three cells were also removed. The resulting quality-control-filtered whole-cell object contained 18,114 cells and 18,936 genes.

A separate pre-quality-control AnnData object, `adata_raw.h5ad`, retained integer UMI counts for all 20,684 Cell Ranger-filtered barcodes and 27,998 features. Count-based analyses selected verified analysis barcodes from this integer matrix and then applied the gene filter specified for the corresponding analysis. For visualization and cell-level summaries, quality-control-filtered counts were normalized to 10,000 counts per cell and log1p-transformed. The complete log-normalized matrix was stored in the processed object's `.raw` attribute. Despite the attribute name, `.raw.X` therefore contains normalized expression rather than integer counts and was not used for count-based inference.

Highly variable genes were selected with Scanpy's Seurat-dispersion procedure using a minimum mean of 0.0125, maximum mean of 3, and minimum normalized dispersion of 0.5, yielding 2,722 genes. Values were scaled to unit variance and clipped at 10. Principal-component analysis generated 50 components. The whole-cell neighborhood graph used 15 neighbors and the first 20 components. UMAP (McInnes et al., 2018) and Leiden clustering (Traag et al., 2019) were run with random state 0, using a whole-cell Leiden resolution of 0.5. Four whole-cell clusters supported by canonical microglial markers, including *P2ry12*, *Tmem119*, *Cx3cr1*, and *Csf1r*, defined the fixed 7,461-cell microglial population. This fixed population and its stored coordinates were the inputs to the partition and annotation diagnostics; those diagnostics did not regenerate cell selection or whole-cell principal components. The stored processed object preserved cell-for-cell continuity because the available executable record did not contain every intermediate output required to regenerate the whole-cell object and reference partition exactly.

Three microglial analysis populations had distinct uses. The reference partition and annotation audit included all 7,461 fixed microglia. Partition-conditioned composition excluded partition 6 and therefore used 7,367 non-Rare cells. Molecular analyses used 7,371 microglia not called as doublets by Scrublet, irrespective of reference-partition identity. The complete 7,461-cell population was also used for cell-inclusion sensitivity analyses where noted.

### SM3. Reference microglial partition, annotation, and stability

The reference within-microglia partition was generated without recomputing microglia-specific highly variable genes or principal components. A 15-neighbor graph was constructed from the first 20 principal components inherited from the whole-cell analysis. Leiden clustering at resolution 0.4 with random state 0 produced seven numeric partitions. Resolution 0.4 was the value encoded in the stored analysis; no prespecified or documented data-driven selection criterion for this value was available. The associated traceability labels were partition 0, Inflammatory; partition 1, Transitional-A; partition 2, Homeostatic-A; partition 3, Homeostatic-B; partition 4, Homeostatic-C; partition 5, Transitional-B; and partition 6, Rare. Because the biological names were not generated by a prespecified reproducible annotation rule, analyses used numeric partition identifiers; the names served only as traceability labels and were not treated as validated biological states.

For the annotation audit, integer counts from all 7,461 reference microglia were mapped to gene symbols. When multiple features mapped to the same symbol, the first feature was retained. Genes without a detected count in this microglial population were removed, leaving 17,878 symbol-deduplicated genes. Counts were normalized once to 10,000 per cell and log1p-transformed. Each numeric partition was compared with all remaining microglia using Scanpy's tie-corrected one-versus-rest Wilcoxon procedure across the complete retained gene set. Benjamini–Hochberg adjustment (Benjamini and Hochberg, 1995) was applied within each partition. Table S1B retains every marker with an adjusted *P* value below 0.05 and a positive Scanpy log2 fold-change estimate. The `top_50_marker` field identifies the first 50 retained markers per partition for inspection and is not an additional retention criterion. These per-cell statistics were descriptive annotation diagnostics and did not constitute biological-replicate inference.

The descriptive dot plot used the author-selected transcripts *P2ry12*, *Tmem119*, *Cx3cr1*, *Tnf*, *Il1b*, *Irf7*, *Ifitm3*, *Nr1d1*, and *Dbp* rather than a set of top differentially expressed genes. Point color was the within-partition mean after total-count normalization to 10,000 and log1p transformation. Point size was the raw-count detection percentage, defined as 100 times the number of partition cells with an integer UMI count greater than zero divided by the total number of cells in that partition. Partition 6 was retained in the marker and raw-count detection audits for traceability but excluded from partition-level biological interpretation and from partition-conditioned composition. Molecular analyses did not use partition identity as an inclusion rule: cells not called as doublets were retained irrespective of numeric partition, and predicted doublets were excluded irrespective of partition.

Partition sensitivity was examined using two distinct computational objects. First, Leiden clustering at resolution 0.4 was repeated with seeds 0, 1, 7, 42, 123, and 2024 while holding fixed the exact neighborhood graph stored in `adata_microglia_subtyped.h5ad`. The number of partitions and adjusted Rand index relative to the reference numeric partition were recorded for each run. All 15 unique pairwise adjusted Rand indices among the six runs were also calculated. This procedure isolated algorithmic seed sensitivity on the stored graph and did not resample cells or animals.

Second, conditional graph-reconstruction diagnostics sampled 5,969 of the 7,461 cells, corresponding to 80% of the population, without replacement while retaining the first 20 inherited principal-component coordinates. Cell sampling used `numpy.random.default_rng(0)`. The random-number generator was reinitialized for the variable-seed and fixed-seed arms, so the first 50 cell samples were shared between them. A 15-neighbor graph was rebuilt for each sample. At Leiden resolution 0.4, 100 runs varied the neighborhood-graph and Leiden seeds together from 0 through 99. A separate 50-run series fixed both graph and Leiden seeds at 0 while otherwise applying the same conditional subsampling framework. Scanpy 1.11.5 defaults were used for `scanpy.pp.neighbors` and `scanpy.tl.leiden` beyond the arguments stated here; the graph and partition engines were leidenalg 0.11.0 and python-igraph 1.0.0. For each reference partition, the best-match Jaccard index was the largest intersection-over-union with any reconstructed partition among the sampled cells. Its empirical mean, sample standard deviation, and median were recorded. Recovery frequency was the proportion of the 100 variable-seed runs in which the best-match Jaccard index was strictly greater than 0.5. A mean-Jaccard indicator of at least 0.75 was retained as a descriptive reference rather than a biological validation threshold.

Resolution sensitivity was assessed separately. A full-cohort 15-neighbor reference graph was constructed from the inherited coordinates with graph seed 0, and Leiden clustering was run with seed 0 at resolutions 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, and 1.0. At each resolution, 30 corresponding 80% cell samples were analyzed, with graph and Leiden seeds varied together from 0 through 29. `numpy.random.default_rng(0)` was reinitialized at each resolution, giving the same sequence of 30 cell samples across resolutions. For every full-cohort partition at a given resolution, the mean best-match Jaccard index across the 30 reconstructed samples was calculated. Table S2F reports these per-partition means. Table S2E reports their across-partition mean, median, minimum, and maximum; the modal number of reconstructed partitions across the 30 samples; and counts of full-cohort partitions with mean Jaccard indices of at least 0.75 and at least 0.60. Reconstructed partition identifiers were local to each resolution and were not treated as corresponding biological states across resolutions. These procedures rebuilt neighborhood graphs but did not repeat cell selection, normalization, highly variable-gene selection, principal-component analysis, or animal sampling. They therefore assessed conditional partition sensitivity rather than end-to-end pipeline reproducibility, biological replication, or biological state stability. No central exposure or interferon analysis was conditioned on distinctions among the six non-Rare reference partitions.

For partition-conditioned composition, the count and percentage of each of the six non-Rare reference labels were calculated separately for every animal. The denominator was the animal's total number of non-Rare microglia. Combined-exposure-minus-control differences were summarized as mean percentage-point differences with Welch–Satterthwaite 95% confidence intervals, Cohen's *d*, and two-sided exact label-permutation *P* values across all 20 allocations of three of the six animals to the exposed group. Perfect animal-level separation was defined descriptively as every value in one group exceeding every value in the other group. These summaries were descriptive and were not used to establish reproducible state expansion.

### SM4. Doublet and contamination assessments

Scrublet (Wolock et al., 2019) was applied separately to each of the six quality-control-filtered integer UMI matrices. The expected doublet rate was 0.06, the random state was 42, and automatically selected sample-specific thresholds ranged from 0.282 to 0.327. Predicted-doublet calls were mapped to the fixed microglial barcode set. Calls were treated as model-based predictions rather than confirmed doublets, and the fixed whole-cell annotations were held unchanged.

The non-Rare reference partitions defined the 7,367-cell composition population, whereas predicted-doublet status defined the 7,371-cell molecular-analysis population. The latter rule excluded 90 predicted doublets from the 7,461-cell reference population without using partition identity. Marker and raw-count detection audits retained all seven partitions for traceability. Partition identity was therefore not used as a surrogate doublet rule for molecular analyses.

Unfiltered droplet matrices were unavailable for GSE267933, precluding estimation of an ambient-RNA profile from empty droplets. No ambient-RNA correction was applied. Raw-count detection was screened descriptively within each numeric partition for microglial-core transcripts (*Csf1r*, *Cx3cr1*, *Hexb*, *P2ry12*, and *Tmem119*), border- or CNS-associated macrophage transcripts (*Cd163*, *F13a1*, *Lyve1*, *Mrc1*, and *Pf4*), the choroid-plexus-associated transcript *Ttr*, and oligodendrocyte-lineage transcripts (*Cldn11*, *Mal*, *Mbp*, *Mobp*, and *Plp1*). For each transcript and partition, detection was 100 times the number of cells with an integer UMI count greater than zero divided by the total number of cells in that partition. The screen included all seven reference partitions. It was not converted into a module score, used for ambient-RNA correction, or interpreted as establishing cellular purity.

### SM5. Pseudobulk construction and animal-influence analyses

Integer UMI counts from `adata_raw.h5ad` were restricted to the 7,371 Scrublet-negative microglial barcodes, comprising 3,323 oxygen-control and 4,048 combined-exposure cells. Counts were summed within animal/library, so the six pseudobulk libraries rather than individual cells were the biological observations. Genes with fewer than 10 total counts across the six libraries were excluded, defining a 13,926-feature primary analysis universe.

Differential expression was fitted with DESeq2 (Love et al., 2014) using `design = ~ group`. Control was the reference level, and the source-labeled Surgery-versus-Control coefficient represented combined sevoflurane-plus-laparotomy exposure relative to oxygen control. Fitted outputs retained unshrunk log2 fold changes, Wald standard errors, nominal *P* values, and Benjamini–Hochberg-adjusted *P* values. Wald 95% confidence intervals were calculated as log2 fold change ± 1.96 times the Wald standard error. This coefficient does not isolate anesthesia, laparotomy, or their interaction.

The complete differential-expression and enrichment workflow was repeated using all 7,461 reference microglia as a cell-inclusion sensitivity analysis. Comparisons between the Scrublet-negative and all-cell analyses included the sign of the direct seven-transcript panel contrast and, for the interferon-alpha-response, interferon-gamma-response, TNFA-signaling-via-NF-κB, and inflammatory-response Hallmarks, the normalized enrichment score sign and whether FDR < 0.05 under the enrichment procedure in SM8. These comparisons were made for the full cohort and each systematic animal omission. The maximum absolute normalized enrichment score difference between cell-inclusion rules was also recorded.

Animal-level influence was evaluated by refitting the DESeq2 model six times, each time omitting one animal. The 13,926-feature universe defined by the full-cohort filter was held fixed for every omission fit rather than refiltered after omission. Median-ratio size factors, model parameters, Wald tests, default DESeq2 independent filtering, and Benjamini–Hochberg adjustment were recalculated separately in each refit. Adjusted *P* values unavailable after independent filtering remained `NA`. These leave-one-animal-out estimates were treated as influence diagnostics rather than independent hypothesis tests, and no animal was removed from the full-cohort analysis. For each library, the audit retained the starting and Scrublet-negative microglial cell counts, the number of predicted doublets removed, retained pseudobulk UMI total, median UMI per retained cell, and percentages of cells and UMIs removed. These quantities characterize potential technical influence but cannot distinguish biological heterogeneity from unrecorded processing or sampling-time variation.

### SM6. Direct-expression and background-adjusted scoring sensitivity

#### Direct expression of the selected transcript panel

A secondary descriptive analysis summarized seven interferon-responsive transcripts selected during exploratory inspection of the primary dataset: *Irf7*, *Ifitm3*, *Isg15*, *Mx1*, *Ifit1*, *Ifit2*, and *Ifit3*. The panel was used to summarize animal-level expression patterns.

The 7,371 Scrublet-negative microglia constituted the primary cell population without conditioning on reference partition or Rare label. Feature rows mapping to the same gene symbol were summed. For gene *g* in animal *j*, direct expression was calculated as

`x_jg = log2[((Y_jg + 0.5)/(N_j + 1)) × 10^6]`,

where `Y_jg` is the animal-level summed UMI count for the gene and `N_j` is the animal's total microglial UMI count across all 27,998 feature rows before gene filtering or symbol collapse. The animal-level panel value was the unweighted mean of the seven gene-level values, preventing more highly expressed genes from receiving greater weight.

The analysis retained all six animal values, the combined-exposure-minus-control mean difference, a Welch–Satterthwaite 95% confidence interval, Cohen's *d*, and the exact two-sided label-permutation *P* value obtained from all 20 allocations of three animals to each group. For each systematic leave-one-animal-out contrast, the mean difference, Welch–Satterthwaite confidence interval, and pooled-standard-deviation Cohen's *d* were recalculated on the retained three-versus-two animals; exact permutation *P* values were left `NA`. These rows were influence diagnostics. The calculation was repeated using all 7,461 reference microglia as a cell-inclusion sensitivity analysis.

#### Background-adjusted module scoring

Background-adjusted module scores were evaluated separately because they estimate expression relative to an algorithmically matched control set rather than direct pseudobulk abundance. Ensembl features were mapped to gene symbols, duplicate symbols were resolved by retaining the first feature, and the resulting 27,933-feature matrix was normalized once to 10,000 counts per cell and log1p-transformed before removal of the Scrublet-predicted doublets.

`scanpy.tl.score_genes` was applied to the 7,371 retained cells without conditioning on reference partition, using `use_raw = False`, the complete symbol-deduplicated matrix as the gene pool, and `ctrl_as_ref = True`. Cell-level scores were averaged within animal before group comparison.

Parameter sensitivity was evaluated for the seven-transcript panel. Seeds 0–19 were evaluated with `ctrl_size = 50` and `n_bins = 25`. Control-set sizes of 25, 50, 100, and 200 were evaluated with seed 42 and `n_bins = 25`. Bin numbers of 10, 25, and 50 were evaluated with seed 42 and `ctrl_size = 50`. The seed-42, `ctrl_size = 50`, `n_bins = 25` tuple appears in both the control-size and bin-number sweep tables by design and does not represent an additional analysis. For every configuration, the full three-versus-three animal values, group means, mean difference, Welch–Satterthwaite interval, Cohen's *d*, and exact permutation *P* value were retained together with the mean difference and Cohen's *d* after omitting S3. The omit-S3 summaries assessed whether the sign change identified by the systematic fixed-score animal-omission analysis was preserved across scoring configurations; they were not treated as independent experiments.

Seed 42, `ctrl_size = 50`, and `n_bins = 25` defined the fixed configuration used for complete animal-influence assessment. Its six animal means were used for the full-cohort contrast and all six systematic animal omissions. Animal-omission contrasts were calculated from the fixed six-animal score vector; neither `score_genes` nor matched-control selection was rerun after omitting an animal. The mean difference, Welch–Satterthwaite interval, and pooled-standard-deviation Cohen's *d* were recalculated on the retained animals, whereas exact permutation *P* values were not assigned to unequal-group omission contrasts.

The exact 249-gene matched-control manifest is provided in Table S4I with matrix indices, expression bins, mean log1p-normalized expression, scoring parameters, and the control-set SHA-256. Cell-level selected-gene and matched-control means were retained separately and averaged within animal, and their difference was checked against the stored Scanpy score. The matched-control component was treated solely as an algorithmic reference term rather than a biological background program and is not on the scale of the direct animal-level expression measure.

Seven fixed-configuration leave-one-transcript-out analyses were also performed. For each six-transcript panel, `score_genes` was rerun and the matched controls were selected for that altered panel using seed 42, `ctrl_size = 50`, and `n_bins = 25`. This rematching differs from the animal-omission calculations, which reused the fixed seven-transcript score vector and matched-control set. Each transcript-omission analysis retained the full-cohort and omit-S3 contrasts. These analyses assessed dependence on panel membership rather than robustness across biological replicates. Direct pseudobulk expression, background-adjusted module scores, and Hallmark enrichment statistics remained distinct estimands and were not combined.

### SM7. GSE289098 same-cohort processing sensitivity

GSE289098 provides a Cell Ranger v3.0.2-aggregated processed-count matrix derived from the same six source libraries as GSE267933. Sample-aware matching used each 10x barcode together with its source-library identifier across the 20,684 Cell Ranger-filtered barcodes. The ordered alignment retained all 27,998 Ensembl feature identifiers and corresponding gene-symbol rows. Because both accessions represent the same animals and libraries, comparisons between them assess an alternative processed-count construction only and do not provide biological replication, a technical replicate, or orthogonal validation.

The exact 7,371 Scrublet-negative microglial barcodes and their fixed animal and group assignments were transferred from the primary analysis by one-to-one sample-aware matching. GSE289098-specific quality-control thresholds, microglial reclassification, clustering, and subtype selection were not applied. The complete ordered 27,998-feature universe was retained for count alignment and for each matrix's all-feature UMI denominators. The exact GSE289098 barcode-suffix-to-animal crosswalk is provided in Table S6B. Per-entry count differences were audited across the exact matched cell-by-feature payload. Per-animal UMI retention was defined as 100 times the GSE289098 all-feature UMI total divided by the GSE267933 all-feature UMI total on the same fixed 7,371-cell set.

For each processed-count matrix, integer counts were summed within animal/library. Gene-level direct-expression values and the equal-weight seven-transcript panel mean were recalculated using the formula specified in SM6. The full-cohort combined-exposure contrast was summarized by the mean difference, Welch–Satterthwaite confidence interval, Cohen's *d*, and the exact two-sided *P* value from all 20 three-versus-three label allocations. For the six leave-one-animal-out scenarios, mean differences, Welch–Satterthwaite intervals, and Cohen's *d* values were recalculated on the retained animals; exact permutation *P* values were left `NA`. Sign agreement was evaluated for the panel and each transcript across the full cohort and the six animal-omission scenarios. Differences between the two matrices were treated as deterministic paired descriptions of alternative processing outputs from the same libraries and were not assigned paired-test *P* values.

As a secondary gene-level sensitivity analysis, DESeq2 was fitted separately to both count matrices with the same `~ group` design and the common 13,926-feature universe defined by the primary GSE267933 total-count filter. Median-ratio size factors and DESeq2 independent filtering were calculated independently within each matrix. Comparisons used unshrunk log2 fold changes and Wald confidence intervals under the common design and feature universe. Matrix-specific Benjamini–Hochberg-adjusted *P* values were not used as criteria for equivalence or concordance. This procedure evaluated sensitivity to processed-count construction within the same biological material and did not test cross-cohort generalizability.

### SM8. Hallmark enrichment procedures

Pathway analyses for GSE267933 and the formal GSE283401 time-by-treatment interaction used the mouse MSigDB 2026.1.Mm Hallmark collection (Liberzon et al., 2015; Castanza et al., 2023), obtained with msigdbr 26.1.0 using `db_species = "MM"`, `species = "Mus musculus"`, and `collection = "MH"`. This MSigDB collection contains mouse-ortholog versions of the Hallmark gene sets; we performed no additional human-to-mouse ortholog mapping. Ensembl identifiers were mapped to gene symbols with org.Mm.eg.db 3.22.0 using `multiVals = "first"`.

For GSE267933, genes were ranked by the signed DESeq2 Wald statistic from the animal-level pseudobulk model. When multiple Ensembl features mapped to the same gene symbol, the feature with the largest absolute Wald statistic was retained. Preranked enrichment (Subramanian et al., 2005) was performed with `clusterProfiler::GSEA` (Wu et al., 2021) using the fgsea engine (fgsea 1.36.2), `minGSSize = 10`, `maxGSSize = 500`, `pvalueCutoff = 1`, `pAdjustMethod = "BH"`, `eps = 1 × 10^−30`, `seed = TRUE`, and `verbose = FALSE`. R's random seed was set to 42 immediately before every call. Benjamini–Hochberg adjustment was applied across all 50 Hallmark sets. Positive normalized enrichment scores denoted ranking toward combined exposure and negative scores toward oxygen control. The four Hallmarks displayed in focused summaries were interferon alpha response, interferon gamma response, TNFA signaling via NF-κB, and inflammatory response; display selection did not alter the 50-set multiple-testing family.

Animal-influence analyses retained the 13,926-feature universe defined by the primary full-cohort count filter. The DESeq2 model and enrichment procedure were repeated after systematically omitting each of the six animals. These leave-one-animal-out estimates were treated as influence diagnostics rather than independent hypothesis tests. The GSE283401 interaction analysis used the same Hallmark collection, identifier-mapping rule, duplicate-symbol rule, enrichment parameters, and 50-set multiple-testing family, with genes ranked by the signed interaction Wald statistic defined in SM9.

### SM9. GSE283401 external-context models

The sample inventory and cohort restrictions are described in SM1. The deposited treatment labels `surgery` and `control` were retained in machine-readable source fields; `control` was the reference level, and text describes `surgery` as isoflurane plus laparotomy. Version suffixes were removed from Ensembl identifiers before mapping with org.Mm.eg.db. For time-specific gene-level analyses, the deposited 6-h and 48-h count matrices were fitted separately with DESeq2 models of the form `~ treatment`. The 6-h analysis retained 26,317 Ensembl rows with at least 10 total counts across its 15 old-animal libraries. The 48-h analysis independently retained 23,732 rows meeting the same criterion across its 17 old-animal libraries. A positive log2 fold change denoted higher expression after isoflurane plus laparotomy. Unshrunk log2 fold changes and Wald 95% intervals, calculated as the log2 fold change ± 1.96 times its standard error, were retained. DESeq2 median-ratio normalization, default independent filtering, Wald testing, and Benjamini–Hochberg adjustment were performed separately within each time-specific fit; adjusted *P* values unavailable after independent filtering remained `NA`. For each of the seven displayed gene symbols, when multiple Ensembl rows mapped to the same symbol, the row with the largest `baseMean` within that time-specific estimand was selected for display. These analyses represent separate cross-sectional contrasts and were not treated as paired or longitudinal comparisons.

For the formal time-by-treatment analysis, Ensembl rows shared by the two deposited matrices were combined across all 32 old-animal libraries. Rows with fewer than 10 total counts across the combined libraries were excluded, leaving 30,280 rows. DESeq2 was fitted as `~ time + treatment + time:treatment`, with 6 h and source-labeled control as the reference levels. Unshrunk coefficients and Wald 95% intervals were retained. Default independent filtering and Benjamini–Hochberg adjustment were performed within this interaction fit, separately from both time-specific fits, and unavailable adjusted *P* values remained `NA`. The interaction coefficient was defined as

`(48-h exposed - 48-h control) - (6-h exposed - 6-h control)`.

A positive coefficient therefore denoted a more positive treatment-associated contrast in the deposited 48-h experiment than in the deposited 6-h experiment. Genes were ranked for Hallmark enrichment by the signed interaction Wald statistic. When multiple Ensembl rows mapped to one gene symbol, the row with the largest absolute Wald statistic was retained. Positive normalized enrichment scores denoted gene-set ranking toward the more positive 48-h treatment contrast. Enrichment and multiple-testing correction followed SM8.

Because different animals contributed at 6 h and 48 h and the two count matrices were deposited separately, sampling time is inseparable from deposited matrix and from unrecorded processing or sequencing factors aligned with those matrices. The interaction was therefore interpreted as a between-experiment difference in treatment-associated contrast. It was not interpreted as a within-animal trajectory or as replication or validation of GSE267933.

### SM10. Statistical definitions, software, and reproducibility

The animal/library was the biological unit for central exposure comparisons. Cells were treated as observations nested within animals. Per-cell marker statistics, score distributions, and correlations were used only for descriptive analyses and did not increase the biological sample size. Gene-level inference used animal-level pseudobulk counts, whereas targeted transcript and module summaries were first calculated within each animal and then compared between groups.

For animal-level comparisons, the reported mean difference was the combined-exposure mean minus the oxygen-control mean. Let `s1^2`, `s0^2`, `n1`, and `n0` denote the exposed and control sample variances and group sizes. The Welch standard error was

`SE = sqrt(s1^2/n1 + s0^2/n0)`,

and the Welch–Satterthwaite degrees of freedom were

`ν = (s1^2/n1 + s0^2/n0)^2 / {[(s1^2/n1)^2/(n1 - 1)] + [(s0^2/n0)^2/(n0 - 1)]}`.

The 95% confidence interval was `difference ± t(0.975, ν) × SE`. Cohen's *d* was calculated as `d = (mean_exposure − mean_control)/sp`, where `sp` was the pooled within-group standard deviation,

`sp = sqrt{[(n1 - 1)s1^2 + (n0 - 1)s0^2]/(n1 + n0 - 2)}`.

Exact two-sided label-permutation tests for the full three-versus-three comparison enumerated all `C(6, 3) = 20` allocations of three of the six animals to the exposed group, with extremeness defined by the absolute between-group mean difference. Permuted values tied with or exceeding the observed absolute difference were counted, using a numerical tolerance of `10^−12`; no +1 correction was applied because the allocation space was fully enumerated. Complementary allocations produce equal-magnitude effects of opposite sign, so the minimum attainable two-sided exact *P* value was 2/20 = 0.1. This limit applies only to exact three-versus-three tests and not to DESeq2 Wald tests, external-cohort models, or descriptive per-cell calculations. Perfect group separation was defined descriptively as every retained value in one group exceeding every retained value in the other group. The source report states that mice were randomly assigned but does not describe the allocation procedure. Label-permutation results were therefore treated as descriptive.

Systematic leave-one-animal-out estimates were used to assess influence. No animal was excluded from the corresponding full-cohort estimate, and omission estimates were not treated as independent experiments or hypothesis tests. Interpretation emphasized complete animal-level values, mean differences, confidence intervals where defined, standardized effect sizes, exact permutation results, and consistency across animal omissions. A nonsignificant estimate was not interpreted as evidence of absence or attributed automatically to insufficient power. Between-group estimates were interpreted as associations with the observed exposure contrast and not as causal, anesthesia-specific, laparotomy-specific, interaction-specific, or aging-specific effects.

The computational environment used Python 3.11 with Scanpy 1.11.5, anndata 0.12.10, Scrublet 0.2.3 through `scanpy.pp.scrublet`, scikit-learn 1.8.0, SciPy 1.17.1, pandas 2.3.3, leidenalg 0.11.0, and python-igraph 1.0.0 for the conditional partition analyses. Retained Python manifests report NumPy distribution metadata version 2.0.2; the fixed-score component environment imported NumPy runtime version 2.0.0, and this discrepancy is preserved in its session manifest. R analyses used R 4.5.2 with DESeq2 1.50.2, clusterProfiler 4.18.4, fgsea 1.36.2, msigdbr 26.1.0, AnnotationDbi 1.72.0, and org.Mm.eg.db 3.22.0. Random seeds were fixed and reported for each stochastic procedure rather than imposed as one global seed. Machine-readable tables retain the complete differential-expression, enrichment, parameter-sensitivity, marker, provenance, and influence outputs. The exact matched-control gene manifest is provided in Table S4I. The accompanying reproducibility package contains analysis and figure-generation scripts together with input/output SHA-256 manifests, execution metadata, and session information.

## Supplementary Figure Legends

### Figure S1. Retrospective annotation and contamination audit of the submitted seven-way microglial partition.

**A**, Expression of nine author-selected
descriptive genes across all 7,461 cells in the submitted numeric partitions.
Point color denotes mean log1p counts per 10,000 total counts (CP10K), and point
size denotes the percentage of cells with at least one raw count. These genes
were selected for descriptive display and are not the result of an independent
marker-selection procedure. **B**, Percentage of cells called as potential
doublets by Scrublet in each submitted partition; labels give predicted
doublets/total cells. Partition 6 contained 82 predicted doublets among 94 cells
(87.2%). **C**, Percentage of cells with raw-count detection of microglia-core,
border-/central-nervous-system-associated macrophage (BAM/CAM),
choroid-plexus-associated, and oligodendrocyte-lineage transcripts. BAM/CAM
markers were detected in at most 1.52% of cells in any partition. Similar *Ttr*
detection across partitions does not support a partition-specific choroid
plexus identity, whereas oligodendrocyte-lineage transcripts were most
frequently detected in doublet-enriched partition 6. The dagger and dashed outline mark
partition 6, which was excluded from biological interpretation. Legacy
biological names are retained only as traceability fields in Table S1; none of
the submitted partitions is treated as a validated microglial state.

### Figure S2. Seed sensitivity of the submitted microglial partition and conditional graph-reconstruction diagnostics.

**A**, Pairwise adjusted Rand indices (ARIs) among Leiden partitions obtained with seeds 0, 1, 7, 42, 123, and 2024 at resolution 0.4 on the exact neighborhood graph stored in the submitted 7,461-cell microglial object. The right-hand columns give the ARI of each run against the submitted numeric partition and the number of partitions recovered (*k*). Seed 0 reproduced the submitted labels exactly (ARI = 1.000); the other five seeds recovered six or seven partitions with ARIs of 0.403–0.586 against the submitted labels. Across the 15 unique seed pairs, ARI ranged from 0.365 to 0.714 (mean, 0.529). **B**, Conditional 80% cell-subsampling diagnostic at resolution 0.4. Each run drew 5,969 of the same 7,461 cells, retained the inherited whole-cell PC1–PC20 coordinates, and rebuilt a 15-neighbor graph. Filled purple points and horizontal bars show the mean ± SD of each submitted partition's best-match Jaccard index across 100 runs in which the graph and Leiden seeds varied together from 0 to 99; open points show the mean across 50 runs in which both seeds were fixed at 0. The SD bars describe run-to-run dispersion and are not confidence intervals. Numeric partitions are primary; names in parentheses are retrospective traceability labels. Partition 6/Rare is doublet-enriched and is not interpreted as a biological state. The dashed 0.75 line is a descriptive reference, not a validation threshold. Lower panel, distribution of the number of partitions recovered by the same two subsampling arms; six or seven partitions were recovered in 92% of variable-seed runs and 88% of fixed-seed runs. **C**, Separate conditional resolution diagnostic. A full-cohort reference graph was rebuilt from the inherited PC coordinates with graph and Leiden seed 0 at each resolution from 0.1 to 1.0; 30 80% subsample graphs were then rebuilt per resolution with graph and Leiden seeds varying together from 0 to 29. Purple and open points show the across-partition mean and median best-match Jaccard index; gray vertical ranges show the minimum and maximum per-partition means and are not uncertainty intervals. Values above each resolution give the number of reconstructed partitions with mean Jaccard ≥ 0.75 divided by the full-cohort partition count; the second line of each x-axis label gives the full-cohort partition count/modal subsample count. At resolutions 0.5 and 0.6, two reconstructed partitions met the 0.75 reference; Table S2F shows that these included 993- and 990-cell partitions, respectively, so the pattern is not confined to a small partition. Panels B–C are conditional diagnostics after graph reconstruction from a fixed inherited coordinate system. They are not tests on the submitted stored graph, do not resample animals or reconstruct the full analysis pipeline, and do not validate discrete biological states or prove a transcriptional continuum.

### Figure S3. Same-cohort GSE289098 count reprocessing yields similar animal-level estimates for the selected seven-transcript panel.

GSE289098 is
a single Cell Ranger v3.0.2-aggregated processed-count matrix derived from the
same six source libraries and 20,684 filtered cell barcodes as GSE267933; it
provides no additional biological replication. The exact 7,371
Scrublet-negative microglia defined in the primary analysis were transferred by
one-to-one sample-aware barcode matching, without applying the GSE289098-
specific quality-control, cell-classification, or clustering workflow. **A**,
Percentage of total UMI counts retained in GSE289098 for each animal/library on
this fixed cell set. All 20,684 source cells, all 27,998 feature identifiers,
and all 7,371 analysis cells showed one-to-one correspondence. Across the analysis cells,
GSE289098 contained 56,056,655 UMI counts versus 57,261,942 in the primary
payload (1,205,287 fewer; 2.10%); no cell-by-gene count increased. Blue and
orange denote oxygen-control and combined anesthesia-plus-laparotomy
libraries, respectively. **B**, Animal-level values for the seven selected
interferon-responsive transcripts under the two count matrices. For each animal
and gene (*Irf7, Ifitm3, Isg15, Mx1, Ifit1, Ifit2,* and *Ifit3*), expression
was calculated as log2[((gene UMI + 0.5)/(all-feature UMI + 1)) × 10^6]; the
panel value is the equal-weight mean of the seven gene-level log2
counts-per-million values.
The dashed line is the identity line. **C**, GSE289098-minus-primary change in
gene-level log2 counts per million for each animal. Values are rounded to two decimal places
for display; unrounded values are provided in Table S6C. **D**, Combined-
exposure-minus-oxygen-control difference in the panel value for the full cohort
and each systematic one-animal omission. Open and purple-filled circles show
primary and GSE289098 estimates, respectively; right-hand values are the paired
GSE289098-minus-primary differences. The full-cohort estimates were +1.069
(95% Welch CI, −1.638 to +3.775; Cohen’s *d* = 1.325; descriptive exact
permutation *P* = 0.20) for the primary counts and +1.098 (95% Welch CI,
−1.619 to +3.815; *d* = 1.367; descriptive *P* = 0.20) for GSE289098, a
paired contrast change of +0.030 log2 counts per million.
All six one-animal-omission panel estimates remained positive under both
matrices. Omission analyses are influence diagnostics and are not treated as
independent hypothesis tests. This analysis evaluates sensitivity to
processed-count construction within the same libraries; it is not independent
replication, orthogonal validation, or evidence of robustness across cohorts.

## Supplementary Table Titles and Notes

### Table S1. Reference-partition summary and one-versus-rest marker statistics for GSE267933 microglia.

#### Part A — Reference-partition summary

Part A reports the numeric reference partition, legacy label, total cell count, numbers contributed by oxygen-control and combined-exposure libraries, Scrublet-predicted doublet count and rate, group-specific doublet rates, and interpretive scope. Partition 6 (legacy label, Rare) contained 82 predicted doublets among 94 cells and was excluded from biological interpretation.

#### Part B — One-versus-rest marker statistics

Part B reports all positive one-versus-rest markers with a Benjamini–Hochberg adjusted *P* value below 0.05 and a positive Scanpy log2 fold-change estimate after normalization to 10,000 counts per cell and log1p transformation. Reported fields include rank, gene symbol, Wilcoxon score, Scanpy log2 fold-change estimate, nominal and adjusted *P* values, and detection percentages within the partition and among all remaining cells. `top_50_marker` identifies the first 50 retained markers for inspection and is not an additional selection criterion.

#### General notes

- Numeric partitions 0–6 define the reference solution. Legacy names identify the original annotation and are not treated as validated biological states.
- Marker statistics use cells as observations and serve as descriptive annotation measures rather than animal-level differential-expression tests.
- These results do not establish partition stability, differential abundance, or biological-state validity.
- BAM/CAM denotes border- or central-nervous-system-associated macrophage.

### Table S2. Reference-partition seed sensitivity and conditional inherited-PC graph-reconstruction analyses for GSE267933 microglia.

#### Parts

- **Part A, fixed stored-graph seed results:** partition count and adjusted Rand index relative to the reference partition for six Leiden seeds at resolution 0.4 on the stored neighborhood graph.
- **Part B, fixed stored-graph pairwise adjusted Rand indices:** all 15 unique pairwise adjusted Rand indices among the same six fixed-graph runs.
- **Part C, reference-partition subsampling Jaccard summaries:** per-partition best-match Jaccard results after 80% cell subsampling and graph reconstruction at resolution 0.4, with variable- and fixed-seed analyses reported separately.
- **Part D, subsampling partition-count distribution:** numbers and percentages of runs yielding five through nine partitions in the 100 variable-seed and 50 fixed-seed graph-reconstruction analyses.
- **Part E, resolution-sweep summary:** full-cohort and modal subsample partition counts and across-partition best-match Jaccard summaries at eight Leiden resolutions.
- **Part F, resolution-sweep per-partition results:** size and mean best-match Jaccard index for each reconstructed full-cohort partition at each resolution.
- **Data dictionary:** row unit and purpose of each part.

#### General notes

1. Parts A–B and Parts C–F evaluate different computational objects. Parts A–B hold fixed the neighborhood graph stored in `adata_microglia_subtyped.h5ad`. Parts C–F rebuild 15-neighbor graphs from the first 20 principal components inherited from the whole-cell analysis.
2. The reference population contains 7,461 cells and seven numeric partitions at Leiden resolution 0.4. The names Inflammatory, Transitional-A/B, Homeostatic-A/B/C, and Rare are legacy labels rather than validated biological states. Partition 6 (Rare) contains 94 cells, is doublet-enriched, and is not interpreted as a biological state.
3. An adjusted Rand index of 1 denotes identical partitions. Parts A–B quantify algorithmic agreement among partitions obtained from the same cells and stored graph; they contain no biological-replicate test or *P* value.
4. Each conditional subsampling run draws 5,969 of the 7,461 cells (80%) without replacement and does not resample animals. Cell selection, features, normalization, highly variable genes, and principal components are held fixed, so the procedure does not represent an end-to-end pipeline bootstrap.
5. In the variable-seed analysis, the neighborhood-graph and Leiden seeds vary together from 0–99 across 100 runs. In the fixed-seed analysis, both seeds equal 0 across 50 runs. For each reference partition, the best-match Jaccard index is the largest intersection-over-union with any reconstructed partition among the sampled cells. Standard deviations and recovery proportions describe the empirical distribution across runs; they are not confidence intervals or independent replication.
6. A mean Jaccard index of at least 0.75 is shown as a descriptive reference and is not a biological validation threshold.
7. For the resolution sweep, the full-cohort reference graph and Leiden solution use seed 0. Thirty 80% subsample graphs are reconstructed at each resolution with graph and Leiden seeds varying together from 0–29. Reconstructed partition identifiers are local to each resolution and do not map to the reference partition numbers.
8. The partitions meeting the descriptive 0.75 reference at resolutions 0.5 and 0.6 included a 993-cell and a 990-cell partition, respectively. Higher Jaccard values were therefore not restricted to the small doublet-enriched partitions. These conditional analyses nevertheless do not establish stable biological states or a transcriptional continuum.
9. All analyses use cells as computational subsamples. They add no animals and do not constitute biological replication or orthogonal validation.

### Table S3. Complete primary animal-level pseudobulk differential expression, animal-influence summaries, and detailed results for the selected seven-transcript panel.

#### General note

The seven interferon-responsive transcripts (*Irf7*, *Ifitm3*, *Isg15*, *Mx1*, *Ifit1*, *Ifit2*, and *Ifit3*) were selected during exploratory inspection of the primary dataset and used to summarize animal-level expression patterns. Primary analyses used the 7,371 Scrublet-negative microglia, with the animal/library as the biological unit. C1–C3 denote oxygen-control animals; S1–S3 denote animals receiving combined anesthesia-plus-laparotomy exposure. Parts D–E provide transcriptome-wide context for the selected-transcript results. Positive log2 fold changes denote higher expression under combined exposure, and negative values denote higher expression in oxygen controls.

#### Part A — Animal values

`TableS3_animal_values` reports each transcript’s UMI count, retained microglial UMI denominator, and direct animal-level expression, calculated as `log2[((gene UMI + 0.5)/(total retained microglial UMI + 1)) × 10^6]`. The panel value is the arithmetic mean of the seven transcript-level values, giving each transcript equal weight.

#### Part B — Full-cohort effects

`TableS3_direct_and_deseq2_effects` reports combined-exposure and oxygen-control means, the combined-exposure-minus-control mean difference, Welch–Satterthwaite 95% confidence interval, Cohen’s *d*, and two-sided exact label-permutation *P* value across all 20 three-versus-three assignments for each direct-expression outcome. These confidence intervals and exact permutation tests are descriptive and are not multiplicity-adjusted.

Transcript rows also report unshrunk DESeq2 log2 fold changes, 95% Wald intervals, Wald statistics, nominal *P* values, and Benjamini–Hochberg adjusted *P* values from animal-level pseudobulk UMI counts. Adjustment was performed across the complete 13,926-feature DESeq2 analysis, not only the seven displayed transcripts; Wald intervals are not multiplicity-adjusted. Adjusted *P* values are `NA` for five selected transcripts after DESeq2 independent filtering and were not replaced by nonsignificant numerical values. DESeq2 and direct-expression columns represent distinct estimands derived from the same six libraries.

#### Part C — Animal-influence summaries

`TableS3_leave_one_animal_out` reports the full-cohort direct-expression estimate and six systematic one-animal-omission estimates for every transcript and for the equal-weight panel. Omission rows are influence summaries rather than independent hypothesis tests; exact permutation *P* values are therefore reported only for the full three-versus-three comparison. `sign_matches_full` indicates whether the raw mean-difference sign matches the corresponding full-cohort sign. `perfect_group_separation` is a descriptive ordering indicator for the retained animal values. No animal was excluded from the full-cohort estimate.

#### Part D — Complete primary transcriptome-wide DESeq2 results

`TableS3D_primary_transcriptome_deseq2` reports all 13,926 Ensembl features tested in the primary `~ group` model. Fields include raw UMI counts for each animal/library, normalized base mean, unshrunk log2 fold change, Wald standard error and 95% interval, Wald statistic, nominal *P* value, and Benjamini–Hochberg adjusted *P* value. DESeq2 independent filtering left 6,632 nonmissing adjusted *P* values. Forty features met adjusted *P* < 0.05: 11 had positive and 29 had negative log2 fold changes. `bh_fdr_lt_0_05` is a reporting indicator rather than an additional test.

The final three columns summarize the six systematic one-animal-omission refits. They report the number of omission fits in which a feature met the false-discovery-rate threshold, the number that also retained the full-fit direction, and whether both conditions held in all six refits. These are influence summaries rather than independent replications. No animal was excluded from the full fit.

#### Part E — Transcriptome-wide animal-omission summary

`TableS3E_primary_transcriptome_loo_summary` reports, for the full fit and each one-animal omission, the group sizes, tested-feature count, number of nonmissing adjusted *P* values, number of false-discovery-rate-threshold features in each direction, and number of the 40 full-fit hits retained. Threshold-hit counts ranged from 10 to 191 across omissions. Five full-fit hits—*Ccl9*, *H1f2*, *Lst1*, *Tnfaip3*, and *Ttr*—met adjusted *P* < 0.05 with the same direction in all six omission refits.

The raw-count detection results for *Ttr* are provided in Figure S1 and Table S1. Because unfiltered droplet matrices were unavailable, biological *Ttr* expression could not be separated from ambient RNA or other sample-processing contributions. Neither the false-discovery-rate-hit counts nor their omission patterns define a coherent exposure-induced program.

### Table S4. Complete Hallmark enrichment, background-adjusted scoring sensitivity, animal-influence summaries, and fixed-score component analyses for GSE267933 microglia.

#### General note

The primary molecular analysis used 7,371 fixed reference microglia after exclusion of 90 Scrublet-predicted doublets. C1–C3 denote oxygen-control animals; S1–S3 denote animals receiving combined anesthesia-plus-laparotomy exposure. The animal/library is the biological unit. Direct-expression differences in log2 counts per million, background-adjusted `score_genes` differences, and normalized enrichment scores (NES) are distinct estimands and are not treated as mutually validating measurements.

Animal omissions assess individual-animal influence. Scoring-parameter sweeps assess estimator sensitivity, and transcript omissions assess dependence on panel definition. These analyses use the same six-animal cohort, provide no additional biological replication or orthogonal validation, and did not remove any animal from the full-cohort estimate.

#### Part A — Complete mouse MSigDB Hallmark analysis

Part A reports all 50 mouse MSigDB 2026.1.Mm Hallmark sets for the full animal-level pseudobulk fit and six systematic one-animal-omission refits. Genes were ranked by the signed DESeq2 Wald statistic. Positive NES values denote ranking toward combined exposure and negative values denote ranking toward oxygen control; NES does not measure protein abundance or pathway activity. Benjamini–Hochberg adjustment was performed across all 50 sets separately within each fit. Animal-omission rows are influence summaries rather than independent hypothesis tests.

#### Part B — Selected-panel `score_genes` parameter sensitivity

Part B reports all one-parameter-at-a-time scoring runs for the seven selected interferon-responsive transcripts. Seeds 0–19 were evaluated with `ctrl_size = 50` and `n_bins = 25`; `ctrl_size` values of 25, 50, 100, and 200 were evaluated with seed 42 and `n_bins = 25`; and `n_bins` values of 10, 25, and 50 were evaluated with seed 42 and `ctrl_size = 50`. The seed-42, `ctrl_size = 50`, `n_bins = 25` configuration appears in both the control-size and bin-number sweeps by design and is not an additional analysis. `n_expression_features` describes the symbol-deduplicated expression universe; `panel_size` and `panel_genes` describe panel membership.

Each configuration reports the full-cohort effect and the contrast after omitting S3. The latter assesses whether the sign change identified by the systematic fixed-score animal-omission analysis persisted across scoring configurations; it is not an independent experiment.

#### Part C — Selected-panel parameter-sweep summaries

Part C summarizes Part B by parameter family. It reports the number of runs and ranges or medians of full-cohort mean differences, Cohen’s *d*, exact permutation *P* values, and omit-S3 contrasts. `n_sign_preserved_omit_s3` counts runs in which the full-cohort direction was retained after S3 was omitted. These ranges describe scoring choices within the same six animals and are not confidence intervals, multiplicity-adjusted inference, or biological replication.

#### Part D — Leave-one-transcript-out panel-definition sensitivity

Part D reports seven fixed-configuration runs in which one transcript was omitted from the seven-transcript panel. `omitted_transcript` identifies the excluded transcript and is distinct from the animal-omission fields in Parts A, E, F, H, and K. Each six-transcript analysis reports the full-cohort score and the contrast after omitting S3. These rows assess dependence on panel membership rather than biological-replicate robustness.

#### Part E — Selected-panel predicted-doublet-inclusion sensitivity

Part E compares direct selected-panel mean differences from the primary 7,371-cell Scrublet-negative population with estimates from all 7,461 fixed reference microglia, including the 90 Scrublet-predicted doublets. It reports the full cohort and all six animal omissions. Fields prefixed `all_7461_microglia_` refer to all fixed reference microglia rather than all brain-cell types. The same six animals contribute under both cell-inclusion rules, so similarity is a within-cohort inclusion-sensitivity result rather than biological replication or validation.

#### Part F — Focused-Hallmark predicted-doublet-inclusion sensitivity

Part F reports paired NES values for the four focused Hallmark sets under the 7,371-cell and 7,461-cell inclusion rules for the full cohort and six animal omissions. `scrublet_negative_adjusted_p_value` and `all_7461_microglia_adjusted_p_value` were each calculated across all 50 Hallmarks within the corresponding fit, although only the four focused sets are included in this part. Agreement in NES direction or false-discovery-rate classification does not establish pathway activity, biological replication, or robustness across independent cohorts.

#### Part G — Fixed selected-panel score animal values

Part G reports the six animal-level values underlying Figure 3A. Values are animal means from `scanpy.tl.score_genes` with seed 42, `ctrl_size = 50`, and `n_bins = 25`. The score represents selected-transcript expression relative to algorithmically matched control-gene expression on the log-normalized matrix. It is distinct from direct expression in log2 counts per million and does not measure pathway activity.

#### Part H — Fixed selected-panel score effects and animal omissions

Part H reports the full-cohort effect and six systematic one-animal-omission contrasts calculated from the fixed animal values in Part G. The full three-versus-three row includes the combined-exposure-minus-oxygen-control mean difference, Welch–Satterthwaite confidence interval, Cohen’s *d*, and descriptive exact label-permutation *P* value. Exact permutation results are not assigned to unequal-group omission rows. Omission rows are influence summaries calculated from fixed scores rather than independent hypothesis tests or score refits. S3 remains in the full-cohort estimate.

#### Part I — Fixed selected-panel matched-control gene manifest

Part I reports the exact 249 algorithm-selected matched-control genes used for the fixed seven-transcript score, with matrix indices, expression bins, mean log1p-normalized expression, scoring configuration, and control-set SHA-256. The set was selected once from the 7,371-cell, 27,933-feature symbol-deduplicated matrix with seed 42, `ctrl_size = 50`, `n_bins = 25`, `use_raw = False`, and `ctrl_as_ref = True`; it was held fixed for the component and animal-omission calculations in Parts J–K. These genes are algorithmic expression-matched references rather than a prespecified negative-control panel or biological background program.

#### Part J — Fixed selected-panel score components by animal

Part J reports the selected-gene mean, fixed matched-control-gene mean, and their difference for each animal. The difference reconstructs the stored Scanpy score within numerical tolerance. Components were calculated at the cell level and averaged within animal; the cell count is reported for each animal. The decomposition is on the log-normalized score scale and is not interchangeable with the direct animal-level expression measure in Figure 2.

#### Part K — Fixed selected-panel component contrasts and animal omissions

Part K reports combined-exposure and oxygen-control means and their differences separately for the selected, fixed matched-control, and final score components in the full cohort and six animal-omission scenarios. In each row, the selected-component difference minus the matched-control-component difference reconstructs the score difference within numerical tolerance. Omission rows reuse the fixed component values and are influence summaries rather than score refits or independent tests. The arithmetic decomposition does not assign biological meaning to the matched-control genes and does not negate the separate direct-expression result.

### Table S5. GSE283401 hippocampal sample provenance, old-animal treatment contrasts, and gene-level and Hallmark time-by-treatment analyses.

#### General note

The hippocampal component of GSE283401 contains 48 bulk RNA-sequencing libraries from FACS-isolated microglia of male C57BL/6 mice. The source design included young mice aged 3–5 months and old mice aged 20–22 months. Exposed animals received 1.2% isoflurane in 30% oxygen balanced with nitrogen for 2 h plus laparotomy. Controls received carrier gas for 2 h and the same bupivacaine and meloxicam regimen, but no anesthesia or laparotomy. The 6- and 48-h endpoints were measured from the start of anesthesia or carrier-gas administration.

Formal analyses in Parts B–E were restricted to 32 old-animal libraries: eight carrier-gas controls and seven exposed animals at 6 h, and eight controls and nine exposed animals at 48 h. Each animal contributed one library and constituted one biological unit. Different animals were sampled at the two endpoints. The accession also contains 48 hemisphere libraries, which were not included in the hippocampal manifest or analyses reported here.

The 6- and 48-h counts were deposited as separate matrices. Sampling time is therefore aligned with deposited matrix and with any unrecorded processing or sequencing differences between the experiments. These analyses provide external cross-experiment context but do not constitute a within-animal time course, replication of GSE267933, or orthogonal validation. The deposited treatment label `surgery` is retained in machine-readable provenance fields; manuscript-facing fields describe the treatment as isoflurane plus laparotomy.

#### Part A — GSE283401 sample manifest

Part A reports GEO, BioSample, and SRA identifiers; source metadata; matrix order; assay information; deposited matrix; total library counts; and analysis-inclusion status for all 48 hippocampal libraries. Young animals document the hippocampal source design but were not included in the old-animal models. `included_in_old_model = TRUE` identifies the 32 libraries used in Parts B–E. `library_total` is the sum of deposited counts for that library and is descriptive rather than a normalization factor or sample-quality exclusion rule.

#### Part B — Old-animal 6-h gene-level treatment contrast

Part B reports the complete DESeq2 `~ treatment` results for the 15 old animals sampled at 6 h. The feature universe contains 26,317 Ensembl rows with at least 10 total counts across these libraries. Positive `log2FoldChange` values denote higher expression under isoflurane plus laparotomy than under carrier-gas control, and negative values denote the reverse. Log2 fold changes are unshrunk; `ci95_low` and `ci95_high` are 95% Wald confidence limits. `padj` is the Benjamini–Hochberg adjusted *P* value returned after DESeq2 independent filtering. An unavailable adjusted value is reported as `NA` and is not evidence of significance or absence of an effect.

Seventy-three Ensembl features met adjusted *P* < 0.05: 40 had positive and 33 had negative exposure-minus-control log2 fold changes. These counts describe the separately deposited 6-h experiment.

#### Part C — Old-animal 48-h gene-level treatment contrast

Part C reports the complete DESeq2 `~ treatment` results for the 17 old animals sampled at 48 h. The feature universe contains 23,732 Ensembl rows with at least 10 total counts across these libraries. Contrast orientation, unshrunk log2 fold changes, Wald confidence intervals, and adjusted-*P*-value handling follow Part B. Parts B and C are separate cross-sectional fits and are not paired measurements.

No Ensembl feature met adjusted *P* < 0.05 in this fit. The 73-versus-0 comparison between Parts B and C involves different animals, matrices, feature universes, and multiple-testing families and does not establish disappearance, emergence, or resolution over time.

#### Part D — Old-animal gene-level time-by-treatment interaction

Part D reports the complete gene-level results from `~ time + treatment + time:treatment` across all 32 old animals. The feature universe contains 30,280 shared Ensembl rows with at least 10 total counts across the combined libraries. The reported interaction log2 fold change is `(48-h exposed − 48-h control) − (6-h exposed − 6-h control)`. A positive value denotes a more positive treatment-associated contrast in the deposited 48-h experiment than in the deposited 6-h experiment; it does not denote within-animal induction or establish a temporal transition. Gene-level confidence intervals and adjusted *P* values follow Parts B–C.

Twelve Ensembl features met adjusted *P* < 0.05 in the interaction fit: five had positive and seven had negative coefficients. The signs describe between-experiment differences in treatment contrast rather than within-animal induction or suppression.

For Parts B–D, version suffixes were removed from Ensembl identifiers before mapping to mouse gene symbols with `org.Mm.eg.db` 3.22.0 and `multiVals = "first"`. `baseMean` is the DESeq2 mean of normalized counts; `lfcSE` is the standard error of the unshrunk log2 fold change; `stat` is the Wald statistic; `pvalue` is the unadjusted Wald-test *P* value; and `padj` is the Benjamini–Hochberg adjusted *P* value. Adjustment was performed separately within each model after independent filtering. Parts B, C, and D contain 12,525, 23,720, and 6,796 nonmissing adjusted *P* values, respectively.

#### Part E — Mouse MSigDB Hallmark time-by-treatment interaction

Part E reports all 50 mouse MSigDB 2026.1.Mm Hallmark gene sets. Genes from Part D were ranked by the signed interaction Wald statistic. When multiple Ensembl rows mapped to one gene symbol, the row with the largest absolute Wald statistic was retained. Enrichment was calculated with `clusterProfiler::GSEA` using `minGSSize = 10`, `maxGSSize = 500`, `pvalueCutoff = 1`, `eps = 1 × 10^−30`, and seed 42. No human-gene-set fallback was used.

`NES` is the normalized enrichment score; positive values denote ranking toward a more positive treatment contrast in the deposited 48-h experiment than in the deposited 6-h experiment. `p.adjust` is the Benjamini–Hochberg false-discovery rate across all 50 sets, not only those displayed in Figure 4. NES is a rank-based transcript-level statistic and does not measure protein abundance, transcription-factor activity, cytokine signaling, or pathway activity.

The interaction ranking contained 21,514 unique mapped mouse symbols after duplicate-symbol handling. `enrichmentScore` is the unnormalized enrichment score; `pvalue`, `p.adjust`, and `qvalue` are values returned by `clusterProfiler::GSEA`; `rank` is the ranked-list position at which the enrichment score was attained; `leading_edge` summarizes tag, list, and signal percentages; and `core_enrichment` lists leading-edge genes. `nes_rank_descending` and `absolute_nes_rank` rank the 50 sets by NES and absolute NES, respectively. `fdr_lt_0_05` indicates whether `p.adjust` is below 0.05. Eleven Hallmark sets met this threshold.

#### Data dictionary

The accompanying data dictionary describes the contents, unit or scale, and interpretation limit for Parts A–E.

### Table S6. GSE289098 same-cohort processed-count correspondence and animal-level sensitivity analyses on the fixed 7,371-cell microglial population.

#### Parts

- **Part A, global payload and analysis summary:** matrix dimensions; one-to-one cell-barcode and feature-identifier correspondence between the primary GSE267933 and alternative GSE289098 processed-count matrices; matrix-wide and fixed-analysis-population count differences; common-estimand checks; and fixed feature universes.
- **Part B, library mapping and counts:** animal/library identifiers, GSE267933 GSM/BioSample/SRA Experiment accessions, GSE289098 barcode-suffix mapping, all-cell and microglial cell counts, and per-library count differences.
- **Part C, paired animal values:** transcript UMI counts, all-feature UMI denominators, transcript-level log2 counts-per-million values, equal-weight seven-transcript values, and GSE289098-minus-primary paired differences for each animal.
- **Part D, paired full-cohort direct effects:** group means, combined-exposure-minus-oxygen-control differences, 95% Welch–Satterthwaite confidence intervals, Cohen’s *d*, and exact permutation *P* values for the seven transcripts and their equal-weight panel.
- **Part E, panel full-cohort and one-animal-omission results:** matched full-cohort and six systematic animal-omission panel estimates under both count matrices, including paired contrast differences and direction concordance.
- **Part F, selected-transcript DESeq2 results:** unshrunk DESeq2 estimates under both count matrices using the same primary-derived 13,926-feature universe.
- **Part G, DESeq2 size factors:** median-ratio size factors for each animal and count matrix.
- **Part H, transcript-level full-cohort and one-animal-omission results:** matched direct-expression effects for all seven transcripts. Forty-eight of 49 full-cohort or omission scenarios retained the same direction. The only boundary crossing was the near-zero *Isg15* estimate after S3 was omitted (primary, −0.031; GSE289098, +0.052 log2 counts per million).
- **Data dictionary:** row unit and purpose of each part.

#### General notes

1. GSE289098 contains a Cell Ranger v3.0.2-aggregated processed-count matrix derived from the same six GSE267933 libraries and 20,684 filtered cell barcodes. It adds no animals and does not provide independent biological replication or orthogonal validation. Machine-readable `source_group` and `dropped_source_group` fields retain the deposited labels `Control` and `Surgery`; manuscript-facing text maps these labels to oxygen control and combined sevoflurane-plus-laparotomy exposure, respectively.
2. The analysis carries forward the exact 7,371 Scrublet-negative microglia by sample-aware barcode matching. No GSE289098-specific quality-control threshold, microglial reclassification, clustering, or subtype definition was applied. Cells remain nested within animals; the biological unit is the animal/library (*n* = 3 per group).
3. Primary counts were obtained from the verified integer-count object corresponding to the six original GSE267933 count matrices. Processed-object `.raw` attributes containing log-normalized expression were not used as counts. GSE289098 counts were obtained from the matched 27,998-feature processed matrix.
4. The seven transcripts selected during exploratory inspection were *Irf7*, *Ifitm3*, *Isg15*, *Mx1*, *Ifit1*, *Ifit2*, and *Ifit3*. For each transcript and animal, expression was calculated as `log2[((transcript UMI + 0.5)/(total UMI across all 27,998 features + 1)) × 10^6]`. The panel value is the arithmetic mean of the seven transcript-level values rather than a value calculated from summed panel counts.
5. Direct effects are combined-exposure-minus-oxygen-control differences in animal-level values. Confidence intervals use Welch–Satterthwaite degrees of freedom, and Cohen’s *d* uses the pooled within-group standard deviation. Descriptive full-cohort exact two-sided permutation *P* values enumerate all 20 possible three-versus-three label assignments without a +1 correction. No permutation *P* values are assigned to unequal-group omission rows.
6. Alternative-minus-primary quantities are deterministic paired descriptions of two processing outputs from the same libraries and were not assigned paired-test *P* values.
7. DESeq2 models used `~ group`, with combined exposure relative to oxygen control. The primary total-count threshold of at least 10 defines the common 13,926-feature universe for both fits; 13,831 features would independently meet that threshold in GSE289098. Wald confidence intervals equal log2 fold change ± 1.96 standard errors. Benjamini–Hochberg adjusted *P* values were calculated separately within each matrix after DESeq2 independent filtering. Because independent filtering is matrix-specific, adjusted *P* values are reported for completeness but are not used as an equivalence or concordance criterion.
8. Full-cohort and all six one-animal-omission panel estimates remained positive under both count matrices. The estimates changed little in this single same-cohort processed-count comparison, which does not establish biological replication or generalizability across cohorts.

