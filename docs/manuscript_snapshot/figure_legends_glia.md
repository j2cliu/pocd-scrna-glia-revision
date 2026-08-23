# Figure legends

## Figure 1. Study design and audit of the submitted microglial partition

**A**, GSE267933 comprised hippocampal single-cell RNA sequencing from
18-month-old male mice assigned to oxygen control (C1–C3) or combined
2.5% sevoflurane-plus-laparotomy exposure (S1–S3), with tissue collected 24 h
after exposure. The frozen submitted microglial set contained 7,461 cells.
Panels B–C audit that submitted set; Panel D excludes the 94-cell
Rare/partition-6 label (7,367 cells), whereas subsequent primary-cohort,
label-independent molecular analyses exclude 90 Scrublet-predicted doublets
(7,371 cells).
Because no anesthesia-only group was available, anesthesia and laparotomy
effects cannot be separated.

**B**, UMAP of all 7,461 frozen microglia colored by the seven submitted numeric
partitions. For traceability to the submitted manuscript, partitions 0–6
correspond to Inflammatory, Transitional-A, Homeostatic-A, Homeostatic-B,
Homeostatic-C, Transitional-B, and Rare, respectively; these names are not
treated as validated biological states.

**C**, adjusted Rand index (ARI) between the submitted partition and Leiden
partitions obtained with seeds 0, 1, 7, 42, 123, and 2024 on the same stored
neighborhood graph at resolution 0.4. The number of recovered partitions
(`k`) is shown for each seed. This is an algorithmic seed-sensitivity
diagnostic, not bootstrap resampling or biological replication.

**D**, per-animal percentage of each non-Rare submitted partition, calculated
relative to that animal’s total non-Rare microglia. The forest plot shows the
combined-exposure-minus-control mean percentage-point difference with
Welch–Satterthwaite 95% confidence intervals and two-sided exact
label-permutation *P* values across all 20 allocations of three of the six
animals to the exposed group. The six tests are unadjusted and descriptive,
and the six partition percentages are compositional. No contrast showed
perfect animal-level separation. These summaries are conditional on the
seed-sensitive submitted partition and are not interpreted as reproducible
state expansion.

## Figure 2. Animal-level transcriptome-wide results and expression of selected interferon-responsive transcripts

**A**, MA plot of all 13,926 features tested by animal-level pseudobulk DESeq2.
The vertical axis shows the unshrunk log2 fold change for combined
anesthesia-plus-laparotomy exposure relative to oxygen control. Gray points
denote features that did not meet Benjamini–Hochberg false discovery rate
(FDR) < 0.05. Blue and orange points denote the 29 control-directed and 11
exposure-directed features, respectively, that met this threshold. Open
diamonds identify the seven transcripts summarized in Panels B–D; none met
FDR < 0.05. Open circles identify the five full-cohort hits that retained the
FDR threshold and direction in all six one-animal-omission refits. The large
control-directed *Ttr* estimate should be interpreted cautiously because
unfiltered droplet matrices were unavailable and biological expression could
not be distinguished from ambient-RNA or other sample-processing
contributions.

**B**, animal-level expression of *Irf7*, *Ifitm3*, *Isg15*, *Mx1*, *Ifit1*,
*Ifit2*, and *Ifit3* among the 7,371 Scrublet-negative microglia. For
transcript *g* in mouse *j*, expression was calculated as
log2[((*Y*<sub>*jg*</sub> + 0.5)/(*N*<sub>*j*</sub> + 1)) × 10<sup>6</sup>],
where *Y*<sub>*jg*</sub> was the summed UMI count for that transcript and
*N*<sub>*j*</sub> was the summed UMI count across all features from that
mouse. Points show individual mice, and short colored ticks show group means.
Within each transcript row, points follow the fixed vertical order shown by
the C1–C3 and S1–S3 key.

**C**, each mouse’s equal-weight seven-transcript value, calculated as the
arithmetic mean of the seven transcript-level values. Horizontal lines show
group means. The combined-exposure-minus-control mean difference was +1.069
log2 counts per million (95% Welch–Satterthwaite confidence interval, −1.638
to +3.775; Cohen’s *d* = 1.325; descriptive two-sided exact permutation
*P* = 0.20).

**D**, the full-cohort estimate and six systematic one-animal-omission
estimates of the same mean difference. These omissions were influence
diagnostics, not independent hypothesis tests, and no mouse was excluded from
the full-cohort analysis. All omission estimates remained positive. Omitting
S3 reduced the difference to +0.449 log2 counts per million, approximately
58% below the full-cohort estimate. C1–C3 denote oxygen-control mice; S1–S3
are source identifiers for mice in the combined anesthesia-plus-laparotomy
group.

## Figure 3. Animal influence and estimator choice alter interferon-related summaries

All panels use the six animals/libraries as the biological units. C1–C3 denote
oxygen-control animals; S1–S3 are source identifiers for animals receiving the
combined anesthesia-plus-laparotomy exposure. Unless otherwise stated, the
primary analysis set contains 7,371 frozen submitted microglia after exclusion
of 90 Scrublet-predicted doublets. Direct-expression differences in log2
counts per million,
background-adjusted `score_genes` differences, and normalized enrichment
scores (NES) are distinct estimands and are displayed on separate axes. Blue
circles denote oxygen-control animals and orange triangles denote
combined-exposure animals.

**A**, fixed-configuration background-adjusted scores for the seven
interferon-responsive transcripts selected during exploratory inspection
(*Irf7*, *Ifitm3*, *Isg15*, *Mx1*, *Ifit1*, *Ifit2*, and *Ifit3*).
`scanpy.tl.score_genes` was applied to the
symbol-deduplicated, log-normalized matrix with seed 42, `ctrl_size = 50`,
`n_bins = 25`, `use_raw = False`, and `ctrl_as_ref = True`; cell-level scores
represent the selected-gene mean relative to matched control-gene expression
and were then averaged within animal. Panel A contains two aligned displays.
The upper display shows all six animal values; short black horizontal ticks
show group means. The lower display shows the full-cohort estimate and six
systematic one-animal-omission estimates. The diamond and horizontal interval
show the full combined-exposure-minus-control difference of +0.0465 (95%
Welch–Satterthwaite confidence interval, −0.1766 to +0.2696; Cohen’s *d* =
0.714; descriptive two-sided exact label-permutation *P* = 0.80). Omission
rows are influence diagnostics rather than independent tests. The estimate
was −0.00609 after omitting S3; S3 was retained in the full-cohort estimate
and was not designated as a technical outlier. The omit-S3 sign change also
occurred in every tested seed, control-set-size, bin-number, and
leave-one-transcript score analysis (Table S4); these assess estimator or
panel-definition sensitivity rather than biological-replicate robustness.
Table S4I–K freezes the exact 249-gene matched-control set and reports the
selected-gene and matched-control components separately. The matched-control
component is an algorithmic reference term, not a validated neutral biological
background program.

**B**, mouse MSigDB 2026.1.Mm Hallmark enrichment after genes were
ranked by the signed DESeq2 Wald statistic. Positive NES values indicate that
genes in the set ranked toward combined exposure, and negative values indicate
ranking toward oxygen control; NES does not measure protein abundance or
pathway activity. Values are shown for the full fit and all six systematic
leave-one-animal-out refits; these are influence diagnostics rather than
independent hypothesis tests. Benjamini–Hochberg false discovery rates were
calculated across all 50 Hallmark sets within each fit, not across only the four
sets displayed here; complete results are provided in Table S4. After S3 was
omitted, IFN-alpha-response NES changed from +2.171 to −1.098 (false discovery
rate from 3.16 × 10^−6 to 0.344), and IFN-gamma-response NES changed from
+1.580 to −1.569 (false discovery rate from 0.00469 to 0.00334). Thus, both
focused interferon Hallmarks reversed direction when S3 was omitted;
IFN-alpha response no longer met the 0.05 false-discovery-rate threshold,
whereas IFN-gamma response remained below that threshold in the opposite
direction. The TNFA signaling via NF-κB and inflammatory-response Hallmark
gene sets ranked toward controls in the full fit and every omission fit. These
were focused secondary results, not unique transcriptome-wide findings: 43 of 50
Hallmarks had negative NES in the full fit, and 27 of 50 remained negative in
the full fit and all six omission refits. The results do not establish lower
TNF abundance, NF-κB activity, or inflammatory-pathway activity.

**C**, sensitivity to retaining the 90 Scrublet-predicted doublets. Filled
points denote estimates after excluding the predicted doublets (7,371 cells),
and open points denote estimates after retaining them (7,461 cells). The same
symbol key applies to both subpanels. The left subpanel shows the
combined-exposure-minus-control difference for the selected
seven-transcript expression summary. This estimate was +1.069 log2 counts per
million after exclusion and +1.066 after retention. The right subpanel shows
the full-cohort NES values for the four focused Hallmark gene sets under the
same two cell-inclusion rules. Across these four Hallmarks and the seven full
or one-animal-omission models, all 28 paired comparisons retained the same NES
direction and FDR-below-0.05 classification; the largest absolute paired NES
difference was 0.0753 (Table S4). This same-cohort sensitivity analysis did
not provide independent biological replication or orthogonal validation.

## Figure 4. An independent bulk-microglial dataset provides bounded 6-h/48-h context for treatment-associated transcript contrasts

**A**, comparison of the primary GSE267933 cohort and the external
GSE283401 dataset. GSE267933 comprised 18-month-old male mice sampled 24 h after
combined anesthesia-plus-laparotomy exposure comprising 2.5% sevoflurane for
30 min plus laparotomy, or assignment
to the oxygen-control group (*n* = 3 animals per group); the primary analysis
used hippocampal microglial single-cell RNA sequencing. GSE283401 used male
C57BL/6 mice, FACS-isolated microglia, and bulk RNA sequencing. The analyzed
hippocampal component compared 1.2% isoflurane in 30% oxygen balanced with
nitrogen for 2 h plus laparotomy with carrier-gas control. Both groups received
the same bupivacaine and meloxicam regimen. The formal analyses shown here were restricted
to the 20–22-month-old subset: 8 carrier-gas controls and 7 exposed animals at
6 h, and 8 carrier-gas controls and 9 exposed animals at
48 h. These endpoints were measured from the start of anesthesia or carrier-gas
administration. Within these 48 hippocampal libraries, the source design also contained
young animals aged 3–5 months (*n* = 4 per treatment group at each time);
these samples are documented in
Table S5 but were not included in the formal models. The biological unit was
the animal/library. Different animals were sampled at 6 h and 48 h, and the
external dataset also differed from the primary cohort in assay, anesthetic,
exposure duration, and experimental cohort.

**B**, separately fitted old-animal treatment contrasts for the same seven
interferon-responsive transcripts displayed for the primary cohort (*Irf7*,
*Ifitm3*, *Isg15*, *Mx1*, *Ifit1*, *Ifit2*, and *Ifit3*).
DESeq2 models of the form `~ treatment` were fitted independently at 6 h and
48 h. Points show unshrunk log2 fold changes and horizontal lines show 95%
Wald confidence intervals. Positive values denote higher expression after
isoflurane plus laparotomy than in the carrier-gas control group;
negative values denote the reverse. Orange and blue point fills indicate
positive and negative estimates, respectively, not statistical significance.
Gene-level Benjamini–Hochberg adjusted *P* values were calculated separately
within each time-specific fit after
DESeq2 independent filtering. None of the seven displayed transcripts had an
adjusted *P* value below 0.05 at either time. The two panels describe separate
cross-sectional experiments and are not paired or longitudinal comparisons.
Across the complete time-specific fits reported in Table S5, 73 Ensembl
features met adjusted *P* < 0.05 at 6 h (40 positive and 33 negative
estimates), whereas none met this threshold at 48 h; these counts arise from
different animals, matrices, feature universes, and multiple-testing families
and do not establish disappearance, emergence, or resolution over time.

**C**, focused mouse MSigDB Hallmark results from the formal old-animal
DESeq2 model `~ time + treatment + time:treatment` across all 32 libraries.
The complete gene-level interaction contained 12 features with adjusted
*P* < 0.05 (five positive and seven negative coefficients); these directions
do not denote within-animal induction or suppression.
Genes were ranked by the signed Wald statistic for the interaction
(48-h exposed-minus-control contrast) minus (6-h exposed-minus-control
contrast), and enrichment was evaluated with the 50 mouse MSigDB
2026.1.Mm Hallmark gene sets. The four manuscript-focused sets are displayed;
complete results for all 50 sets are provided in Table S5. Positive normalized
enrichment scores (NES) indicate that genes in a set ranked toward a more
positive treatment contrast
in the deposited 48-h experiment than in the deposited 6-h experiment; NES
does not measure protein abundance or pathway activity. Filled points denote
Benjamini–Hochberg false-discovery rates below 0.05 after adjustment across
all 50 Hallmark sets, and the open point denotes a false-discovery rate of at
least 0.05. The displayed NES (false-discovery rate) values were +2.624
(7.31 × 10^−11) for interferon-alpha response, +2.365
(2.36 × 10^−11) for interferon-gamma response, +1.343 (0.114) for TNFA
signaling via NF-κB, and +1.709 (0.00129) for inflammatory response. Because
different animals were used at the two times and time is inseparable from the
two deposited count matrices and any unrecorded time-aligned processing or
sequencing effects, this interaction is a between-experiment difference in
treatment-associated contrast. It does not establish within-animal change, a
temporal trajectory, independent replication of GSE267933, orthogonal
validation, or causal attribution to sampling time.

## Figure S1. Retrospective annotation and contamination audit of the submitted seven-way microglial partition.

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

## Figure S2. Seed sensitivity of the submitted microglial partition and conditional graph-reconstruction diagnostics.

**A**, Pairwise adjusted Rand indices (ARIs) among Leiden partitions obtained with seeds 0, 1, 7, 42, 123, and 2024 at resolution 0.4 on the exact neighborhood graph stored in the submitted 7,461-cell microglial object. The right-hand columns give the ARI of each run against the submitted numeric partition and the number of partitions recovered (*k*). Seed 0 reproduced the submitted labels exactly (ARI = 1.000); the other five seeds recovered six or seven partitions with ARIs of 0.403–0.586 against the submitted labels. Across the 15 unique seed pairs, ARI ranged from 0.365 to 0.714 (mean, 0.529). **B**, Conditional 80% cell-subsampling diagnostic at resolution 0.4. Each run drew 5,969 of the same 7,461 cells, retained the inherited whole-cell PC1–PC20 coordinates, and rebuilt a 15-neighbor graph. Filled purple points and horizontal bars show the mean ± SD of each submitted partition's best-match Jaccard index across 100 runs in which the graph and Leiden seeds varied together from 0 to 99; open points show the mean across 50 runs in which both seeds were fixed at 0. The SD bars describe run-to-run dispersion and are not confidence intervals. Numeric partitions are primary; names in parentheses are retrospective traceability labels. Partition 6/Rare is doublet-enriched and is not interpreted as a biological state. The dashed 0.75 line is a descriptive reference, not a validation threshold. Lower panel, distribution of the number of partitions recovered by the same two subsampling arms; six or seven partitions were recovered in 92% of variable-seed runs and 88% of fixed-seed runs. **C**, Separate conditional resolution diagnostic. A full-cohort reference graph was rebuilt from the inherited PC coordinates with graph and Leiden seed 0 at each resolution from 0.1 to 1.0; 30 80% subsample graphs were then rebuilt per resolution with graph and Leiden seeds varying together from 0 to 29. Purple and open points show the across-partition mean and median best-match Jaccard index; gray vertical ranges show the minimum and maximum per-partition means and are not uncertainty intervals. Values above each resolution give the number of reconstructed partitions with mean Jaccard ≥ 0.75 divided by the full-cohort partition count; the second line of each x-axis label gives the full-cohort partition count/modal subsample count. At resolutions 0.5 and 0.6, two reconstructed partitions met the 0.75 reference; Table S2F shows that these included 993- and 990-cell partitions, respectively, so the pattern is not confined to a small partition. Panels B–C are conditional diagnostics after graph reconstruction from a fixed inherited coordinate system. They are not tests on the submitted stored graph, do not resample animals or reconstruct the full analysis pipeline, and do not validate discrete biological states or prove a transcriptional continuum.

## Figure S3. Same-cohort GSE289098 count reprocessing yields similar animal-level estimates for the selected seven-transcript panel.

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
