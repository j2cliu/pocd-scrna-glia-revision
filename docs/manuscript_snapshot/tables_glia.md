# Tables

## Table 1. Per-animal Scrublet exclusion and UMI summary for the GSE267933 reference microglial population

| Animal/library | Exposure group | Cells before exclusion, *n* | Predicted doublets removed, *n* (%) | Cells retained, *n* | Microglial UMIs before exclusion, *n* | Microglial UMIs removed, *n* (%) | Microglial UMIs retained, *n* | Median UMI per retained cell |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| C1 | Oxygen control | 1,147 | 13 (1.13) | 1,134 | 8,180,581 | 152,125 (1.86) | 8,028,456 | 7,024.5 |
| C2 | Oxygen control | 932 | 15 (1.61) | 917 | 6,335,382 | 160,820 (2.54) | 6,174,562 | 6,582.0 |
| C3 | Oxygen control | 1,282 | 10 (0.78) | 1,272 | 9,846,616 | 140,412 (1.43) | 9,706,204 | 7,518.5 |
| S1 | Combined exposure | 1,321 | 16 (1.21) | 1,305 | 10,344,043 | 247,871 (2.40) | 10,096,172 | 7,653.0 |
| S2 | Combined exposure | 1,622 | 21 (1.29) | 1,601 | 12,247,978 | 244,331 (1.99) | 12,003,647 | 7,399.0 |
| S3 | Combined exposure | 1,157 | 15 (1.30) | 1,142 | 11,476,010 | 223,109 (1.94) | 11,252,901 | 9,793.5 |
| **Total** | **All animals** | **7,461** | **90 (1.21)** | **7,371** | **58,430,610** | **1,168,668 (2.00)** | **57,261,942** | — |

*Note.* Values refer to the fixed reference microglial population, not all barcodes in each sequencing library. Scrublet was applied separately to each animal/library; predicted doublets are model-based calls and were not confirmed experimentally. Percentages removed were calculated relative to the corresponding cell or UMI total before exclusion. The retained cells constitute the 7,371-cell label-independent molecular-analysis population. C1–C3 were oxygen controls, whereas S1–S3 received combined exposure to 2.5% sevoflurane in 50% O₂ for 30 min plus laparotomy. No anesthesia-only group was available. UMI, unique molecular identifier.

## Table 2. Per-animal composition of the six non-Rare reference partitions and conditional exposure contrasts (GSE267933)

### Part A. Per-animal composition

| Reference partition (legacy label) | C1 | C2 | C3 | S1 | S2 | S3 |
|---|---:|---:|---:|---:|---:|---:|
| 0 (Inflammatory) | 176 (15.5) | 323 (35.3) | 210 (16.5) | 254 (19.5) | 208 (13.0) | 479 (41.8) |
| 1 (Transitional-A) | 205 (18.1) | 102 (11.1) | 192 (15.1) | 287 (22.1) | 455 (28.5) | 140 (12.2) |
| 2 (Homeostatic-A) | 206 (18.2) | 207 (22.6) | 228 (17.9) | 263 (20.2) | 99 (6.2) | 282 (24.6) |
| 3 (Homeostatic-B) | 206 (18.2) | 98 (10.7) | 258 (20.3) | 257 (19.8) | 275 (17.2) | 130 (11.3) |
| 4 (Homeostatic-C) | 252 (22.2) | 92 (10.1) | 249 (19.6) | 167 (12.8) | 432 (27.0) | 17 (1.5) |
| 5 (Transitional-B) | 89 (7.8) | 93 (10.2) | 134 (10.5) | 73 (5.6) | 130 (8.1) | 99 (8.6) |
| **Non-Rare total** | **1,134 (100.0)** | **915 (100.0)** | **1,271 (100.0)** | **1,301 (100.0)** | **1,599 (100.0)** | **1,147 (100.0)** |

### Part B. Animal-level group contrasts

| Reference partition (legacy label) | Oxygen control, mean % | Combined exposure, mean % | Difference (95% CI), percentage points | Cohen’s *d* | Exact *P* |
|---|---:|---:|---:|---:|---:|
| 0 (Inflammatory) | 22.4 | 24.8 | 2.3 (−28.8 to 33.4) | 0.17 | 0.70 |
| 1 (Transitional-A) | 14.8 | 20.9 | 6.1 (−11.3 to 23.6) | 0.97 | 0.30 |
| 2 (Homeostatic-A) | 19.6 | 17.0 | −2.6 (−24.5 to 19.3) | −0.37 | 0.90 |
| 3 (Homeostatic-B) | 16.4 | 16.1 | −0.3 (−11.0 to 10.4) | −0.06 | 0.90 |
| 4 (Homeostatic-C) | 17.3 | 13.8 | −3.5 (−30.1 to 23.1) | −0.35 | 0.70 |
| 5 (Transitional-B) | 9.5 | 7.5 | −2.1 (−5.6 to 1.4) | −1.34 | 0.30 |

*Note.* Entries in Part A are cell count (percentage). Percentages were calculated within each animal after excluding reference partition 6 (Rare; 94 cells) from the fixed 7,461-cell population, yielding 7,367 cells. Scrublet status was not used to define this composition population, which therefore differs from the 7,371-cell molecular-analysis population in Table 1. Legacy labels in parentheses identify the original annotation and are not treated as validated biological states. Part B uses the animal/library as the biological unit. Differences are the unweighted combined-exposure-minus-oxygen-control mean differences. Confidence intervals are Welch–Satterthwaite 95% intervals; Cohen’s *d* is a point estimate calculated using the pooled within-group standard deviation. Exact two-sided *P* values were obtained from all 20 assignments of three of the six animals to the combined-exposure group. These unadjusted tests are descriptive; the six partition percentages are compositional and not independent. CI, confidence interval.
