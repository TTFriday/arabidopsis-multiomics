# Data profile


## emqtl  (`emqtl.csv`)

- rows: **265,776**  cols: 5  (~38.2 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| Gene | str | 0.0 | 2,645 | AT5G32460 (12143); AT5G34810 (11240) |
| Methylation | str | 0.0 | 109,260 | 5_22601552_CG (55); 3_13587694_CG (34) |
| Pearson | float64 | 0.0 | 265,655 | min=-0.8108 med=0.3243 max=0.8836 |
| Pvalue | float64 | 0.0 | 34,195 | min=8.88e-206 med=2.56e-20 max=2.32e-14 |
| FDR | float64 | 0.0 | 34,990 | min=7.3e-200 med=2.3e-17 max=1.87e-08 |

## epathway  (`epathway.csv`)

- rows: **122,344**  cols: 8  (~20.6 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| Pathway | str | 0.0 | 375 | PWY-5267 (12291); PWYQT-4471 (6005) |
| Methylation | str | 0.0 | 31,637 | 2_14402954_CHG (78); 2_14402737_CG (71) |
| Pearson | float64 | 0.0 | 111,436 | min=-0.4364 med=-0.2026 max=0.4116 |
| Pvalue_Pearson | float64 | 0.0 | 6,419 | min=3.27e-30 med=2.64e-08 max=5.15e-07 |
| FDR_Pearson | float64 | 0.0 | 24,367 | min=5.64e-24 med=3.02e-05 max=0.04968 |
| Spearman | float64 | 0.0 | 111,452 | min=-0.4065 med=-0.2045 max=0.4139 |
| Pvalue_Spearman | float64 | 0.0 | 8,206 | min=4.63e-27 med=1.75e-08 max=5.15e-07 |
| FDR_Spearman | float64 | 0.0 | 26,623 | min=8e-21 med=2.36e-05 max=0.04898 |

## gwas  (`gwas.csv`)

- rows: **89,884**  cols: 14  (~51.19 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| SNP | str | 0.0 | 49,974 | _ (77); chr1_12760855 (38) |
| Chr | str | 0.0857 | 5 | chr5 (29339); chr1 (21676) |
| Position | float64 | 0.0857 | 49,950 | min=570 med=1.239e+07 max=3.041e+07 |
| Score | float64 | 0.0 | 30,696 | min=6.274 med=9.976 max=119.6 |
| MAF | float64 | 0.0 | 1,181 | min=0.01029 med=0.02196 max=0.5 |
| Study | int64 | 0.0 | 254 | min=4 med=587 max=768 |
| Phenotype | str | 0.0 | 248 | LRDpBZ75 (11527); clim-gs10 (9666) |
| Phenotype ontology | str | 1.4029 | 44 | other miscellaneous trait (60985); root morphology trait (13502) |
| Reference | str | 0.0857 | 4 | G (23402); C (22981) |
| Alteration | str | 0.0857 | 4 | A (28290); T (27312) |
| Genename | str | 58.5499 | 8,376 | AT3G20910 (235); AT4G02930 (192) |
| Function | str | 80.4492 | 3 | MISSENSE (10103); SILENT (7247) |
| Impact | str | 8.2462 | 4 | MODIFIER (64899); MODERATE (10063) |
| Effect | str | 8.2462 | 12 | INTERGENIC (45215); INTRON (14262) |

## gwaseqtl  (`gwaseqtl.csv`)

- rows: **286,032**  cols: 10  (~149.04 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| eQTL loci | str | 0.0 | 2,661 | chr1_14119557 (1784); chr3_446110 (1745) |
| eQTL Allels | str | 0.0 | 74 | G_A (62462); C_T (59640) |
| eQTL Alteration | str | 0.0 | 316 | 42 (19352); 43 (16123) |
| GWAS loci | str | 0.0 | 6,968 | chr1_15650600 (1747); chr1_14056466 (1627) |
| GWAS Allels | str | 0.0 | 4 | C (81834); G (76955) |
| GWAS Alteration | str | 0.0 | 4 | T (95916); A (86503) |
| LD | float64 | 0.0 | 22,501 | min=0.5 med=0.5504 max=1 |
| GWAS phenotype | str | 0.0 | 94 | clim-gs10 (77184); clim-tmin2 (45042) |
| Type | str | 0.0 | 3 | trans (207996); cis (61367) |
| Phenotype ontology | str | 0.0 | 18 | other miscellaneous trait (258942); leaf chlorosis (8858) |

## phenotypepathway  (`phenotypepathway.csv`)

- rows: **1,354**  cols: 10  (~0.39 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| Phenotype | str | 0.0 | 88 | Se82 (182); Co59 (179) |
| Phenotype ontology | str | 1.6987 | 22 | other miscellaneous trait (322); selenium concentration (182) |
| Pathway | str | 0.0 | 403 | PWYQT-4450 (42); PWY-7618 (39) |
| Pearson | float64 | 0.0 | 1,275 | min=-0.6571 med=-0.3047 max=0.7522 |
| Pvalue_Pearson | float64 | 0.0 | 1,233 | min=1.21e-72 med=1.77e-17 max=0.001215 |
| FDR_Pearson | float64 | 0.0 | 1,087 | min=7.1e-70 med=1.195e-16 max=0.04845 |
| Spearman | float64 | 0.0 | 1,275 | min=-0.6892 med=-0.3137 max=0.7516 |
| Pvalue_Spearman | float64 | 0.0 | 1,239 | min=7.24e-82 med=1.03e-24 max=0.002492 |
| FDR_Spearman | float64 | 0.0 | 1,145 | min=4.24e-79 med=4.57e-24 max=0.04899 |
| Study | str | 0.0 | 7 | study16 (969); study34 (327) |

## twas  (`twas.csv`)

- rows: **62,754**  cols: 10  (~17.95 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| Phenotype | str | 0.0 | 145 | Na23 (10689); Co59 (9648) |
| Phenotype ontology | str | 2.7249 | 29 | sodium concentration (10689); cobalt concentration (9648) |
| Gene | str | 0.0 | 16,897 | AT3G21690 (45); AT4G23670 (43) |
| Pearson | float64 | 0.0 | 62,701 | min=-0.7678 med=-0.305 max=0.8618 |
| Pvalue_Pearson | float64 | 0.0 | 29,477 | min=2.67e-160 med=7.755e-21 max=0.002299 |
| FDR_Pearson | float64 | 0.0 | 27,582 | min=6.45e-156 med=8.06e-20 max=0.04999 |
| Spearman | float64 | 0.0 | 62,704 | min=-0.8205 med=-0.3148 max=0.8519 |
| Pvalue_Spearman | float64 | 0.0 | 34,381 | min=2.79e-156 med=6.345e-29 max=0.003225 |
| FDR_Spearman | float64 | 0.0 | 32,844 | min=6.75e-152 med=3.97e-28 max=0.04993 |
| Study | str | 0.0 | 11 | study16 (53237); study34 (5252) |

## ewas  (`ewas.csv`)

- rows: **3,993,378**  cols: 10  (~1153.17 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| Phenotype | str | 0.0 | 114 | Na23 (771710); Co59 (705051) |
| Phenotype ontology | str | 0.1804 | 38 | sodium concentration (771710); cobalt concentration (705051) |
| Methylation | str | 0.0 | 1,222,388 | 5_7532747_CG (31); 5_8612498_CG (29) |
| Pearson | float64 | 0.0 | 3,993,313 | min=-1 med=-0.3042 max=0.9996 |
| Pvalue_Pearson | float64 | 0.0 | 3,838,582 | min=8.389e-161 med=1.214e-24 max=0.001046 |
| FDR_Pearson | float64 | 0.0 | 3,662,110 | min=1.448e-154 med=1.479e-23 max=0.05 |
| Spearman | float64 | 0.0 | 3,992,697 | min=-1 med=-0.3113 max=0.9883 |
| Pvalue_Spearman | float64 | 0.0 | 3,992,631 | min=0 med=7.521e-29 max=0.001507 |
| FDR_Spearman | float64 | 0.0 | 3,803,781 | min=0 med=6.519e-28 max=0.05 |
| Study | str | 0.0 | 14 | study16 (3803909); study34 (138965) |

## ciseqtl  (`ciseqtl.xls`)

- rows: **11,796**  cols: 8  (~3.02 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| SNP | str | 0.0 | 10,509 | chr4_16921004 (9); chr4_11443907 (9) |
| Alleles | str | 0.0 | 452 | G_A (1622); C_T (1488) |
| Alteration number | str | 0.0 | 1,248 | 42 (183); 45 (178) |
| Gene | str | 0.0 | 2,643 | AT5G65687 (338); AT5G25130 (314) |
| Beta | float64 | 0.0 | 11,470 | min=-1.438 med=0.4685 max=1.909 |
| t-stat | float64 | 0.0 | 11,471 | min=-7.841 med=5.183 max=10.85 |
| p-value | float64 | 0.0 | 2,822 | min=3.66e-25 med=2.43e-07 max=9.7e-07 |
| FDR | float64 | 0.0 | 6,373 | min=2.22e-16 med=0.02502 max=0.04997 |

## transeqtl  (`transeqtl.xls`)

- rows: **10,118**  cols: 8  (~2.59 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| SNP | str | 0.0 | 6,776 | chr1_11301651 (101); chr1_11301693 (101) |
| Alleles | str | 0.0 | 681 | G_A (1642); C_T (1540) |
| Alteration number | str | 0.0 | 422 | 9 (1164); 8 (1015) |
| Gene | str | 0.0 | 580 | AT3G32092 (1421); AT1G13310 (1386) |
| Beta | float64 | 0.0 | 8,593 | min=-1.171 med=0.5397 max=2.315 |
| t-stat | float64 | 0.0 | 8,593 | min=-12.35 med=7.732 max=18.1 |
| p-value | float64 | 0.0 | 3,850 | min=1.06e-58 med=3.47e-14 max=1.18e-12 |
| FDR | float64 | 0.0 | 3,693 | min=4.51e-48 med=0.00147 max=0.0499 |

## environmenteqtl  (`environmenteqtl.xlsx`)

- rows: **68,403**  cols: 8  (~20.28 MB in memory)

| column | dtype | missing % | unique | notes |
|---|---|---|---|---|
| eQTL | str | 0.0 | 11,593 | chr5_25983733 (19); chr5_26030954 (19) |
| Alleles | str | 0.0 | 468 | G_A (9812); C_T (9509) |
| Alteration | str | 0.0 | 1,140 | 40 (1499); 43 (1320) |
| Environment | str | 0.0 | 22 | bio3 (6191); Longitude (5795) |
| Beta | float64 | 0.0 | 66,466 | min=-2.075 med=0.2072 max=2.092 |
| p-value | float64 | 0.0 | 49,998 | min=5.29e-20 med=0.0006636 max=0.009768 |
| FDR | float64 | 0.0 | 45,348 | min=2.4e-14 med=0.00664 max=0.04999 |
| type | str | 0.0 | 3 | cis (44300); trans (20264) |