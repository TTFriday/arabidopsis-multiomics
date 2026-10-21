# Cross-layer integration audit


## Methylation feature space

| set | n sites |
|---|---|
| ewas_unique | 1,222,388 |
| emqtl_unique | 109,260 |
| epathway_unique | 31,637 |
| ewas_and_emqtl | 91,258 |
| ewas_and_epathway | 27,347 |
| emqtl_and_epathway | 7,728 |
| all_three | 7,247 |
| jaccard_ewas_epathway | 0.02229354402703888 |

### methylation site chrom/context distribution

- parseable as `chrom_pos_CONTEXT`: 1,244,199 / 1,244,199
- chrom counts: {'1': 260225, '2': 244145, '3': 250406, '4': 208086, '5': 281337}
- context counts: {'CG': 1043262, 'CHG': 155647, 'CHH': 45290}

## Phenotype axis overlap

| key | value |
|---|---|
| n_twas | 145 |
| n_ewas | 114 |
| n_ppathway | 88 |
| n_gwas | 248 |
| n_gwaseqtl | 94 |
| twas_and_ewas | 79 |
| twas_and_ppathway | 82 |
| ewas_and_ppathway | 42 |
| twas_ewas_ppathway | 40 |
| twas_and_gwas | 119 |
| twas_and_gwaseqtl | 55 |
| ppathway_and_gwas | 70 |
| env_in_gwas | 0 |
| env_in_twas | 0 |
| env_in_ewas | 0 |
| env_in_gwaseqtl | 0 |
| env_vars | Elevation, Latitude, Longitude, bio1, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19, bio2 |
| example_twas_ewas | 0W GH FT, 2W, As75, B11, CL, Cd111, Co59, Cu65, DTF spain 2008 (2nd experiment), DTF spain 2009 (1st experiment), DTF spain 2009 (2nd experiment), DTF sweden 2008 (1st experiment), DTF sweden 2008 (2nd experiment), DTF sweden 2009 (1st experiment), DTF sweden 2009 (2nd experiment) |
| example_twas_ewas_ppath | 2W, As75, Cd111, Co59, Cu65, DTF spain 2008 (2nd experiment), DTF spain 2009 (1st experiment), DTF spain 2009 (2nd experiment), DTF sweden 2009 (1st experiment), DTFmainEffect2009, DTFplantingSummer2009, GR21, K39, Mg25, Mn55 |

## Gene overlap

| key | value |
|---|---|
| sizes | {'emqtl': 2645, 'ciseqtl': 2643, 'transeqtl': 580, 'twas': 16897, 'gwas_annot': 8376} |
| emqtl_and_ciseqtl | 387 |
| emqtl_and_twas | 1799 |
| twas_and_ciseqtl | 1694 |
| twas_and_transeqtl | 220 |
| emqtl_and_twas_and_ciseqtl | 257 |
| gwas_annot_and_twas | 5004 |
| gwas_annot_and_emqtl | 737 |

## Pathway overlap

```
{
  "n_epathway": 375,
  "n_ppathway": 403,
  "shared": 286,
  "shared_examples": [
    "ARG-PRO-PWY",
    "ARGASEDEG-PWY",
    "ARGDEG-V-PWY",
    "ARGSPECAT-PWY",
    "ARGSYN-PWY",
    "ASPARTATE-DEG1-PWY",
    "ASPARTATESYN-PWY",
    "ASPASN-ARA-PWY",
    "ASPSYNII-PWY",
    "BSUBPOLYAMSYN-PWY"
  ]
}
```

## Locus parsing and positional joins

| key | value |
|---|---|
| environmenteqtl_eQTL | 68,403 |
| ciseqtl_SNP | 11,796 |
| transeqtl_SNP | 10,118 |
| gwaseqtl_eQTL | 286,032 |
| gwaseqtl_GWAS | 286,032 |
| gwas_SNP | 89,807 |
| gwas_SNP_total | 89,884 |
| environmenteqtl_x_gwaseqtl_eQTL | 2,374 |
| n_environmenteqtl | 11,593 |
| n_gwaseqtl_eqlt | 2,661 |
| gwaseqtl GWAS loci also in gwas.csv | 6,968 |

## Study provenance

```
{
  "ewas": {
    "study16": 3803909,
    "study34": 138965,
    "study19": 27405,
    "study2": 11478,
    "study22": 5834,
    "study1": 2501,
    "study38": 1320,
    "study36": 713,
    "study12": 543,
    "study31": 372,
    "study4": 209,
    "study46": 125,
    "study29": 2,
    "study30": 2
  },
  "twas": {
    "study16": 53237,
    "study34": 5252,
    "study2": 3586,
    "study38": 232,
    "study1": 215,
    "study12": 115,
    "study19": 57,
    "study31": 34,
    "study43": 9,
    "study46": 9,
    "study29": 8
  },
  "ppathway": {
    "study16": 969,
    "study34": 327,
    "study2": 27,
    "study1": 13,
    "study19": 13,
    "study22": 3,
    "study18": 2
  },
  "gwas_numeric_study_n": 254
}
```