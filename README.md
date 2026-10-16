# Arabidopsis Multi-Omics Association Analysis

Analysis pipeline for integrative re-analysis of the AtMAD (Arabidopsis thaliana Multi-omics Association Database) release, examining the layered architecture of environmental adaptation and trait association.

## Overview

This repository contains the complete computational pipeline for analyzing multi-omics association data from Arabidopsis thaliana, including:

- **Data audit**: Multiple-testing structure and cohort design validation
- **Environmental adaptation**: Architecture of climate and geography-associated loci across 22 environmental variables
- **Layer comparison**: Concordance analysis between methylation, expression and pathway layers
- **Robustness testing**: Null models and matched controls for regulatory enrichment

## Project Structure

```
├── src/                          # Analysis scripts
│   ├── alib.py                   # Shared utilities
│   ├── 00_profile.py             # Data profiling and inventory
│   ├── 01_overlap.py             # Cross-layer overlap analysis
│   ├── 02_probe.py               # Test-set structure investigation
│   ├── 03_design.py              # Cohort and design matrix analysis
│   ├── 04_gene_axis.py           # Gene coordinate frame recovery
│   ├── 05_annot.py               # Functional annotation analysis
│   ├── 06_gene_pos.py            # Gene position calibration
│   ├── 07_coords_and_tests.py    # Statistical testing framework
│   ├── 08_enrich.py              # Regulatory enrichment analysis
│   ├── 09_validate.py            # Validation and robustness checks
│   ├── 10_offset.py              # Offset and normalization analysis
│   ├── 11_env_architecture.py    # Environmental locus architecture
│   ├── 12_trait_architecture.py  # Trait-level architecture analysis
│   ├── 13_convergence.py         # Cross-layer convergence analysis
│   ├── 21_breadth_metric.py      # Environmental breadth metrics
│   └── 22_null_models.py         # Null models and permutation tests
├── tables/                       # Intermediate analysis tables (CSV)
├── results/                      # Summary statistics and metrics (JSON/CSV)
└── figures/                      # Analysis figures (PDF/EPS)
```

## Data Source

All analyses use data from the [Arabidopsis thaliana Multi-omics Association Database (AtMAD)](https://plantomics.ndbc.ac.cn/AtMAD/), which aggregates:

- Genome-wide association studies (GWAS)
- Expression quantitative trait loci (eQTL): cis and trans
- Environmental eQTL (environment-associated loci)
- Epigenome-wide association studies (EWAS)
- Transcriptome-wide association studies (TWAS)
- Expression-methylation QTL (EMQTL)
- Pathway-level associations

## Analysis Pipeline

The scripts are numbered to indicate execution order:

1. **Data profiling** (00-03): Inventory, overlap analysis, test-set structure, cohort design
2. **Coordinate systems** (04-06): Gene axis recovery, annotation, position calibration
3. **Statistical analysis** (07-10): Testing framework, enrichment, validation, normalization
4. **Architecture analysis** (11-13): Environmental and trait-level architecture, convergence
5. **Extended analysis** (21-22): Breadth metrics, null models

## Key Findings

- All released tables are pre-filtered discovery sets (FDR < 0.05)
- Row counts reflect testing power, not association architecture
- Environmental-eQTL loci are strongly clustered (57% within 1 kb vs 35% expected)
- Environment-associated loci enrich for cis-regulatory variation (1.54-fold, P = 5×10⁻⁴)
- Methylation and expression layers show strong effect size concordance (ρ = 0.62-0.87)

## Requirements

```bash
pip install numpy scipy pandas statsmodels matplotlib
```

## Usage

Run analyses in numerical order:

```bash
cd src
python 00_profile.py
python 01_overlap.py
# ... continue through the pipeline
```

Each script reads from the AtMAD data files and writes intermediate results to `tables/` and summary statistics to `results/`.

## Citation

If you use this analysis pipeline, please cite:

> Tu, L. & Peng, Q. (2026). Layered architecture of environmental adaptation and trait association in Arabidopsis thaliana: an integrative audit and multi-layer analysis of the AtMAD multi-omics association release.

## License

This analysis code is provided for research and educational purposes.


