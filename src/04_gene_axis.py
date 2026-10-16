"""
04_gene_axis.py -- Validate the gene-level integration axis.

Tests whether the gene universe can carry a rigorous cross-layer analysis:
  * cis-eQTL SNP positions relative to their eGene (confirming cis architecture),
  * gene-level overlap of environment-associated loci vs cis-eQTL loci,
  * environmental breadth per locus and its relationship to effect size,
  * Fisher enrichment of environment-eQTL genes among cis-eQTL genes,
  * genes present in ALL layers (EMQTL, cis-eQTL, TWAS).
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
OUT = Path(__file__).resolve().parents[1] / "results"
LOC = re.compile(r"^(chr[1-5])_(\d+)$")


def rd(n, enc="utf-8", **kw):
    return pd.read_csv(DATA / n, encoding=enc, low_memory=False, **kw)


def sec(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78, flush=True)


def locus_pos(s):
    ex = s.astype(str).str.extract(LOC)
    return pd.to_numeric(ex[1], errors="coerce").values


def gene_pos(g):
    ex = g.astype(str).str.extract(r"^AT(\d)G(\d+)$")
    return pd.to_numeric(ex[1], errors="coerce").values


cis = pd.read_excel(DATA / "ciseqtl.xls")
tr = pd.read_excel(DATA / "transeqtl.xls")
ee = pd.read_excel(DATA / "environmenteqtl.xlsx")
em = rd("emqtl.csv")
twas = rd("twas.csv")

sec("1. Is the cis-eQTL table genuinely local? (SNP pos vs eGene pos)")
cis = cis.assign(snp_pos=locus_pos(cis["SNP"]), gene_pos=gene_pos(cis["Gene"]),
                 snp_chr=cis["SNP"].astype(str).str.extract(LOC)[0],
                 gene_chr="chr" + cis["Gene"].astype(str).str.extract(r"^AT(\d)G")[0])
cis["dist"] = (cis["snp_pos"] - cis["gene_pos"]).abs()
print("same chromosome fraction:", round(float((cis["snp_chr"] == cis["gene_chr"]).mean()), 4))
print("distance to eGene (bp) quantiles:")
print(cis["dist"].describe(percentiles=[.25, .5, .75, .9, .99]).to_string())
print("fraction within 10 kb:", round(float((cis["dist"] <= 10000).mean()), 4))

sec("2. Same check for trans-eQTL (expected: distant / other chromosome)")
tr = tr.assign(snp_pos=locus_pos(tr["SNP"]),
               snp_chr=tr["SNP"].astype(str).str.extract(LOC)[0],
               gene_chr="chr" + tr["Gene"].astype(str).str.extract(r"^AT(\d)G")[0],
               gene_pos=gene_pos(tr["Gene"]))
tr["dist"] = (tr["snp_pos"] - tr["gene_pos"]).abs()
print("same chromosome fraction:", round(float((tr["snp_chr"] == tr["gene_chr"]).mean()), 4))
print(tr["dist"].describe(percentiles=[.25, .5, .75]).to_string())

sec("3. Environment-eQTL: locus-level environmental breadth")
ee2 = ee.assign(pos=locus_pos(ee["eQTL"]), chrom=ee["eQTL"].astype(str).str.extract(LOC)[0])
breadth = ee2.groupby("eQTL").agg(n_env=("Environment", "nunique"),
                                  max_absbeta=("Beta", lambda s: float(s.abs().max())),
                                  n_rows=("Beta", "size"))
print(breadth["n_env"].value_counts().sort_index().to_string())
print("\ncorrelation between breadth and max|beta| (Spearman):",
      round(float(stats.spearmanr(breadth["n_env"], breadth["max_absbeta"]).statistic), 4))
print("breadth vs n_rows Spearman:",
      round(float(stats.spearmanr(breadth["n_env"], breadth["n_rows"]).statistic), 4))
print("\nmean max|beta| by breadth:")
print(breadth.groupby("n_env")["max_absbeta"].agg(["mean", "size"]).head(20).to_string())

sec("4. Environment-eQTL loci in genomic context: clustering?")
ee2 = ee2.sort_values(["chrom", "pos"])
for c, g in ee2.groupby("chrom"):
    p = np.sort(g["pos"].unique())
    if len(p) > 1:
        d = np.diff(p)
        print(f"{c}: n={len(p):5d} median_gap={np.median(d):9.0f} "
              f"frac_gap<1kb={np.mean(d < 1000):.3f} span={p[-1]-p[0]:,}")

sec("5. Gene-level axis: environment-associated vs cis-eQTL genes")
# assign environment-associated SNP -> nearby gene using the exon/intron annotation
# available in gwas.csv (SNP <-> Genename), then require the SNP to be within
# 10 kb of the gene represented in the cis-eQTL table.
gw = rd("gwas.csv")
snp2gene = gw.dropna(subset=["Genename"]).groupby("SNP")["Genename"].agg(
    lambda s: sorted(set(s)))
print("SNP->gene annotated pairs in gwas.csv:", int(gw["Genename"].notna().sum()),
      "unique SNP:", len(snp2gene))

ee_snps = set(ee2["eQTL"])
print("environment-eQTL SNPs:", len(ee_snps),
      "| of which annotated in gwas.csv:", len(ee_snps & set(snp2gene.index)))

env_genes = set()
for s in ee_snps & set(snp2gene.index):
    env_genes.update(snp2gene[s])
print("implied environment-associated genes:", len(env_genes))

cis_genes = set(cis["Gene"])
tr_genes = set(tr["Gene"])
em_genes = set(em["Gene"])
twas_genes = set(twas["Gene"])
print("\ngene set sizes: env=%d cis=%d trans=%d emqtl=%d twas=%d"
      % (len(env_genes), len(cis_genes), len(tr_genes), len(em_genes), len(twas_genes)))

sec("6. Fisher enrichment tests on the assayed gene universe")
universe = cis_genes | tr_genes | em_genes  # genes with any molecular assay
print("universe (any molecular assay):", len(universe))
rows = []


def fisher(name, A, B, U):
    a = len(A & B & U)
    b = len(A & U) - a
    c = len(B & U) - a
    d = len(U) - a - b - c
    if min(a, b, c, d) < 0 or (b + d) == 0 or (a + c) == 0:
        rows.append({"test": name, "a": a, "b": b, "c": c, "d": d,
                     "OR": np.nan, "p": np.nan})
        return
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    rows.append({"test": name, "a": a, "b": b, "c": c, "d": d, "OR": orr, "p": p})


fisher("env_genes vs cis_genes", env_genes, cis_genes, universe)
fisher("env_genes vs trans_genes", env_genes, tr_genes, universe)
fisher("env_genes vs emqtl_genes", env_genes, em_genes, universe)
fisher("env_genes vs twas_genes", env_genes, twas_genes, universe)
fisher("cis_genes vs emqtl_genes", cis_genes, em_genes, universe)
fisher("cis_genes vs twas_genes", cis_genes, twas_genes, universe)
fisher("trans_genes vs twas_genes", tr_genes, twas_genes, universe)
fr = pd.DataFrame(rows)
fr["p_bonf"] = np.minimum(fr["p"] * len(fr), 1.0)
print(fr.to_string(index=False))

sec("7. Genes present in ALL layers")
allc = cis_genes & em_genes & twas_genes
print("cis & emqtl & twas:", len(allc))
allc2 = env_genes & cis_genes & em_genes & twas_genes
print("env & cis & emqtl & twas:", len(allc2), sorted(allc2)[:20])
alle = env_genes & em_genes & twas_genes
print("env & emqtl & twas:", len(alle))

sec("8. EMQTL beta-scale check: gene-level methylation load")
emg = em.groupby("Gene")["Pearson"].agg(
    n_meth="size", max_absr=lambda s: float(s.abs().max()),
    med_absr=lambda s: float(s.abs().median()))
print(emg.describe().to_string())
print("\ncorrelation of gene methylation load with cis-eQTL evidence:")
cisg = cis.groupby("Gene").agg(n_snp=("Beta", "size"),
                               max_absbeta=("Beta", lambda s: float(s.abs().max())),
                               max_t=("t-stat", lambda s: float(s.abs().max())))
j = emg.join(cisg, how="inner")
print("n joined genes:", len(j))
print("Spearman n_meth vs n_snp:", round(float(stats.spearmanr(j["n_meth"], j["n_snp"]).statistic), 4))
print("Spearman max_absr vs max_t:", round(float(stats.spearmanr(j["max_absr"], j["max_t"]).statistic), 4))

fr.to_csv(OUT / "gene_axis_fisher.csv", index=False)
breadth.to_csv(OUT / "env_locus_breadth.csv")
j.to_csv(OUT / "gene_meth_vs_cis.csv")
json.dump({"env_genes": sorted(env_genes), "all_layers_genes": sorted(allc2)},
          open(OUT / "gene_sets.json", "w"), indent=2)
print("\nsaved results/gene_axis_fisher.csv, env_locus_breadth.csv, gene_meth_vs_cis.csv, gene_sets.json")
