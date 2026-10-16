"""
06_gene_pos.py -- Attempt to recover gene coordinates from the functional
annotation embedded in gwas.csv, so that positional cross-layer integration
becomes possible without an external GFF.

Logic: rows whose Effect is a coding/UTR/intron consequence must lie inside the
gene they name, so the distribution of their SNP positions brackets that gene.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
OUT = Path(__file__).resolve().parents[1] / "results"
LOC = re.compile(r"^(chr[1-5])_(\d+)$")

INTRAGENIC = {"INTRON", "EXON", "SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING",
              "STOP_GAINED", "STOP_LOST", "START_LOST", "SPLICE_SITE_ACCEPTOR",
              "SPLICE_SITE_DONOR", "SPLICE_SITE_REGION", "UTR_3_PRIME",
              "UTR_5_PRIME", "SPLICE_SITE"}

gw = pd.read_csv(DATA / "gwas.csv", low_memory=False)
print("Effect categories:")
print(gw["Effect"].value_counts(dropna=False).to_string())

g = gw.dropna(subset=["Genename"]).copy()
g["Effect"] = g["Effect"].astype(str).str.upper()
ex = g["SNP"].astype(str).str.extract(LOC)
g["chrom"] = ex[0]
g["snp_pos"] = pd.to_numeric(ex[1], errors="coerce")
gx = g["Genename"].astype(str).str.extract(r"^AT(\d)G(\d+)$")
g["gene_num"] = pd.to_numeric(gx[1], errors="coerce")
g["gene_chr"] = gx[0]

print("\nrows with parsable chrom+pos:", int(g["snp_pos"].notna().sum()), "/", len(g))

strict = g[g["Effect"].isin(INTRAGENIC)].copy()
print("strictly intragenic-annotated rows:", len(strict))

est = strict.groupby(["Genename", "chrom", "gene_chr", "gene_num"]).agg(
    n=("snp_pos", "size"),
    pos_min=("snp_pos", "min"),
    pos_max=("snp_pos", "max"),
    pos_med=("snp_pos", "median"),
).reset_index()
est["same_chr"] = est["chrom"] == ("chr" + est["gene_chr"])
est["span"] = est["pos_max"] - est["pos_min"]
print("\nestimates:", len(est), "| on matching chromosome:", int(est["same_chr"].sum()))
est = est[est["same_chr"]].copy()
print("gene-position spans (bp):")
print(est["span"].describe(percentiles=[.25, .5, .75, .9]).to_string())
print("\nspans <= 20 kb:", int((est["span"] <= 20000).sum()), "/", len(est))

# validate monotonicity: within a chromosome, gene number should track position
print("\nValidation - Spearman(gene_num, estimated position) within chromosome:")
for c, sub in est.groupby("chrom"):
    if len(sub) > 20:
        r = sub[["gene_num", "pos_med"]].corr(method="spearman").iloc[0, 1]
        print(f"  {c}: n={len(sub):5d} spearman={r:+.4f} "
              f"num_range=({int(sub.gene_num.min())},{int(sub.gene_num.max())}) "
              f"pos_range=({int(sub.pos_med.min()):,},{int(sub.pos_med.max()):,})")

# gene-number -> position regression per chromosome (for imputing all genes)
print("\nLinear fit position ~ gene_num per chromosome (for imputation):")
fits = {}
for c, sub in est.groupby("chrom"):
    if len(sub) < 30:
        continue
    x = sub["gene_num"].values.astype(float)
    y = sub["pos_med"].values.astype(float)
    A = np.polyfit(x, y, 1)
    pred = np.polyval(A, x)
    resid = y - pred
    r2 = 1 - resid.var() / y.var()
    fits[c] = A
    print(f"  {c}: n={len(sub):5d} slope={A[0]:.1f} intercept={A[1]:.0f} "
          f"R2={r2:.4f} resid_sd={resid.std():.0f}")

# use the fit to sanity-check the strict annotation set
print("\nHow well does the monotone mapping place known genes? (residual quantiles)")
for c, sub in est.groupby("chrom"):
    if c not in fits:
        continue
    pred = np.polyval(fits[c], sub["gene_num"].values.astype(float))
    r = sub["pos_med"].values - pred
    print(f"  {c}: |resid| median={np.median(np.abs(r)):,.0f} p90={np.percentile(np.abs(r),90):,.0f}")

est.to_csv(OUT / "gene_position_estimates.csv", index=False)
print("\nsaved results/gene_position_estimates.csv")
