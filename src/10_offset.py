"""
10_offset.py -- Test whether the released eQTL gene assignments are consistent
with a constant per-chromosome coordinate offset.

If SNP and gene coordinates came from different assembly builds, the error
d(SNP, reported gene) would be dominated by a per-chromosome constant that can
be estimated and removed. A successful fit would let the eQTL tables be used
positionally; a failed fit means the gene assignments must be treated as
identifier-only evidence.
"""
from __future__ import annotations

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


def split(s):
    ex = s.astype(str).str.extract(LOC)
    return ex[0].values, pd.to_numeric(ex[1], errors="coerce").values


GF = pd.read_csv(OUT / "gene_frame.csv")
gpos = GF.dropna(subset=["pos_med"]).set_index("Gene")["pos_med"]

cis = pd.read_excel(DATA / "ciseqtl.xls")
tr = pd.read_excel(DATA / "transeqtl.xls")

print("=" * 78, "\nOffset estimation per chromosome (cis pairs)\n", "=" * 78)
ch, po = split(cis["SNP"])
gchrom = "chr" + cis["Gene"].astype(str).str.extract(r"^AT(\d)G")[0].values
gp = cis["Gene"].map(gpos).values
res = pd.DataFrame({"snp_chrom": ch, "snp_pos": po, "gene_chrom": gchrom, "gene_pos": gp,
                    "gene": cis["Gene"].values})
res = res.dropna(subset=["snp_pos", "gene_pos"])
res["delta"] = res["snp_pos"] - res["gene_pos"]
print(f"usable pairs: {len(res):,}")
for c, sub in res.groupby("snp_chrom"):
    d = sub["delta"]
    # circular-mean style estimate: offset that would make |delta| smallest
    print(f"  {c}: n={len(sub):5d} delta median={d.median():>12,.0f} mean={d.mean():>13,.0f} "
          f"sd={d.std():>12,.0f} |sd/median|={abs(d.std() / d.median()) if d.median() else np.nan:.3f}")

print("\nIf an offset exists, |delta - median(delta)| should be small.")
for c, sub in res.groupby("snp_chrom"):
    d = sub["delta"]
    resid = (d - d.median()).abs()
    print(f"  {c}: median |resid| = {resid.median():>12,.0f} bp   "
          f"frac<100kb = {(resid < 1e5).mean():.3f}  frac<1Mb = {(resid < 1e6).mean():.3f}")

print("\n" + "=" * 78, "\nInteraction test: does delta depend on gene position?\n", "=" * 78)
for c, sub in res.groupby("snp_chrom"):
    if len(sub) > 50:
        r = stats.spearmanr(sub["gene_pos"], sub["delta"])
        print(f"  {c}: Spearman(gene_pos, delta) = {r.statistic:+.3f} (p={r.pvalue:.1e})")
        r2 = stats.spearmanr(sub["snp_pos"], sub["gene_pos"])
        print(f"      Spearman(snp_pos, gene_pos) = {r2.statistic:+.3f}")

print("\n" + "=" * 78, "\nIndependent check: does the gene-number ordering explain it?\n", "=" * 78)
gn = pd.to_numeric(cis["Gene"].astype(str).str.extract(r"^AT\dG(\d+)")[0], errors="coerce")
sub = pd.DataFrame({"gene_num": gn.values, "snp_pos": po, "gene_pos": gp}).dropna()
sub = sub[sub["snp_pos"].notna()]
print(f"n={len(sub):,}")
print("Spearman(gene_num, snp_pos) =", round(float(stats.spearmanr(sub["gene_num"], sub["snp_pos"]).statistic), 4))
print("Spearman(gene_num, gene_pos) =", round(float(stats.spearmanr(sub["gene_num"], sub["gene_pos"]).statistic), 4))

print("\n" + "=" * 78, "\nConclusion inputs\n", "=" * 78)
frac_near = ((res["delta"].abs()) < 1e5).mean()
print(f"fraction of cis pairs with |SNP - reported gene| < 100 kb : {frac_near:.4f}")
print(f"fraction with < 10 kb                                     : "
      f"{(res['delta'].abs() < 1e4).mean():.4f}")
print(f"fraction where SNP chromosome == reported gene chromosome  : "
      f"{(res['snp_chrom'] == res['gene_chrom']).mean():.4f}")
