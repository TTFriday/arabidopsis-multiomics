"""
05_annot.py -- Resolve the SNP annotation table and the true meaning of the
"cis"/"trans" labels, so that gene-level integration is done correctly.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
LOC = re.compile(r"^(chr[1-5])_(\d+)$")

gw = pd.read_csv(DATA / "gwas.csv", low_memory=False)
cis = pd.read_excel(DATA / "ciseqtl.xls")
ee = pd.read_excel(DATA / "environmenteqtl.xlsx")
ge = pd.read_csv(DATA / "gwaseqtl.csv", encoding="gb18030", low_memory=False)

print("=" * 78, "\nA. gwas.csv SNP->gene annotation: local or global?\n", "=" * 78)
g = gw.dropna(subset=["Genename"]).copy()
ex = g["SNP"].astype(str).str.extract(LOC)
g["snp_pos"] = pd.to_numeric(ex[1], errors="coerce")
gx = g["Genename"].astype(str).str.extract(r"^AT(\d)G(\d+)$")
g["gene_chr"] = "chr" + gx[0]
g["gene_pos"] = pd.to_numeric(gx[1], errors="coerce")
g["same_chr"] = g["SNP"].astype(str).str.extract(LOC)[0] == g["gene_chr"]
g["dist"] = (g["snp_pos"] - g["gene_pos"]).abs()
print("annotated rows:", len(g))
print("same chromosome fraction:", round(float(g["same_chr"].mean()), 4))
print("distance quantiles:")
print(g["dist"].describe(percentiles=[.25, .5, .75, .9]).to_string())
print("fraction within 5 kb:", round(float((g["dist"] <= 5000).mean()), 4))
print("fraction within 50 kb:", round(float((g["dist"] <= 50000).mean()), 4))
print("\nsample rows:")
print(g[["SNP", "Genename", "snp_pos", "gene_pos", "dist", "Effect", "Impact"]].head(8).to_string())

print("\n" + "=" * 78, "\nB. Do environment-eQTL loci overlap GWAS SNP annotation?\n", "=" * 78)
ee_snps = set(ee["eQTL"].astype(str))
annot_snps = set(g["SNP"].astype(str))
print("env-eQTL SNPs:", len(ee_snps), "| annotated SNPs in gwas.csv:", len(annot_snps),
      "| overlap:", len(ee_snps & annot_snps))
inter = sorted(ee_snps & annot_snps)
if inter:
    sub = g[g["SNP"].isin(inter)]
    print("\ndistance stats for the overlapping, annotated env-eQTL SNPs:")
    print(sub["dist"].describe(percentiles=[.25, .5, .75]).to_string())
    print("\nannotation types:")
    print(sub["Effect"].value_counts().to_string())
    print(sub[["SNP", "Genename", "dist", "Effect", "Impact"]].head(15).to_string())

print("\n" + "=" * 78, "\nC. gwas-eqtl: relationship between the joined eQTL and GWAS loci\n", "=" * 78)
lx = ge["eQTL loci"].astype(str).str.extract(LOC)
lg = ge["GWAS loci"].astype(str).str.extract(LOC)
same = lx[0].values == lg[0].values
dist = np.abs(pd.to_numeric(lx[1], errors="coerce") - pd.to_numeric(lg[1], errors="coerce"))
print("same chromosome fraction:", round(float(np.mean(same)), 4))
print("distance quantiles (bp):")
print(pd.Series(dist).describe(percentiles=[.25, .5, .75, .9]).to_string())
print("\nType x same-chromosome:")
print(pd.crosstab(ge["Type"], same).to_string())
print("\nType counts:", ge["Type"].value_counts().to_dict())
print("\nLD min/max:", ge["LD"].min(), ge["LD"].max(),
      "| n LD==0.5:", int((ge["LD"] <= 0.5001).sum()))

print("\n" + "=" * 78, "\nD. 'Alteration number' semantics\n", "=" * 78)
for nm, s in (("cis", cis["Alteration number"]), ("env", ee["Alteration"]),
              ("gwaseqtl-eQTL", ge["eQTL Alteration"])):
    v = pd.to_numeric(s, errors="coerce")
    print(f"{nm}: min={v.min()} max={v.max()} median={v.median()} "
          f"n_unique={v.nunique()}")

print("\n" + "=" * 78, "\nE. How many env-eQTL SNPs appear in the cis/trans eQTL tables?\n", "=" * 78)
tr = pd.read_excel(DATA / "transeqtl.xls")
print("env SNPs in cis table:", len(ee_snps & set(cis['SNP'].astype(str))))
print("env SNPs in trans table:", len(ee_snps & set(tr['SNP'].astype(str))))

print("\n" + "=" * 78, "\nF. Positional overlap: env-eQTL loci vs cis/trans eQTL SNPs (<1kb)\n", "=" * 78)
def pos(s):
    e = s.astype(str).str.extract(LOC)
    return pd.DataFrame({"c": e[0], "p": pd.to_numeric(e[1], errors="coerce")})

ep = pos(ee["eQTL"]).dropna().drop_duplicates()
cp = pos(cis["SNP"]).dropna().drop_duplicates()
tp = pos(tr["SNP"]).dropna().drop_duplicates()
from scipy.spatial import cKDTree  # noqa: E402
for nm, t in (("cis", cp), ("trans", tp)):
    hits = 0
    for c, grp in ep.groupby("c"):
        sub = t[t["c"] == c]
        if len(sub) == 0:
            continue
        tree = cKDTree(sub["p"].values.reshape(-1, 1))
        d, _ = tree.query(grp["p"].values.reshape(-1, 1), k=1)
        hits += int((d <= 1000).sum())
    print(f"env-eQTL loci within 1 kb of a {nm}-eQTL SNP: {hits} / {len(ep)}")
