"""
02_probe.py -- Decisive feasibility probes.

Checks the details that determine which paper is actually writable:
  * per-phenotype significant-hit counts per omics layer (is n>=3 per phen?),
  * empirical null from extreme p-values (do we need a global FDR redo?),
  * environment-eQTL structure (SNP->environment, or SNP->gene->environment?),
  * GWAS "Score" semantics and p-value feasibility,
  * whether per-phenotype effect sizes are comparable across layers.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
OUT = Path(__file__).resolve().parents[1] / "results"

FDR = 0.05


def rd(name, enc="utf-8", **kw):
    return pd.read_csv(DATA / name, encoding=enc, low_memory=False, **kw)


def section(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78, flush=True)


# ------------------------------------------------------------------ 1. EWAS
section("1. EWAS: per-phenotype hit counts and p-value floor")
ewas = rd("ewas.csv", usecols=["Phenotype", "Phenotype ontology", "Pearson",
                               "Pvalue_Pearson", "FDR_Pearson", "Study"])
print("min p-value overall:", ewas["Pvalue_Pearson"].min())
print("max p-value overall:", ewas["Pvalue_Pearson"].max())
print("n rows with FDR<=0.05:", int((ewas["FDR_Pearson"] <= FDR).sum()))
g = ewas.groupby("Phenotype").agg(
    n_tests=("Pearson", "size"),
    n_sig=("FDR_Pearson", lambda s: int((s <= FDR).sum())),
    med_r=("Pearson", "median"),
    max_absr=("Pearson", lambda s: float(s.abs().max())),
)
print("\nTests per phenotype: describe")
print(g["n_tests"].describe())
print("\nSignificant hits per phenotype: describe")
print(g["n_sig"].describe())
print("\nphenotypes with n_sig >= 30:", int((g["n_sig"] >= 30).sum()), "/", len(g))
print("\nTop 12 phenotypes by hit count")
print(g.sort_values("n_sig", ascending=False).head(12))

# empirical p-value floor / granularity
print("\nmost frequent p-values in EWAS (top 5):")
print(ewas["Pvalue_Pearson"].round(12).value_counts().head(5))

# ------------------------------------------------------------------ 2. TWAS
section("2. TWAS: per-phenotype hit counts")
twas = rd("twas.csv")
print("n rows:", len(twas), "| min p:", twas["Pvalue_Pearson"].min(),
      "| max p:", twas["Pvalue_Pearson"].max())
gt = twas.groupby("Phenotype").agg(
    n_tests=("Pearson", "size"),
    n_sig=("FDR_Pearson", lambda s: int((s <= FDR).sum())),
    max_absr=("Pearson", lambda s: float(s.abs().max())),
)
print(gt["n_tests"].describe())
print(gt["n_sig"].describe())
print("phenotypes with n_sig>=10:", int((gt["n_sig"] >= 10).sum()), "/", len(gt))
print(gt.sort_values("n_sig", ascending=False).head(10))

# --------------------------------------------------- 3. shared phenotype panel
section("3. Shared phenotype panel across EWAS / TWAS / phenotype-pathway")
ppath = rd("phenotypepathway.csv")
pe = set(ewas["Phenotype"].unique())
pt = set(twas["Phenotype"].unique())
pp = set(ppath["Phenotype"].unique())
core = pe & pt & pp
print(f"|EWAS|={len(pe)} |TWAS|={len(pt)} |Ppath|={len(pp)} -> three-way core = {len(core)}")
d = g.loc[sorted(core), ["n_sig"]].rename(columns={"n_sig": "ewas_sig"}).join(
    gt.loc[sorted(core), ["n_sig"]].rename(columns={"n_sig": "twas_sig"}))
print("\nCore-panel hit counts (head):")
print(d.sort_values("ewas_sig", ascending=False).head(20))
print("\nCore panel: phenotypes where BOTH layers have >=30 hits:",
      int(((d["ewas_sig"] >= 30) & (d["twas_sig"] >= 30)).sum()))
print("Core panel: EWAS-only (>=30, twas<10):",
      int(((d["ewas_sig"] >= 30) & (d["twas_sig"] < 10)).sum()))
print("Core panel: TWAS-only (>=10, ewas<30):",
      int(((d["twas_sig"] >= 10) & (d["ewas_sig"] < 30)).sum()))

# ------------------------------------------------------ 4. environment eQTL
section("4. environmenteqtl structure")
ee = pd.read_excel(DATA / "environmenteqtl.xlsx")
print(ee.head(3).to_string())
print("\nEnvironment levels:", sorted(ee["Environment"].unique()))
print("\ntype x Environment (top):")
print(pd.crosstab(ee["Environment"], ee["type"]).head(25))
print("\np-value range:", ee["p-value"].min(), ee["p-value"].max())
print("FDR range:", ee["FDR"].min(), ee["FDR"].max())
print("n sig FDR<=0.05:", int((ee["FDR"] <= FDR).sum()), "/", len(ee))
print("\nHow many distinct SNPs per environment:")
print(ee.groupby("Environment")["eQTL"].nunique().sort_values(ascending=False))
print("\nHow many environments per SNP (describe):")
print(ee.groupby("eQTL")["Environment"].nunique().describe())
print("\nBeta distribution by type:")
print(ee.groupby("type")["Beta"].describe())

# ------------------------------------------------------------------- 5. GWAS
section("5. GWAS Score semantics")
gw = rd("gwas.csv")
print(gw.head(4).to_string())
print("\nScore distribution:")
print(gw["Score"].describe())
print("\nScore vs MAF corr:", gw[["Score", "MAF"]].corr().iloc[0, 1])
print("Study values: min", gw["Study"].min(), "max", gw["Study"].max(),
      "n unique", gw["Study"].nunique())
print("\nPhenotype x study (n unique studies per phenotype):")
print(gw.groupby("Phenotype")["Study"].nunique().describe())
print("\nCan we recover n from MAF? MAF unique values:", gw["MAF"].nunique(),
      "min", gw["MAF"].min())
print("MAF granularity examples:", sorted(gw["MAF"].unique())[:8])
print("\nSNP rows with missing Chr (malformed SNP id):")
print(gw.loc[gw["Chr"].isna(), "SNP"].head(10).tolist())
print("\nPhenotype ontology counts:")
print(gw["Phenotype ontology"].value_counts().head(20))

# ------------------------------------------------- 6. ciseqtl / transeqtl
section("6. cis/trans eQTL tables")
cis = pd.read_excel(DATA / "ciseqtl.xls")
tr = pd.read_excel(DATA / "transeqtl.xls")
for nm, df in (("cis", cis), ("trans", tr)):
    print(f"\n--- {nm} ---")
    print(df.head(3).to_string())
    print("genes:", df["Gene"].nunique(), "snps:", df["SNP"].nunique(),
          "| Beta range", round(df["Beta"].min(), 3), round(df["Beta"].max(), 3))
    print("p range:", df["p-value"].min(), df["p-value"].max())
    print("FDR<=0.05:", int((df["FDR"] <= FDR).sum()), "/", len(df))
    print("'Alteration number' describe:", df["Alteration number"].describe().to_dict())
    print("beta sign frac positive:", round(float((df["Beta"] > 0).mean()), 3))

section("7. EMQTL details")
em = rd("emqtl.csv")
print("Gene->methylation rows:", len(em), "genes:", em["Gene"].nunique(),
      "meth:", em["Methylation"].nunique())
print("Pearson range:", em["Pearson"].min(), em["Pearson"].max())
print("sign positive frac:", round(float((em["Pearson"] > 0).mean()), 3))
print("rows per gene describe:")
print(em.groupby("Gene").size().describe())
print("rows per methylation describe:")
print(em.groupby("Methylation").size().describe())

section("8. methyl context distribution across tables")
import re
MRE = re.compile(r"^(\d+)_(\d+)_(CG|CHG|CHH)$")
for nm, s in (("ewas", ewas["Phenotype"]),):
    pass
ewas_full = rd("ewas.csv", usecols=["Methylation"])
for nm, ser in (("ewas", ewas_full["Methylation"]), ("emqtl", em["Methylation"]),
                ("epathway", rd("epathway.csv", usecols=["Methylation"])["Methylation"])):
    ex = ser.astype(str).str.extract(MRE)
    print(f"{nm}: n={len(ser):,} contexts={ex[2].value_counts().to_dict()}")
