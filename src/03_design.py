"""
03_design.py -- Understand the experimental design / cohort structure.

The release contains only FDR-passing summary statistics, so hit counts are
driven by per-trait test-set size (a power proxy), not by architecture. This
script recovers the cohort structure, the per-trait test-set sizes, and the
overlap of the three-way core panel, so that downstream analyses can be
stratified on design rather than confounded by it.
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
OUT = Path(__file__).resolve().parents[1] / "results"


def rd(name, enc="utf-8", **kw):
    return pd.read_csv(DATA / name, encoding=enc, low_memory=False, **kw)


def sec(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78, flush=True)


ewas = rd("ewas.csv", usecols=["Phenotype", "Phenotype ontology", "Methylation", "Pearson", "Study"])
twas = rd("twas.csv")
ppath = rd("phenotypepathway.csv")
gwas = rd("gwas.csv")

sec("A. Study -> phenotype mapping per layer")
for nm, df in (("EWAS", ewas), ("TWAS", twas), ("Ppath", ppath)):
    print(f"\n--- {nm} ---")
    m = df.groupby("Study")["Phenotype"].agg(["nunique", "size"])
    print(m.sort_values("size", ascending=False).to_string())

sec("B. Do core-panel phenotypes all come from a single study?")
core = sorted(set(ewas["Phenotype"]) & set(twas["Phenotype"]) & set(ppath["Phenotype"]))
print("core n =", len(core))
st = pd.DataFrame({
    "ewas_study": ewas.groupby("Phenotype")["Study"].agg(lambda s: ",".join(sorted(set(s)))),
    "twas_study": twas.groupby("Phenotype")["Study"].agg(lambda s: ",".join(sorted(set(s)))),
    "ppath_study": ppath.groupby("Phenotype")["Study"].agg(lambda s: ",".join(sorted(set(s)))),
}).loc[core]
print(st.to_string())
print("\ncore phenotypes with identical study string in EWAS and TWAS:",
      int((st["ewas_study"] == st["twas_study"]).sum()), "/", len(core))

sec("C. Per-trait test-set size = cohort methylome/transcriptome coverage")
tn = ewas.groupby("Phenotype").size().rename("ewas_rows")
tt = twas.groupby("Phenotype").size().rename("twas_rows")
tot = pd.concat([tn, tt], axis=1).fillna(0).astype(int)
tot["ratio"] = tot["ewas_rows"] / tot["twas_rows"].replace(0, np.nan)
print(tot.describe())
print("\nMost common ewas_rows values (cohort coverage classes):")
print(tn.value_counts().head(20).to_string())

sec("D. Are the 40 core traits a balanced ion/rosette panel?")
ont = ewas.groupby("Phenotype")["Phenotype ontology"].first()
coreont = ont.loc[core]
print(coreont.value_counts().to_string())
print("\ncore traits list:")
print(core)

sec("E. Classify core traits into biological groups")
ION = re.compile(r"^[A-Z][a-z]?\d+$")
ions = [p for p in core if ION.match(p)]
climate = [p for p in core if p.lower().startswith("clim") or "DTF" in p or "GH" in p]
other = [p for p in core if p not in ions and p not in climate]
print("ions   :", len(ions), ions)
print("climate:", len(climate), climate)
print("other  :", len(other), other)

sec("F. Expression/methylation effect-size distributions for core traits")
e = ewas[ewas["Phenotype"].isin(core)]
t = twas[twas["Phenotype"].isin(core)]
print("EWAS |r| quantiles (core):", np.round(e["Pearson"].abs().quantile([.5, .75, .9, .99]).values, 3))
print("TWAS |r| quantiles (core):", np.round(t["Pearson"].abs().quantile([.5, .75, .9, .99]).values, 3))
print("\nmedian |r| per trait, EWAS vs TWAS (core panel):")
me = e.groupby("Phenotype")["Pearson"].apply(lambda s: s.abs().median()).rename("ewas_med_absr")
mt = t.groupby("Phenotype")["Pearson"].apply(lambda s: s.abs().median()).rename("twas_med_absr")
cmp = pd.concat([me, mt], axis=1)
cmp["diff"] = cmp["ewas_med_absr"] - cmp["twas_med_absr"]
print(cmp.sort_values("diff", ascending=False).to_string())
print("\nSummary: n traits where |r|_methylation > |r|_expression:", int((cmp["diff"] > 0).sum()), "/", len(cmp))

sec("G. Which phenotype-axis traits exist in GWAS but not in molecular layers")
gp = set(gwas["Phenotype"])
print("GWAS-only traits:", len(gp - set(ewas['Phenotype']) - set(twas['Phenotype'])))
print("molecular-layer traits absent from GWAS:", len((set(ewas['Phenotype']) | set(twas['Phenotype'])) - gp))

sec("H. Trait-name level cross-layer coherence: methylation<->expression rank corr")
sub = cmp.dropna()
print("Spearman between ewas_med_absr and twas_med_absr:", round(float(sub.corr(method="spearman").iloc[0, 1]), 3))
print("Pearson  :", round(float(sub.corr().iloc[0, 1]), 3))
print("log-count Spearman:",
      round(float(np.log10(tot.loc[sub.index, "ewas_rows"]).corr(
          np.log10(tot.loc[sub.index, "twas_rows"]), method="spearman")), 3))

OUT.mkdir(parents=True, exist_ok=True)
cmp.to_csv(OUT / "core_panel_effect_sizes.csv")
tot.to_csv(OUT / "per_trait_testset_sizes.csv")
json.dump({"core_traits": core, "ions": ions, "climate": climate, "other": other},
          open(OUT / "core_panel_groups.json", "w"), indent=2)
print("\nsaved results/core_panel_effect_sizes.csv, per_trait_testset_sizes.csv, core_panel_groups.json")
