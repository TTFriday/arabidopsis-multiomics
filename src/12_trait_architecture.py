"""
12_trait_architecture.py -- Analysis 2: trait-level architecture of methylation
versus expression associations.

The released tables are pre-filtered (every row already passes FDR <= 0.05), so
their row counts measure power, not architecture. This script therefore

  1. recovers FDR from the released p-values and re-derives a common
     significance threshold in every layer and every study, giving a
     like-for-like significant-feature count;
  2. restricts the methylation/expression contrast to traits evaluated in the
     same cohort (same `Study` in both EWAS and TWAS), which removes the
     cohort-size confound;
  3. tests whether the relative contribution of the two molecular layers is
     trait-class dependent, and whether the two layers agree on which traits
     are strongly associated.

Outputs: tables/trait_*.csv and results/trait_architecture.json
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import alib as A  # noqa: E402

warnings.filterwarnings("ignore")


def sec(t):
    print("\n" + "=" * 84 + f"\n{t}\n" + "=" * 84, flush=True)


ewas = A.load("ewas")
twas = A.load("twas")
epath = A.load("epathway")
ppath = A.load("phenotypepathway")
ewas = ewas.rename(columns={"Pvalue_Pearson": "p", "FDR_Pearson": "fdr"})
twas = twas.rename(columns={"Pvalue_Pearson": "p", "FDR_Pearson": "fdr"})
epath = epath.rename(columns={"Pvalue_Pearson": "p", "FDR_Pearson": "fdr"})
ppath = ppath.rename(columns={"Pvalue_Pearson": "p", "FDR_Pearson": "fdr"})

sec("1. Is the released FDR reproducible from the released p-values?")
repro = {}
for nm, df in (("EWAS", ewas), ("TWAS", twas), ("epathway", epath), ("ppath", ppath)):
    sub = df[["p", "fdr"]].dropna()
    # reconstruct per-phenotype / per-study BH and compare to the released value
    key = "Study" if "Study" in df.columns else None
    rec = {}
    if nm == "EWAS":
        grp = df.groupby("Phenotype")
    elif nm == "TWAS":
        grp = df.groupby("Phenotype")
    elif nm == "epathway":
        grp = df.groupby("Methylation")
    else:
        grp = df.groupby("Phenotype")
    sample_keys = list(grp.groups.keys())[:400]
    errs = []
    for k in sample_keys:
        g = grp.get_group(k)
        pv = g["p"].values
        q = A.bh_fdr(pv)
        errs.append(np.nanmax(np.abs(q - g["fdr"].values)))
    errs = np.array(errs)
    print(f"{nm}: released rows={len(df):,} | max |BH(p) - released FDR| over "
          f"{len(sample_keys)} groups: median={np.median(errs):.3g} p90={np.percentile(errs,90):.3g}")
    repro[nm] = dict(n_rows=int(len(df)), median_abs_diff=float(np.median(errs)),
                     p90_abs_diff=float(np.percentile(errs, 90)))

sec("2. Common-threshold significant-feature counts (global BH within layer)")
sig = {}
for nm, df, feat in (("EWAS", ewas, "Methylation"), ("TWAS", twas, "Gene"),
                     ("epathway", epath, "Methylation"), ("ppath", ppath, "Pathway")):
    q = A.bh_fdr(df["p"].values)
    df = df.assign(q_global=q)
    sig[nm] = df
    print(f"{nm}: p-value range [{df['p'].min():.3g}, {df['p'].max():.3g}] | "
          f"BH(p)<0.05: {(q < 0.05).sum():,}/{len(df):,} | released FDR<0.05: "
          f"{(df['fdr'] < 0.05).sum():,}/{len(df):,}")

sec("3. Cohort-matched methylation/expression trait panel")
# traits measured in the same study in both layers
se = ewas.groupby("Phenotype")["Study"].agg(lambda s: sorted(set(s)))
st = twas.groupby("Phenotype")["Study"].agg(lambda s: sorted(set(s)))
shared = []
for ph in sorted(set(se.index) & set(st.index)):
    common = set(se[ph]) & set(st[ph])
    if common:
        shared.append(ph)
print(f"traits with a shared cohort in EWAS and TWAS: {len(shared)}")

rows = []
for ph in shared:
    e = sig["EWAS"][sig["EWAS"]["Phenotype"] == ph]
    t = sig["TWAS"][sig["TWAS"]["Phenotype"] == ph]
    st_common = sorted(set(e["Study"]) & set(t["Study"]))
    e = e[e["Study"].isin(st_common)]
    t = t[t["Study"].isin(st_common)]
    if len(e) < 20 or len(t) < 20:
        continue
    e08 = e["q_global"] < 0.05
    t08 = t["q_global"] < 0.05
    rows.append(dict(
        Phenotype=ph, study=",".join(st_common), trait_class=A.trait_class(ph),
        ewas_n=int(len(e)), twas_n=int(len(t)),
        ewas_sig=int(e08.sum()), twas_sig=int(t08.sum()),
        ewas_frac_sig=float(e08.mean()), twas_frac_sig=float(t08.mean()),
        ewas_med_absr=float(e["Pearson"].abs().median()),
        twas_med_absr=float(t["Pearson"].abs().median()),
        ewas_sig_med_absr=float(e.loc[e08, "Pearson"].abs().median()) if e08.any() else np.nan,
        twas_sig_med_absr=float(t.loc[t08, "Pearson"].abs().median()) if t08.any() else np.nan,
        ewas_max_absr=float(e["Pearson"].abs().max()),
        twas_max_absr=float(t["Pearson"].abs().max()),
    ))
T = pd.DataFrame(rows)
T["log2_ratio_n"] = np.log2((T["ewas_sig"] + 1) / (T["twas_sig"] + 1))
T["diff_med_absr"] = T["ewas_med_absr"] - T["twas_med_absr"]
T["diff_sig_med_absr"] = T["ewas_sig_med_absr"] - T["twas_sig_med_absr"]
T["meth_dominant"] = T["diff_med_absr"] > 0
print(f"analysable cohort-matched trait pairs: {len(T)}")
print(T[["Phenotype", "study", "trait_class", "ewas_n", "twas_n", "ewas_sig", "twas_sig",
         "diff_med_absr"]].sort_values("diff_med_absr", ascending=False).to_string(index=False))
T.to_csv(A.TAB / "trait_layer_panel.csv", index=False)

sec("4. Primary contrast: significant features per layer")
print("significant-feature counts (global BH<0.05):")
print(T[["ewas_sig", "twas_sig"]].describe().to_string())
w = stats.wilcoxon(T["ewas_sig"], T["twas_sig"])
print(f"\nWilcoxon signed-rank EWAS vs TWAS significant features: "
      f"W={w.statistic:.0f} p={w.pvalue:.3e}")
print(f"traits where methylation has more significant features: "
      f"{int((T['ewas_sig'] > T['twas_sig']).sum())}/{len(T)}")
print("\nby trait class:")
byc = T.groupby("trait_class").agg(
    n=("Phenotype", "size"), med_ewas_sig=("ewas_sig", "median"),
    med_twas_sig=("twas_sig", "median"),
    n_meth_dom=("meth_dominant", "sum"),
    med_diff_absr=("diff_med_absr", "median"))
print(byc.to_string())
byc.to_csv(A.TAB / "trait_class_summary.csv")

sec("5. Effect-size contrast among significant features")
print(T[["ewas_sig_med_absr", "twas_sig_med_absr"]].describe().to_string())
w2 = stats.wilcoxon(T["ewas_sig_med_absr"].dropna(), T.loc[T["ewas_sig_med_absr"].notna(),
                                                          "twas_sig_med_absr"])
print(f"\nWilcoxon on median |r| among significant features: W={w2.statistic:.0f} p={w2.pvalue:.3e}")
print(f"traits where methylation |r| exceeds expression |r|: "
      f"{int((T['diff_sig_med_absr'] > 0).sum())}/{T['diff_sig_med_absr'].notna().sum()}")

sec("6. Do the two layers agree on which traits are strongly associated?")
for x, y, lab in (("ewas_n", "twas_n", "test-set size"),
                  ("ewas_med_absr", "twas_med_absr", "median |r|"),
                  ("ewas_sig", "twas_sig", "significant features")):
    sp = stats.spearmanr(T[x], T[y])
    print(f"Spearman({x}, {y}) [{lab}] = {sp.statistic:+.3f}  p = {sp.pvalue:.3e}")

sec("7. Covariate check: is the methylation excess explained by cohort size?")
T["log_ewas_n"] = np.log10(T["ewas_n"])
T["log_twas_n"] = np.log10(T["twas_n"])
T["log2_n_ratio"] = T["log_ewas_n"] - T["log_twas_n"]
sp1 = stats.spearmanr(T["log2_n_ratio"], T["log2_ratio_n"])
print(f"Spearman(log2 test-set ratio, log2 significant-feature ratio) = "
      f"{sp1.statistic:+.3f} (p={sp1.pvalue:.3e})")
partial = None
try:
    import statsmodels.formula.api as smf
    m = smf.ols("log2_ratio_n ~ log2_n_ratio + C(trait_class)", data=T).fit()
    print("\nOLS log2(sig ratio) ~ log2(test-set ratio) + trait class:")
    print(f"  log2 test-set ratio slope = {m.params['log2_n_ratio']:+.3f} "
          f"(p={m.pvalues['log2_n_ratio']:.3e})")
    print(f"  intercept = {m.params['Intercept']:+.3f} (p={m.pvalues['Intercept']:.3e})")
    print(f"  R2 = {m.rsquared:.3f}")
    partial = dict(slope=float(m.params["log2_n_ratio"]), p=float(m.pvalues["log2_n_ratio"]),
                   intercept=float(m.params["Intercept"]),
                   intercept_p=float(m.pvalues["Intercept"]), r2=float(m.rsquared))
except Exception as exc:  # pragma: no cover
    print("statsmodels unavailable:", exc)

sec("8. Replication across independent cohorts")
# traits present in >1 study in both layers
rep_rows = []
for layer, df, feat in (("EWAS", ewas, "Methylation"), ("TWAS", twas, "Gene")):
    sub = df[df["Phenotype"].isin(A.ION_TRAITS + A.CLIMATE_TRAITS)]
    g = sub.groupby(["Phenotype", "Study"]).agg(
        n=("Pearson", "size"), med_absr=("Pearson", lambda s: float(s.abs().median())),
        max_absr=("Pearson", lambda s: float(s.abs().max()))).reset_index()
    multi = g.groupby("Phenotype").filter(lambda x: len(x) > 1)
    for ph, gg in multi.groupby("Phenotype"):
        rep_rows.append(dict(layer=layer, Phenotype=ph, n_studies=len(gg),
                             studies=",".join(gg["Study"]),
                             med_absr_spread=float(gg["med_absr"].max() - gg["med_absr"].min())))
print(f"cross-cohort replication available for {len(set(r['Phenotype'] for r in rep_rows))} traits")
if rep_rows:
    R = pd.DataFrame(rep_rows)
    print(R.to_string(index=False))
    R.to_csv(A.TAB / "trait_cross_cohort_replication.csv", index=False)
else:
    print("no trait was evaluated in more than one cohort in both layers")

sec("9. Pathway layer: methylation-pathway and trait-pathway convergence")
shared_pw = set(epath["Pathway"]) & set(ppath["Pathway"])
print(f"pathways with both methylation and trait associations: {len(shared_pw)}")
trait_pw = ppath.groupby("Phenotype")["Pathway"].apply(set)
meth_pw = epath.groupby("Methylation")["Pathway"].apply(set)
print(f"phenotypes with pathway associations: {ppath['Phenotype'].nunique()}")
print(f"methylation sites with pathway associations: {epath['Methylation'].nunique():,}")
pw_stats = (ppath.groupby("Phenotype")
            .agg(n_pathway=("Pathway", "nunique"),
                 med_absr=("Pearson", lambda s: float(s.abs().median())),
                 max_absr=("Pearson", lambda s: float(s.abs().max())))
            .sort_values("n_pathway", ascending=False))
print(pw_stats.head(15).to_string())
pw_stats.to_csv(A.TAB / "phenotype_pathway_summary.csv")

json.dump({
    "fdr_reproducibility": repro,
    "n_rows_significant_global_bh": {k: int((v["q_global"] < 0.05).sum()) for k, v in sig.items()},
    "n_rows_released": {k: int(len(v)) for k, v in sig.items()},
    "cohort_matched_pairs": int(len(T)),
    "wilcoxon_sig_features": {"W": float(w.statistic), "p": float(w.pvalue)},
    "wilcoxon_sig_effect_size": {"W": float(w2.statistic), "p": float(w2.pvalue)},
    "n_traits_methylation_more_features": int((T["ewas_sig"] > T["twas_sig"]).sum()),
    "n_traits_methylation_larger_effect": int((T["diff_sig_med_absr"] > 0).sum()),
    "spearman_testset_size": {"rho": float(stats.spearmanr(T["ewas_n"], T["twas_n"]).statistic),
                              "p": float(stats.spearmanr(T["ewas_n"], T["twas_n"]).pvalue)},
    "spearman_effect_size": {"rho": float(stats.spearmanr(T["ewas_med_absr"], T["twas_med_absr"]).statistic),
                             "p": float(stats.spearmanr(T["ewas_med_absr"], T["twas_med_absr"]).pvalue)},
    "cohort_ratio_vs_sig_ratio": {"rho": float(sp1.statistic), "p": float(sp1.pvalue)},
    "ols_adjusted": partial,
    "trait_class_summary": byc.reset_index().to_dict(orient="records"),
    "n_shared_pathways": len(shared_pw),
}, open(A.RES / "trait_architecture.json", "w"), indent=2, default=str)
print("\nsaved tables/trait_*.csv and results/trait_architecture.json")
