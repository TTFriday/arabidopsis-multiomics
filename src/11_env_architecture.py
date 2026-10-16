"""
11_env_architecture.py -- Analysis 1: architecture of environmental adaptation
signals.

Characterises the landscape of environment-associated loci (climate,
geography) and quantifies three architectural properties:

  A1  genomic clustering of environment-associated loci, benchmarked against a
      matched random sample of the assayed variant panel;
  A2  the environmental-breadth/effect-size relationship: do loci that associate
      with more environmental variables carry larger allele effects;
  A3  regulatory architecture: environment-associated loci are enriched for
      local (cis-eQTL) membership and depleted for distant (trans-eQTL)
      membership, tested against chromosome-matched random locus sets and
      repeated on a linkage-pruned locus set.

Outputs: tables/env_*.csv and results/env_architecture.json
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import alib as A  # noqa: E402

warnings.filterwarnings("ignore")

ENV_CATEGORY = {}
for v in ["bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9",
          "bio10", "bio11"]:
    ENV_CATEGORY[v] = "temperature"
for v in ["bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19"]:
    ENV_CATEGORY[v] = "precipitation"
for v in ["Latitude", "Longitude", "Elevation"]:
    ENV_CATEGORY[v] = "geography"


def sec(t):
    print("\n" + "=" * 84 + f"\n{t}\n" + "=" * 84, flush=True)


# --------------------------------------------------------------------- inputs
env = A.load("environmenteqtl")
cis = A.load("ciseqtl")
trans = A.load("transeqtl")
env["Category"] = env["Environment"].map(ENV_CATEGORY)
env["abs_beta"] = env["Beta"].abs()

cis_loci = set(cis["SNP"].astype(str))
trans_loci = set(trans["SNP"].astype(str))
env_loci_all = set(env["eQTL"].astype(str))
universe = env_loci_all | cis_loci | trans_loci

sec("0. Input summary")
print(f"environment-eQTL associations : {len(env):,}")
print(f"distinct environment loci     : {len(env_loci_all):,}")
print(f"environmental variables       : {env['Environment'].nunique()}")
print(f"cis-eQTL loci / trans-eQTL loci: {len(cis_loci):,} / {len(trans_loci):,}")
print(f"shared assayed locus universe  : {len(universe):,}")
print("\nassociations per environmental variable:")
per_env = (env.groupby("Environment")
           .agg(n_pairs=("eQTL", "size"), n_loci=("eQTL", "nunique"),
                mean_abs_beta=("abs_beta", "mean"), max_abs_beta=("abs_beta", "max"),
                median_p=("p_value", "median"), category=("Category", "first"))
           .sort_values("n_loci", ascending=False))
per_env["cis_fraction"] = (env.assign(is_cis=env["eQTL"].isin(cis_loci))
                           .groupby("Environment")["is_cis"].mean())
per_env["trans_fraction"] = (env.assign(is_trans=env["eQTL"].isin(trans_loci))
                             .groupby("Environment")["is_trans"].mean())
print(per_env.to_string())
per_env.to_csv(A.TAB / "env_variable_summary.csv")

sec("A1. Genomic clustering of environment-associated loci")
env_pos = A.locus_pos(env.drop_duplicates("eQTL"), "eQTL").dropna()
env_pos = env_pos.sort_values(["chrom", "pos"])
gap_rows = []
nn_env = []
for c, sub in env_pos.groupby("chrom"):
    p = np.sort(sub["pos"].values)
    if len(p) < 3:
        continue
    d = np.diff(p)
    nn = np.minimum(np.r_[np.inf, d], np.r_[d, np.inf])
    nn_env.append(pd.DataFrame({"chrom": c, "pos": p, "nn": nn}))
    gap_rows.append(dict(chrom=c, n_loci=len(p), median_gap=np.median(d),
                         mean_gap=d.mean(), frac_gap_lt1kb=(d < 1000).mean(),
                         frac_gap_lt5kb=(d < 5000).mean(),
                         span=p[-1] - p[0]))
nn_env = pd.concat(nn_env, ignore_index=True)
gaps = pd.DataFrame(gap_rows)
print(gaps.to_string(index=False))
print(f"\noverall: median nearest-neighbour distance = {nn_env['nn'].median():,.0f} bp; "
      f"fraction < 1 kb = {(nn_env['nn'] < 1000).mean():.3f}; "
      f"< 5 kb = {(nn_env['nn'] < 5000).mean():.3f}")

# matched null: random draws from the assayed variant panel
pool = pd.DataFrame({"locus": sorted(universe)})
pch, ppo = A.split_locus(pool["locus"])
pool["chrom"], pool["pos"] = pch, ppo
pool = pool.dropna()
pool_by_chrom = {c: np.sort(g["pos"].values) for c, g in pool.groupby("chrom")}
rng = np.random.default_rng(11)
n_loci = len(env_pos)
null_frac_1kb, null_med_nn = [], []
for _ in range(300):
    draw = []
    for c, g in env_pos.groupby("chrom"):
        arr = pool_by_chrom.get(c)
        if arr is None:
            continue
        draw.append(pd.DataFrame({"chrom": c,
                                  "pos": rng.choice(arr, size=len(g), replace=False)}))
    d = pd.concat(draw).sort_values(["chrom", "pos"])
    fr, med = [], []
    for c, g in d.groupby("chrom"):
        p = np.sort(g["pos"].values)
        if len(p) < 3:
            continue
        df_ = np.diff(p)
        fr.append((df_ < 1000).mean())
        nn = np.minimum(np.r_[np.inf, df_], np.r_[df_, np.inf])
        med.append(np.median(nn))
    null_frac_1kb.append(np.mean(fr))
    null_med_nn.append(np.median(med))
null_frac_1kb = np.array(null_frac_1kb)
null_med_nn = np.array(null_med_nn)
obs_frac_1kb = float((nn_env["nn"] < 1000).mean())
obs_med_nn = float(nn_env["nn"].median())
print(f"\nmatched random panel (n=300 draws):")
print(f"  fraction of loci with NN < 1 kb : observed {obs_frac_1kb:.3f} vs "
      f"null {null_frac_1kb.mean():.3f} +- {null_frac_1kb.std():.3f}  "
      f"z = {(obs_frac_1kb - null_frac_1kb.mean()) / null_frac_1kb.std():.1f}")
print(f"  median NN distance             : observed {obs_med_nn:,.0f} bp vs "
      f"null {null_med_nn.mean():,.0f} +- {null_med_nn.std():,.0f} bp")

# spatial autocorrelation of the environment association signal along chromosomes
sec("A1b. Autocorrelation of effect size along the chromosome")
acf_rows = []
for c, sub in env.sort_values("eQTL").groupby(env["eQTL"].map(lambda x: str(x).split("_")[0])):
    pos = A.split_locus(sub["eQTL"])[1]
    d = pd.DataFrame({"pos": pos, "ab": sub["abs_beta"].values}).dropna()
    d = d.groupby("pos")["ab"].max().sort_index()
    if len(d) < 50:
        continue
    for lag_kb in (1, 5, 20, 100, 500):
        x, y = [], []
        arr_p, arr_v = d.index.values, d.values
        j = np.searchsorted(arr_p, arr_p + lag_kb * 1000)
        ok = j < len(arr_p)
        if ok.sum() < 30:
            continue
        x, y = arr_v[ok], arr_v[j[ok]]
        acf_rows.append(dict(chrom=c, lag_kb=lag_kb, n=int(ok.sum()),
                             spearman=stats.spearmanr(x, y).statistic))
acf = pd.DataFrame(acf_rows)
print(acf.groupby("lag_kb")["spearman"].agg(["mean", "size"]).to_string())
acf.to_csv(A.TAB / "env_autocorrelation.csv", index=False)

sec("A2. Environmental breadth and effect size")
locus = (env.groupby("eQTL")
         .agg(n_env=("Environment", "nunique"), max_abs_beta=("abs_beta", "max"),
              mean_abs_beta=("abs_beta", "mean"), n_geo=("Category", lambda s: (s == "geography").sum()),
              n_temp=("Category", lambda s: (s == "temperature").sum()),
              n_prec=("Category", lambda s: (s == "precipitation").sum())))
locus["chrom"] = A.split_locus(locus.reset_index()["eQTL"])[0]
locus["is_cis"] = locus.index.isin(cis_loci)
locus["is_trans"] = locus.index.isin(trans_loci)
print(locus[["n_env", "max_abs_beta", "mean_abs_beta"]].describe().to_string())

sp = stats.spearmanr(locus["n_env"], locus["max_abs_beta"])
pr = stats.pearsonr(locus["n_env"], locus["max_abs_beta"])
print(f"\nSpearman(n_env, max|beta|) = {sp.statistic:+.4f}  p = {sp.pvalue:.3e}")
print(f"Pearson (n_env, max|beta|) = {pr.statistic:+.4f}  p = {pr.pvalue:.3e}")

grp = locus.groupby("n_env")["max_abs_beta"].agg(["size", "median", "mean"])
grp["sem"] = locus.groupby("n_env")["max_abs_beta"].sem()
print("\nmax|beta| by environmental breadth:")
print(grp.to_string())
grp.to_csv(A.TAB / "env_breadth_effect.csv")

# partial correlation controlling for chromosome and for the number of pairs
locus["n_pairs"] = env.groupby("eQTL").size()
def partial_spearman(x, y, z):
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    rz = pd.Series(z).rank().values
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return stats.pearsonr(ex, ey)

r_ctrl, p_ctrl = partial_spearman(locus["n_env"], locus["max_abs_beta"], locus["n_pairs"])
print(f"\npartial Spearman controlling for the number of tested environments: "
      f"r = {r_ctrl:.4f}  p = {p_ctrl:.3e}")

# slope from a linear model with chromosome fixed effects
import statsmodels.formula.api as smf  # noqa: E402
lm = smf.ols("max_abs_beta ~ n_env + C(chrom)", data=locus).fit()
print("\nOLS max|beta| ~ n_env + chromosome:")
print(f"  slope(n_env) = {lm.params['n_env']:+.4f}  SE = {lm.bse['n_env']:.4f}  "
      f"t = {lm.tvalues['n_env']:.1f}  p = {lm.pvalues['n_env']:.3e}  R2 = {lm.rsquared:.4f}")

sec("A3. Regulatory architecture of environment-associated loci")
rows = []


def arch_test(name, target, partner, partner_name, pruned=None):
    obs, null = A.matched_locus_null(target, partner, pool_by_chrom, n_perm=2000, seed=13)
    z, p = A.z_empirical(obs, null)
    rec = dict(target=name, partner=partner_name, n_target=len(target),
               observed=obs, expected=float(null.mean()), sd=float(null.std()),
               z=z, p_emp=p, ratio=obs / null.mean() if null.mean() else np.nan)
    if pruned is not None:
        obs2, null2 = A.matched_locus_null(pruned, partner, pool_by_chrom, n_perm=2000, seed=17)
        z2, p2 = A.z_empirical(obs2, null2)
        rec.update(pruned_observed=obs2, pruned_expected=float(null2.mean()),
                   pruned_z=z2, pruned_p=p2,
                   pruned_ratio=obs2 / null2.mean() if null2.mean() else np.nan)
    rows.append(rec)
    return rec


# linkage-pruned environment locus set (>=5 kb apart)
env_p = env.drop_duplicates("eQTL").sort_values("eQTL")
pch, ppo = A.split_locus(env_p["eQTL"])
env_p = pd.DataFrame({"locus": env_p["eQTL"].values, "chrom": pch, "pos": ppo}).dropna()
pruned = []
for c, g in env_p.sort_values("pos").groupby("chrom"):
    last = -np.inf
    for x in g["pos"].values:
        if x - last >= 5000:
            pruned.append(f"{c}_{int(x)}")
            last = x
pruned = set(pruned)

r1 = arch_test("environment", env_loci_all, cis_loci, "cis", pruned=pruned)
r2 = arch_test("environment", env_loci_all, trans_loci, "trans", pruned=pruned)
res = pd.DataFrame(rows)
print(res.to_string(index=False))

print("\nInterpretation anchors:")
print(f"  observed/expected for cis   : {r1['ratio']:.2f} (z={r1['z']:.1f}, p={r1['p_emp']:.4f})")
print(f"  observed/expected for trans : {r2['ratio']:.2f} (z={r2['z']:.1f}, p={r2['p_emp']:.4f})")
print(f"  after linkage pruning the cis excess changes from "
      f"{r1['ratio']:.2f}x to {r1['pruned_ratio']:.2f}x")

sec("A3b. Effect size and breadth by regulatory class")
cls = np.where(locus["is_cis"] & locus["is_trans"], "cis+trans",
               np.where(locus["is_cis"], "cis only",
                        np.where(locus["is_trans"], "trans only", "neither")))
locus["reg_class"] = cls
summary = (locus.groupby("reg_class")
           .agg(n=("n_env", "size"), median_n_env=("n_env", "median"),
                mean_n_env=("n_env", "mean"),
                median_max_abs_beta=("max_abs_beta", "median"),
                mean_max_abs_beta=("max_abs_beta", "mean")))
print(summary.to_string())

kw_b = stats.kruskal(*[g["max_abs_beta"].values for _, g in locus.groupby("reg_class")])
kw_n = stats.kruskal(*[g["n_env"].values for _, g in locus.groupby("reg_class")])
print(f"\nKruskal-Wallis max|beta| across regulatory classes: H={kw_b.statistic:.1f} p={kw_b.pvalue:.3e}")
print(f"Kruskal-Wallis n_env  across regulatory classes: H={kw_n.statistic:.1f} p={kw_n.pvalue:.3e}")
for a, b in (("cis only", "neither"), ("trans only", "neither"), ("cis only", "trans only")):
    x = locus.loc[locus["reg_class"] == a, "max_abs_beta"]
    y = locus.loc[locus["reg_class"] == b, "max_abs_beta"]
    u = stats.mannwhitneyu(x, y)
    print(f"  max|beta| {a} (n={len(x)}) vs {b} (n={len(y)}): "
          f"median {x.median():.3f} vs {y.median():.3f}, p={u.pvalue:.3e}")

locus.to_csv(A.TAB / "env_locus_architecture.csv")
res.to_csv(A.TAB / "env_architecture_tests.csv", index=False)
gaps.to_csv(A.TAB / "env_clustering_by_chromosome.csv", index=False)
summary.to_csv(A.TAB / "env_regulatory_class_summary.csv")

json.dump({
    "n_associations": int(len(env)),
    "n_env_loci": int(len(env_loci_all)),
    "n_env_variables": int(env["Environment"].nunique()),
    "n_universe": int(len(universe)),
    "clustering": {
        "observed_frac_nn_lt_1kb": obs_frac_1kb,
        "null_frac_nn_lt_1kb_mean": float(null_frac_1kb.mean()),
        "null_frac_nn_lt_1kb_sd": float(null_frac_1kb.std()),
        "observed_median_nn_bp": obs_med_nn,
        "null_median_nn_bp_mean": float(null_med_nn.mean()),
    },
    "breadth_effect": {
        "spearman_rho": float(sp.statistic), "spearman_p": float(sp.pvalue),
        "pearson_r": float(pr.statistic), "pearson_p": float(pr.pvalue),
        "partial_spearman_r": float(r_ctrl), "partial_spearman_p": float(p_ctrl),
        "ols_slope": float(lm.params["n_env"]), "ols_se": float(lm.bse["n_env"]),
        "ols_p": float(lm.pvalues["n_env"]), "ols_r2": float(lm.rsquared),
        "min_breadth": int(locus["n_env"].min()), "max_breadth": int(locus["n_env"].max()),
        "median_absbeta_at_min": float(grp.loc[locus["n_env"].min(), "median"]),
        "median_absbeta_at_max": float(grp.loc[locus["n_env"].max(), "median"]),
    },
    "architecture": rows,
    "pruned_loci": len(pruned),
}, open(A.RES / "env_architecture.json", "w"), indent=2, default=str)
print("\nsaved tables/env_*.csv and results/env_architecture.json")
