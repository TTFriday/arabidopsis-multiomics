"""
21_breadth_metric.py -- Disambiguate the breadth/effect-size relationship.

Two distinct quantities were being conflated:
  (i)  max|beta| per locus, which is taken over a variable number of tested
       environments and therefore grows mechanically with breadth;
  (ii) mean|beta| per locus, which is a per-locus average and does not.
The partial correlation must control for the number of tests in (i) but not in
(ii). This script computes both, plus a stratified analysis holding the number
of tested environments exactly fixed.
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

env = A.load("environmenteqtl")
env["abs_beta"] = env["Beta"].abs()

locus = env.groupby("eQTL").agg(
    n_env=("Environment", "nunique"),
    n_pairs=("abs_beta", "size"),
    max_abs_beta=("abs_beta", "max"),
    mean_abs_beta=("abs_beta", "mean"),
    median_abs_beta=("abs_beta", "median"),
)
locus["chrom"] = [str(x).split("_")[0] for x in locus.index]
print(f"loci: {len(locus):,}")
print("n_env == n_pairs for all loci:", bool((locus["n_env"] == locus["n_pairs"]).all()))


def partial_spearman(x, y, z):
    rx = pd.Series(np.asarray(x)).rank().values
    ry = pd.Series(np.asarray(y)).rank().values
    rz = pd.Series(np.asarray(z)).rank().values
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r, p = stats.pearsonr(ex, ey)
    return float(r), float(p)


def residualise_chrom(y, chrom):
    """Remove per-chromosome means and return residuals."""
    s = pd.Series(np.asarray(y, dtype=float))
    return (s - s.groupby(pd.Series(np.asarray(chrom)).values).transform("mean")).values


print("\n=== (i) max|beta| (extreme order statistic) ===")
raw = stats.spearmanr(locus["n_env"], locus["max_abs_beta"])
pr_npairs = partial_spearman(locus["n_env"], locus["max_abs_beta"], locus["n_pairs"])
print(f"raw Spearman(n_env, max|beta|)                 = {raw.statistic:+.4f} (p={raw.pvalue:.2e})")
print(f"partial, controlling for n_pairs               = {pr_npairs[0]:+.4f} (p={pr_npairs[1]:.3f})")
print("  -> n_pairs is identical to n_env, so this control removes the predictor;")
print("     the appropriate control is the number of DRAWs, which equals breadth itself.")
print("  -> correct phrasing: max|beta| is an extreme order statistic whose expectation")
print("     rises with breadth by construction.")

print("\n=== (ii) mean|beta| (per-locus average) ===")
raw2 = stats.spearmanr(locus["n_env"], locus["mean_abs_beta"])
# chromosome-adjusted, since loci are spatially clustered
res_b = residualise_chrom(locus["mean_abs_beta"], locus["chrom"])
res_n = residualise_chrom(locus["n_env"], locus["chrom"])
r_chrom, p_chrom = stats.pearsonr(res_n, res_b)
print(f"raw Spearman(n_env, mean|beta|)                = {raw2.statistic:+.4f} (p={raw2.pvalue:.2e})")
print(f"chromosome-adjusted Pearson                    = {r_chrom:+.4f} (p={p_chrom:.2e})")

print("\n=== (iii) stratified: breadth held exactly fixed is impossible for the max,")
print("          so compare mean|beta| across breadth strata ===")
grp = locus.groupby("n_env")["mean_abs_beta"].agg(["size", "mean", "median", "sem"])
grp["mean_minus_overall"] = grp["mean"] - locus["mean_abs_beta"].mean()
print(grp.to_string())
fit = np.polyfit(grp.index.values.astype(float), grp["mean"].values, 1)
print(f"\nOLS slope of stratum mean(|beta|) on breadth: {fit[0]:+.5f} per environment")
lo = grp.loc[grp.index.min(), "mean"]
hi = grp.loc[grp.index.max(), "mean"]
print(f"stratum mean at breadth {grp.index.min()}: {lo:.4f} | "
      f"at breadth {grp.index.max()}: {hi:.4f} | ratio {hi / lo:.3f}")

# How much of the max|beta| rise is predicted by order statistics alone?
print("\n=== (iv) order-statistic prediction for max|beta| ===")
print("If per-locus |beta| values were drawn from a common distribution, the expected")
print("maximum over k draws would rise with k. Compare observed vs a permutation in")
print("which each locus keeps its number of draws but the values are pooled.")
pool = locus["max_abs_beta"].values
rng = np.random.default_rng(5)
k = locus["n_env"].values.astype(int)
# fit a distribution-free approximation: sample k values from the pooled
# per-association |beta| distribution
vals = env["abs_beta"].values
pred = []
for _ in range(400):
    m = np.array([vals[rng.integers(0, len(vals), size=kk)].max() for kk in k])
    pred.append([m[k == b].mean() for b in sorted(set(k))])
pred = np.array(pred).mean(axis=0)
obs = [locus.loc[locus["n_env"] == b, "max_abs_beta"].mean() for b in sorted(set(k))]
comp = pd.DataFrame({"breadth": sorted(set(k)), "observed_mean_max": obs,
                     "expected_if_exchangeable": pred})
comp["ratio"] = comp["observed_mean_max"] / comp["expected_if_exchangeable"]
print(comp.to_string(index=False))
sp = stats.spearmanr(comp["expected_if_exchangeable"], comp["observed_mean_max"])
print(f"\nSpearman(expected under exchangeability, observed) = {sp.statistic:+.3f} (p={sp.pvalue:.2e})")

out = {
    "max_beta": {"raw_spearman": float(raw.statistic), "raw_p": float(raw.pvalue),
                 "partial_controlling_npairs": pr_npairs[0],
                 "partial_p": pr_npairs[1]},
    "mean_beta": {"raw_spearman": float(raw2.statistic), "raw_p": float(raw2.pvalue),
                  "chromosome_adjusted_pearson": float(r_chrom),
                  "chromosome_adjusted_p": float(p_chrom),
                  "stratum_slope": float(fit[0]),
                  "stratum_mean_min": float(lo), "stratum_mean_max": float(hi),
                  "stratum_ratio": float(hi / lo)},
    "order_statistic": {"spearman_expected_vs_observed": float(sp.statistic),
                        "p": float(sp.pvalue),
                        "table": comp.to_dict(orient="records")},
}
json.dump(out, open(A.RES / "breadth_metric.json", "w"), indent=2, default=str)
grp.to_csv(A.TAB / "env_breadth_mean_effect.csv")
comp.to_csv(A.TAB / "env_breadth_orderstat.csv", index=False)
print("\nsaved results/breadth_metric.json, tables/env_breadth_mean_effect.csv, "
      "tables/env_breadth_orderstat.csv")
