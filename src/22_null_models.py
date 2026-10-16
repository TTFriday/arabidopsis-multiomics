"""
22_null_models.py -- Locus-level null models and robustness analyses.

All environment-associated loci in the release are themselves cis- or
trans-eQTL loci, so the natural reference set is the eQTL locus panel
(cis | trans). This script re-tests the locus-level questions against that
panel with sampling *without* replacement and with a label-rotation null that
preserves the physical clustering of loci along each chromosome:

  N1  composition of environment-associated loci (cis / trans) versus all
      other eQTL loci: Fisher tests, label-rotation null, 5-kb block collapse;
  N2  overlap with GWAS association loci and GWAS-eQTL partner loci;
  N3  genomic clustering, using the same nearest-neighbour statistic for the
      observed set and for random subsets of the panel;
  N4  effect size and breadth by regulatory class (maximum and mean |beta|);
  N5  gene-level proximity analysis restricted to genes whose coordinates were
      measured directly (no imputed positions);
  N6  trait-level hit rates normalised by the number of distinct features each
      layer contributes to the release;
  N7  descriptive checks of the GWAS-eQTL join and the three-layer trait panel;
  N8  the correlation floor applied to the trait tables, and effect sizes
      measured above that floor.

Outputs: tables/null_*.csv and results/null_models.json
"""
from __future__ import annotations

import itertools
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
ROOT = Path(__file__).resolve().parents[1]
RES, TAB = ROOT / "results", ROOT / "tables"
LOCUS_RE = re.compile(r"^(chr[1-5])_(\d+)$")
N_ROT = 2000
N_SUB = 1000
BLOCK_GAP = 5000
WINDOW = 10000


def sec(t):
    print("\n" + "=" * 80 + f"\n{t}\n" + "=" * 80, flush=True)


def read_xls(name):
    """Legacy BIFF workbooks; read with xlrd directly to avoid version pins."""
    import xlrd
    sh = xlrd.open_workbook(str(DATA / name)).sheet_by_index(0)
    rows = [sh.row_values(i) for i in range(sh.nrows)]
    return pd.DataFrame(rows[1:], columns=[str(c).strip() for c in rows[0]])


def split_locus(s):
    ex = pd.Series(s).astype(str).str.extract(LOCUS_RE)
    return ex[0].values, pd.to_numeric(ex[1], errors="coerce").values


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(r, 0, 1)
    return out


def nn_stats(df):
    """Fraction of loci with a nearest neighbour < 1 kb, and median NN (bp)."""
    nn = []
    for _, g in df.groupby("chrom"):
        p = np.sort(g["pos"].values)
        if len(p) < 2:
            continue
        d = np.diff(p).astype(float)
        nn.append(np.minimum(np.r_[np.inf, d], np.r_[d, np.inf]))
    nn = np.concatenate(nn)
    return float((nn < 1000).mean()), float(np.median(nn))


# ---------------------------------------------------------------- inputs
env = pd.read_excel(DATA / "environmenteqtl.xlsx")
env.columns = [str(c).strip() for c in env.columns]
cis = read_xls("ciseqtl.xls")
trans = read_xls("transeqtl.xls")
gwas = pd.read_csv(DATA / "gwas.csv", low_memory=False)
ge = pd.read_csv(DATA / "gwaseqtl.csv", encoding="gb18030", low_memory=False)

E = set(env["eQTL"].astype(str))
C = set(cis["SNP"].astype(str))
T = set(trans["SNP"].astype(str))
U = C | T
out = {}

sec("0. Panel")
print(f"environment loci {len(E):,}; cis {len(C):,}; trans {len(T):,}; "
      f"cis&trans {len(C & T):,}; eQTL panel {len(U):,}; env outside panel {len(E - U):,}")
out["panel"] = dict(n_env=len(E), n_cis=len(C), n_trans=len(T), n_both=len(C & T),
                    n_panel=len(U), n_env_outside_panel=len(E - U))

ch, po = split_locus(sorted(U))
P = pd.DataFrame({"locus": sorted(U), "chrom": ch, "pos": po}).dropna()
P["pos"] = P["pos"].astype(int)
P = P.sort_values(["chrom", "pos"]).reset_index(drop=True)
P["env"] = P["locus"].isin(E)
P["cis"] = P["locus"].isin(C)
P["trans"] = P["locus"].isin(T)
P["gwas"] = P["locus"].isin(set(gwas["SNP"].astype(str)))
P["gwas_eqtl"] = P["locus"].isin(set(ge["eQTL loci"].astype(str)))

# ------------------------------------------------------------ N1 composition
sec("N1. Composition of environment-associated loci within the eQTL panel")


def fisher_row(df, flag, partner):
    a = int((df[flag] & df[partner]).sum())
    b = int((df[flag] & ~df[partner]).sum())
    c = int((~df[flag] & df[partner]).sum())
    d = int((~df[flag] & ~df[partner]).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    exp = (a + b) * (a + c) / (a + b + c + d)
    return dict(observed=a, expected=exp, ratio=a / exp, frac_target=a / (a + b),
                frac_other=c / (c + d), OR=orr, p=p, n11=a, n10=b, n01=c, n00=d)


rng = np.random.default_rng(2026)
chrom_idx = {c: g.index.values for c, g in P.groupby("chrom")}


def rotation_null(df, flag, partner, n=N_ROT):
    """Circularly shift the target labels along each chromosome."""
    lab = df[flag].values
    par = df[partner].values
    obs = int((lab & par).sum())
    null = np.empty(n)
    for k in range(n):
        tot = 0
        for _, ix in chrom_idx.items():
            s = rng.integers(1, len(ix))
            tot += int((np.roll(lab[ix], s) & par[ix]).sum())
        null[k] = tot
    hi = (np.sum(null >= obs) + 1) / (n + 1)
    lo = (np.sum(null <= obs) + 1) / (n + 1)
    return obs, float(null.mean()), float(null.std()), float(min(1.0, 2 * min(hi, lo)))


rows = []
for partner in ("cis", "trans", "gwas", "gwas_eqtl"):
    r = fisher_row(P, "env", partner)
    obs, mu, sd, prot = rotation_null(P, "env", partner)
    r.update(partner=partner, rot_mean=mu, rot_sd=sd, rot_ratio=obs / mu,
             rot_z=(obs - mu) / sd, rot_p=prot)
    rows.append(r)

# 5-kb blocks: single-linkage clustering of panel loci
blk = np.zeros(len(P), dtype=int)
bid = 0
for c, ix in chrom_idx.items():
    pos = P.loc[ix, "pos"].values
    b = np.r_[0, np.cumsum(np.diff(pos) >= BLOCK_GAP)] + bid
    blk[ix] = b
    bid = b.max() + 1
P["block"] = blk
B = P.groupby("block").agg(env=("env", "any"), cis=("cis", "any"), trans=("trans", "any"),
                           gwas=("gwas", "any"), gwas_eqtl=("gwas_eqtl", "any"))
print(f"5-kb blocks: {len(B):,} (environment blocks {int(B['env'].sum()):,})")
for r in rows:
    rb = fisher_row(B, "env", r["partner"])
    r.update(block_observed=rb["observed"], block_expected=rb["expected"],
             block_ratio=rb["ratio"], block_OR=rb["OR"], block_p=rb["p"])
comp = pd.DataFrame(rows)
comp["p_bh"] = bh(comp["p"].values)
comp["block_p_bh"] = bh(comp["block_p"].values)
cols = ["partner", "observed", "expected", "ratio", "frac_target", "frac_other", "OR", "p",
        "p_bh", "rot_mean", "rot_ratio", "rot_z", "rot_p", "block_observed", "block_expected",
        "block_ratio", "block_OR", "block_p_bh"]
print(comp[cols].to_string(index=False))
comp[cols].to_csv(TAB / "null_locus_composition.csv", index=False)
out["composition"] = comp[cols].to_dict(orient="records")
out["n_blocks"] = int(len(B))
out["n_env_blocks"] = int(B["env"].sum())

cls = np.where(P["cis"] & P["trans"], "cis+trans", np.where(P["cis"], "cis only", "trans only"))
tab3 = pd.crosstab(P["env"], cls)
chi = stats.chi2_contingency(tab3.values)
print("\nregulatory class by environment status:\n", tab3.to_string())
print(f"chi-square = {chi[0]:.1f}, df = {chi[2]}, p = {chi[1]:.2e}")
out["class_table"] = {str(k): {str(c): int(v) for c, v in r.items()} for k, r in tab3.iterrows()}
out["class_chi2"] = dict(stat=float(chi[0]), df=int(chi[2]), p=float(chi[1]))

# ------------------------------------------------------------ N3 clustering
sec("N3. Clustering with an identical statistic for observed and null sets")
envP = P[P["env"]][["chrom", "pos"]]
obs_f, obs_m = nn_stats(envP)
panel_f, panel_m = nn_stats(P[["chrom", "pos"]])
nf, nm = [], []
sizes = envP.groupby("chrom").size()
for _ in range(N_SUB):
    draw = []
    for c, n in sizes.items():
        pos = P.loc[chrom_idx[c], "pos"].values
        draw.append(pd.DataFrame({"chrom": c, "pos": rng.choice(pos, size=n, replace=False)}))
    f, m = nn_stats(pd.concat(draw))
    nf.append(f)
    nm.append(m)
nf, nm = np.array(nf), np.array(nm)
clus = dict(observed_frac_lt1kb=obs_f, observed_median_nn=obs_m,
            null_frac_mean=float(nf.mean()), null_frac_sd=float(nf.std()),
            null_median_nn_mean=float(nm.mean()), null_median_nn_sd=float(nm.std()),
            z=float((obs_f - nf.mean()) / nf.std()),
            p_emp=float((np.sum(nf >= obs_f) + 1) / (N_SUB + 1)),
            panel_frac_lt1kb=panel_f, panel_median_nn=panel_m, n_draws=N_SUB)
print(json.dumps(clus, indent=1))
out["clustering"] = clus
pd.DataFrame({"null_frac_lt1kb": nf, "null_median_nn": nm}).to_csv(
    TAB / "null_clustering_draws.csv", index=False)

# ------------------------------------------------------- N4 effect by class
sec("N4. Effect size and breadth by regulatory class")
env["abs_beta"] = pd.to_numeric(env["Beta"], errors="coerce").abs()
L = env.groupby("eQTL").agg(n_env=("Environment", "nunique"), n_rows=("Environment", "size"),
                            max_abs_beta=("abs_beta", "max"), mean_abs_beta=("abs_beta", "mean"))
L["cls"] = np.where(L.index.isin(C) & L.index.isin(T), "cis+trans",
                    np.where(L.index.isin(C), "cis only", "trans only"))
summ = L.groupby("cls").agg(n=("n_env", "size"), median_breadth=("n_env", "median"),
                            mean_breadth=("n_env", "mean"),
                            median_max=("max_abs_beta", "median"),
                            median_mean=("mean_abs_beta", "median"),
                            mean_mean=("mean_abs_beta", "mean"))
print(summ.to_string())
summ.to_csv(TAB / "null_regulatory_class.csv")
is_cis = L.index.isin(C)
mw_max = stats.mannwhitneyu(L.loc[~is_cis, "max_abs_beta"], L.loc[is_cis, "max_abs_beta"])
mw_mean = stats.mannwhitneyu(L.loc[~is_cis, "mean_abs_beta"], L.loc[is_cis, "mean_abs_beta"])
kw_b = stats.kruskal(*[g["n_env"].values for _, g in L.groupby("cls")])
kw_mean = stats.kruskal(*[g["mean_abs_beta"].values for _, g in L.groupby("cls")])
out["regulatory_class"] = dict(
    table=summ.reset_index().to_dict(orient="records"),
    n_rows_equals_n_env=bool((L["n_rows"] == L["n_env"]).all()),
    breadth_median=float(L["n_env"].median()), breadth_max=int(L["n_env"].max()),
    noncis_median_max=float(L.loc[~is_cis, "max_abs_beta"].median()),
    cis_median_max=float(L.loc[is_cis, "max_abs_beta"].median()),
    noncis_median_mean=float(L.loc[~is_cis, "mean_abs_beta"].median()),
    cis_median_mean=float(L.loc[is_cis, "mean_abs_beta"].median()),
    mw_max_p=float(mw_max.pvalue), mw_mean_p=float(mw_mean.pvalue),
    kw_breadth_H=float(kw_b.statistic), kw_breadth_p=float(kw_b.pvalue),
    kw_mean_H=float(kw_mean.statistic), kw_mean_p=float(kw_mean.pvalue))
print(json.dumps({k: v for k, v in out["regulatory_class"].items() if k != "table"}, indent=1))

# -------------------------------------------- N5 gene level, measured frame
sec("N5. Gene-level proximity restricted to directly measured gene positions")
GF = pd.read_csv(RES / "gene_frame.csv")
obsG = GF[GF["src"] == "observed"].dropna(subset=["pos_med"]).copy()
env_genes_obs = set()
dist_obs = []
for c, g in envP.groupby("chrom"):
    sub = obsG[obsG["gene_chr"] == c].sort_values("pos_med")
    if sub.empty:
        continue
    gp = sub["pos_med"].values
    j = np.searchsorted(gp, g["pos"].values)
    for x, k in zip(g["pos"].values, j):
        cand = [i for i in (k - 1, k) if 0 <= i < len(gp)]
        i = min(cand, key=lambda i: abs(gp[i] - x))
        dist_obs.append(abs(gp[i] - x))
        if abs(gp[i] - x) <= WINDOW:
            env_genes_obs.add(sub["Gene"].iloc[i])
dist_obs = np.array(dist_obs)
em = pd.read_csv(DATA / "emqtl.csv", usecols=["Gene"])
twas = pd.read_csv(DATA / "twas.csv", usecols=["Phenotype", "Gene", "Pearson", "Study"])
Ug = set(obsG["Gene"])
gsets = {"environment": env_genes_obs, "cis_eGene": set(cis["Gene"].astype(str)),
         "trans_eGene": set(trans["Gene"].astype(str)), "methylation": set(em["Gene"].astype(str)),
         "trait_TWAS": set(twas["Gene"].astype(str))}
grows = []
for a, b in itertools.combinations(gsets, 2):
    Aa, Bb = gsets[a] & Ug, gsets[b] & Ug
    n11 = len(Aa & Bb)
    n10, n01 = len(Aa) - n11, len(Bb) - n11
    n00 = len(Ug) - n11 - n10 - n01
    orr, p = stats.fisher_exact([[n11, n10], [n01, n00]])
    grows.append(dict(layer_a=a, layer_b=b, n_both=n11, n_a=len(Aa), n_b=len(Bb), OR=orr, p=p))
gq = pd.DataFrame(grows)
gq["p_bh"] = bh(gq["p"].values)
gq = gq.sort_values("p")
print(f"measured-coordinate genes: {len(Ug):,}; environment genes within 10 kb: "
      f"{len(env_genes_obs):,}; median locus-gene distance {np.median(dist_obs):,.0f} bp; "
      f"{(dist_obs <= WINDOW).mean():.1%} of loci within 10 kb")
print(gq.to_string(index=False))
gq.to_csv(TAB / "null_gene_measured_frame.csv", index=False)
membership = pd.DataFrame({k: [g in v for g in sorted(Ug)] for k, v in gsets.items()},
                          index=pd.Index(sorted(Ug), name="Gene"))
membership["n_layers"] = membership.sum(axis=1)
membership.to_csv(TAB / "null_gene_layers.csv")
out["gene_layer_counts"] = membership["n_layers"].value_counts().sort_index().to_dict()
print("genes by number of layers:", out["gene_layer_counts"])
out["gene_measured"] = dict(n_genes=len(Ug), n_env_genes=len(env_genes_obs),
                            median_dist=float(np.median(dist_obs)),
                            frac_within=float((dist_obs <= WINDOW).mean()),
                            set_sizes={k: len(v & Ug) for k, v in gsets.items()},
                            tests=gq.to_dict(orient="records"))
n_imp = int((GF["src"] == "imputed").sum())
out["gene_frame"] = dict(n_frame=int(len(GF)), n_observed=int(len(obsG)), n_imputed=n_imp)

# ------------------------------------------------ N6 normalised hit rates
sec("N6. Hit rates normalised by the number of distinct features per layer")
ew = pd.read_csv(DATA / "ewas.csv", usecols=["Methylation"])
n_sites = int(ew["Methylation"].nunique())
n_genes = int(twas["Gene"].nunique())
del ew
TP = pd.read_csv(TAB / "trait_layer_panel.csv")
TP["hr_meth"] = TP["ewas_sig"] / n_sites
TP["hr_expr"] = TP["twas_sig"] / n_genes
w = stats.wilcoxon(TP["hr_meth"], TP["hr_expr"])
n_meth_higher = int((TP["hr_meth"] > TP["hr_expr"]).sum())
ratio = np.log2(TP["hr_meth"] / TP["hr_expr"])
print(f"distinct methylation sites {n_sites:,}; distinct genes {n_genes:,}")
print(f"traits with higher methylation hit rate: {n_meth_higher}/{len(TP)}; "
      f"Wilcoxon W = {w.statistic:.0f}, p = {w.pvalue:.2e}; "
      f"median log2 ratio = {ratio.median():+.2f}")
print("by class:\n", TP.assign(l2=ratio).groupby("trait_class")["l2"].agg(["size", "median"]).to_string())
same = bool((TP["ewas_n"] == TP["ewas_sig"]).all() and (TP["twas_n"] == TP["twas_sig"]).all())
cohorts = TP.groupby("study").size().sort_values(ascending=False)
out["hit_rate"] = dict(n_sites=n_sites, n_genes=n_genes, n_traits=int(len(TP)),
                       n_meth_higher=n_meth_higher, W=float(w.statistic), p=float(w.pvalue),
                       median_log2_ratio=float(ratio.median()),
                       by_class=TP.assign(l2=ratio).groupby("trait_class")["l2"].median().to_dict(),
                       rows_equal_significant=same, cohorts=cohorts.to_dict())
TP[["Phenotype", "study", "trait_class", "ewas_sig", "twas_sig", "hr_meth", "hr_expr"]].to_csv(
    TAB / "null_trait_hit_rate.csv", index=False)

# ------------------------------------------------------ N7 descriptive
sec("N7. GWAS-eQTL join and three-layer panel")
ge["env"] = ge["eQTL loci"].astype(str).isin(E)
ld_env = ge.loc[ge["env"], "LD"]
ld_non = ge.loc[~ge["env"], "LD"]
mw_ld = stats.mannwhitneyu(ld_env, ld_non)
out["gwas_eqtl"] = dict(n_pairs=int(len(ge)), min_ld=float(ge["LD"].min()),
                        n_eqtl_loci=int(ge["eQTL loci"].nunique()),
                        n_env_eqtl_loci=int(ge.loc[ge["env"], "eQTL loci"].nunique()),
                        n_gwas_loci=int(ge["GWAS loci"].nunique()),
                        gwas_loci_in_gwas=int(ge["GWAS loci"].astype(str).isin(
                            set(gwas["SNP"].astype(str))).sum() == len(ge)),
                        median_ld_env=float(ld_env.median()),
                        median_ld_non=float(ld_non.median()), mw_p=float(mw_ld.pvalue),
                        n_env_pairs=int(len(ld_env)), n_non_pairs=int(len(ld_non)))
print(json.dumps(out["gwas_eqtl"], indent=1))
TT = pd.read_csv(TAB / "conv_trait_panel.csv")
sp = {}
for a, b in (("med_meth_absr", "med_gene_absr"), ("med_meth_absr", "med_path_absr"),
             ("med_gene_absr", "med_path_absr"), ("n_meth_sig", "n_gene_sig"),
             ("n_meth_sig", "n_pathway_sig"), ("n_gene_sig", "n_pathway_sig")):
    r = stats.spearmanr(TT[a], TT[b])
    sp[f"{a}|{b}"] = dict(rho=float(r.statistic), p=float(r.pvalue))
out["three_layer"] = dict(n=int(len(TT)), spearman=sp)
print(json.dumps(sp, indent=1))

# ------------------------------------- N8 reporting floor and effect sizes
sec("N8. Correlation floor and threshold-adjusted effect sizes")
ewr = pd.read_csv(DATA / "ewas.csv", usecols=["Phenotype", "Pearson", "Study"])
floor = {}
for nm, df in (("EWAS", ewr), ("TWAS", twas)):
    a = df["Pearson"].abs()
    by_study = df.assign(a=a).groupby("Study")["a"].min()
    floor[nm] = dict(min_abs_r=float(a.min()), frac_lt_030=float((a < 0.2995).mean()),
                     n_studies=int(len(by_study)),
                     n_studies_floor_030=int(((by_study >= 0.2995) & (by_study < 0.3005)).sum()))
print(json.dumps(floor, indent=1))
rows = []
for _, r in TP.iterrows():
    st = r["study"].split(",")
    e = ewr.loc[(ewr["Phenotype"] == r["Phenotype"]) & ewr["Study"].isin(st), "Pearson"].abs()
    t = twas.loc[(twas["Phenotype"] == r["Phenotype"]) & twas["Study"].isin(st), "Pearson"].abs()
    rows.append(dict(Phenotype=r["Phenotype"], study=r["study"], trait_class=r["trait_class"],
                     meth_min=e.min(), expr_min=t.min(), meth_med=e.median(),
                     expr_med=t.median()))
X = pd.DataFrame(rows)
X["meth_excess"] = X["meth_med"] - X["meth_min"]
X["expr_excess"] = X["expr_med"] - X["expr_min"]
w_min = stats.wilcoxon(X["meth_min"], X["expr_min"])
w_exc = stats.wilcoxon(X["meth_excess"], X["expr_excess"])
by = X.assign(larger=X["meth_excess"] > X["expr_excess"],
              d=X["meth_excess"] - X["expr_excess"]).groupby("trait_class").agg(
    n=("larger", "size"), meth_larger=("larger", "sum"), median_diff=("d", "median"))
print(by.to_string())
out["threshold"] = dict(
    floor=floor,
    min_meth_gt_expr=int((X["meth_min"] > X["expr_min"]).sum()),
    w_min=float(w_min.statistic), p_min=float(w_min.pvalue),
    excess_meth_larger=int((X["meth_excess"] > X["expr_excess"]).sum()),
    w_excess=float(w_exc.statistic), p_excess=float(w_exc.pvalue),
    rho_min=float(stats.spearmanr(X["meth_min"], X["expr_min"]).statistic),
    rho_median=float(stats.spearmanr(X["meth_med"], X["expr_med"]).statistic),
    rho_excess=float(stats.spearmanr(X["meth_excess"], X["expr_excess"]).statistic),
    p_rho_excess=float(stats.spearmanr(X["meth_excess"], X["expr_excess"]).pvalue),
    by_class=by.reset_index().to_dict(orient="records"))
print(json.dumps({k: v for k, v in out["threshold"].items() if k != "by_class"}, indent=1))
X.to_csv(TAB / "null_trait_threshold.csv", index=False)

json.dump(out, open(RES / "null_models.json", "w"), indent=2, default=float)
print("\nsaved tables/null_*.csv and results/null_models.json")