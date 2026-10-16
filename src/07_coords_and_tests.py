"""
07_coords_and_tests.py -- Establish a common coordinate frame and run the
decisive cross-layer statistical tests.

Coordinate frame
    SNPs      : chr<1-5>_<position> as released.
    Genes     : positions recovered from the functional annotation in gwas.csv
                (rows whose Effect is an intragenic consequence must lie inside
                the gene they name), supplemented by a per-chromosome linear
                model  position ~ ATxGnnnnn  (R2 > 0.98) for genes lacking
                intragenic SNPs.

Tests
    T1  Coordinate consistency: are the released "cis" eQTL pairs local?
    T2  Environment-association vs cis-/trans-regulatory architecture
        (Fisher on the assayed gene universe).
    T3  Environmental breadth of a locus scales with its regulatory
        connectivity and effect size.
    T4  Multi-layer convergence: genes carrying environment, cis-/trans-eQTL,
        methylation and trait associations.
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
OUT.mkdir(parents=True, exist_ok=True)
LOC = re.compile(r"^(chr[1-5])_(\d+)$")
INTRAGENIC = {"INTRON", "EXON", "SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING",
              "STOP_GAINED", "STOP_LOST", "START_LOST", "SPLICE_SITE_ACCEPTOR",
              "SPLICE_SITE_DONOR", "SPLICE_SITE_REGION", "UTR_3_PRIME",
              "UTR_5_PRIME", "SYNONYMOUS_STOP", "NON_SYNONYMOUS_START"}


def sec(t):
    print("\n" + "=" * 80 + f"\n{t}\n" + "=" * 80, flush=True)


def rd(n, enc="utf-8", **kw):
    return pd.read_csv(DATA / n, encoding=enc, low_memory=False, **kw)


def split_locus(s):
    ex = s.astype(str).str.extract(LOC)
    return ex[0].values, pd.to_numeric(ex[1], errors="coerce").values


# ---------------------------------------------------------------- gene universe
gw = rd("gwas.csv")
g = gw.dropna(subset=["Genename"]).copy()
g["Effect"] = g["Effect"].astype(str).str.upper()
g["chrom"], g["snp_pos"] = split_locus(g["SNP"])
gx = g["Genename"].astype(str).str.extract(r"^AT(\d)G(\d+)$")
g["gene_num"] = pd.to_numeric(gx[1], errors="coerce")
g["gene_chr"] = "chr" + gx[0]

sec("Coordinate frame construction")
strict = g[g["Effect"].isin(INTRAGENIC)]
obs = strict.groupby(["Genename", "gene_chr", "gene_num"]).agg(
    n_snp=("snp_pos", "size"), pos_lo=("snp_pos", "min"),
    pos_hi=("snp_pos", "max"), pos_med=("snp_pos", "median")).reset_index()
obs["span"] = obs["pos_hi"] - obs["pos_lo"]
print(f"genes with an observed intragenic SNP: {len(obs):,}")
print(f"median intragenic span: {obs['span'].median():.0f} bp "
      f"(p90 = {obs['span'].quantile(.9):,.0f} bp)")

fits = {}
for c, sub in obs.groupby("gene_chr"):
    A = np.polyfit(sub["gene_num"].values.astype(float),
                   sub["pos_med"].values.astype(float), 1)
    pred = np.polyval(A, sub["gene_num"].values.astype(float))
    fits[c] = A
    print(f"  {c}: n={len(sub):5d} slope={A[0]:7.2f} R2="
          f"{1 - (sub['pos_med'].values - pred).var() / sub['pos_med'].values.var():.4f} "
          f"resid_sd={np.std(sub['pos_med'].values - pred):>9,.0f}")

# all genes mentioned anywhere in the release
gene_frames = {
    "emqtl": rd("emqtl.csv", usecols=["Gene"]),
    "twas": rd("twas.csv", usecols=["Gene"]),
    "ciseqtl": pd.read_excel(DATA / "ciseqtl.xls", usecols=["Gene"]),
    "transeqtl": pd.read_excel(DATA / "transeqtl.xls", usecols=["Gene"]),
}
all_genes = pd.unique(np.concatenate([f["Gene"].astype(str).values
                                      for f in gene_frames.values()]))
G = pd.DataFrame({"Gene": all_genes})
gx2 = G["Gene"].str.extract(r"^AT(\d)G(\d+)$")
G["gene_chr"] = "chr" + gx2[0]
G["gene_num"] = pd.to_numeric(gx2[1], errors="coerce")
G = G[G["gene_num"].notna() & G["gene_chr"].notna()].copy()
G = G.merge(obs[["Genename", "pos_lo", "pos_hi", "pos_med", "n_snp"]],
            left_on="Gene", right_on="Genename", how="left").drop(columns=["Genename"])
G["src"] = np.where(G["pos_med"].notna(), "observed", "imputed")
for c, sub in G[G["pos_med"].isna()].groupby("gene_chr"):
    if c in fits:
        G.loc[sub.index, "pos_med"] = np.polyval(fits[c], sub["gene_num"].values.astype(float))
print(f"\ngene frame: {len(G):,} genes "
      f"({int((G['src'] == 'observed').sum()):,} observed, "
      f"{int((G['src'] == 'imputed').sum()):,} imputed)")
G.to_csv(OUT / "gene_frame.csv", index=False)

gpos = G.dropna(subset=["pos_med"]).set_index("Gene")[["gene_chr", "pos_med", "src"]]

# ------------------------------------------------------------------- T1: cis?
sec("T1. Are the released 'cis' eQTL pairs local? (coordinate consistency)")
cis = pd.read_excel(DATA / "ciseqtl.xls")
tr = pd.read_excel(DATA / "transeqtl.xls")
tr.rename(columns={"p-value": "p_value"}, inplace=True)
cis.rename(columns={"p-value": "p_value"}, inplace=True)


def snp_gene_dist(df, snpcol="SNP", genecol="Gene"):
    chrom, pos = split_locus(df[snpcol])
    out = pd.DataFrame({"chrom": chrom, "snp_pos": pos})
    out = out.join(gpos, on=df[genecol].astype(str).values)
    out["dist"] = np.abs(out["snp_pos"] - out["pos_med"])
    out["same_chr"] = out["chrom"] == out["gene_chr"]
    return out


for nm, df in (("ciseqtl", cis), ("transeqtl", tr)):
    d = snp_gene_dist(df)
    m = d["same_chr"]
    print(f"\n{nm}: n={len(d):,} | same chromosome = {m.mean():.3f}")
    print("  distance to eGene, same-chromosome pairs (bp):")
    print("   ", d.loc[m, "dist"].describe(percentiles=[.25, .5, .75, .9]).round(0).to_dict())
    print(f"    within 10 kb: {(d.loc[m, 'dist'] <= 10_000).mean():.3f} | "
          f"within 100 kb: {(d.loc[m, 'dist'] <= 100_000).mean():.3f}")

sec("T1b. Same test using the gwas.csv annotation as an independent frame")
d = snp_gene_dist(cis)
print("cis pairs whose SNP is intragenic for the named gene (gwas.csv):",
      int(cis["SNP"].isin(set(g["SNP"])).sum()))

# ------------------------------------------------------- T2/T3: env vs eQTL
sec("T2. Environment-associated loci vs cis/trans eQTL -- SNP-level join")
ee = pd.read_excel(DATA / "environmenteqtl.xlsx")
ee_alt = pd.to_numeric(ee["Alteration"], errors="coerce")
ee["env_beta_abs"] = ee["Beta"].abs()
env_loci = set(ee["eQTL"].astype(str))
cis_loci = set(cis["SNP"].astype(str))
tr_loci = set(tr["SNP"].astype(str))
print(f"environment-eQTL loci : {len(env_loci):,}")
print(f"cis-eQTL loci         : {len(cis_loci):,}  overlap with env loci: "
      f"{len(env_loci & cis_loci):,} ({len(env_loci & cis_loci) / len(env_loci):.1%} of env)")
print(f"trans-eQTL loci       : {len(tr_loci):,}  overlap with env loci: "
      f"{len(env_loci & tr_loci):,} ({len(env_loci & tr_loci) / len(env_loci):.1%} of env)")
print(f"cis trans overlap     : {len(cis_loci & tr_loci):,}")

# universe = all loci tested in any of the three scans
universe = env_loci | cis_loci | tr_loci
print(f"\nunion locus universe  : {len(universe):,}")
rows = []


def fisher2(name, A, B, U):
    a = len(A & B & U)
    b = len(A & U) - a
    c = len(B & U) - a
    d = len(U) - a - b - c
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    rows.append(dict(test=name, a=a, b=b, c=c, d=d, odds_ratio=orr, p=p))


fisher2("env vs cis", env_loci, cis_loci, universe)
fisher2("env vs trans", env_loci, tr_loci, universe)
fisher2("cis vs trans", cis_loci, tr_loci, universe)
res = pd.DataFrame(rows)
res["p_bonf"] = np.minimum(res["p"] * len(res), 1)
print(res.to_string(index=False))

sec("T2b. Gene-level version of the same test")
gx_env = set(g.loc[g["SNP"].isin(env_loci), "Genename"].dropna().astype(str))
gx_any = set(g["Genename"].dropna().astype(str))
print("genes with an annotated SNP:", len(gx_any),
      "| of which environment-associated (any annotated SNP):", len(gx_env))
GU = gx_any  # gene universe = genes covered by the functional annotation
rows = []


def fisher_g(name, A, B, U):
    a = len(A & B & U)
    b = len(A & U) - a
    c = len(B & U) - a
    d = len(U) - a - b - c
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    rows.append(dict(test=name, a=a, b=b, c=c, d=d, odds_ratio=orr, p=p))


gx_cis = set(cis["Gene"].astype(str))
gx_tr = set(tr["Gene"].astype(str))
gx_em = set(gene_frames["emqtl"]["Gene"].astype(str))
gx_tw = set(gene_frames["twas"]["Gene"].astype(str))
fisher_g("env-gene vs cis-gene", gx_env, gx_cis, GU)
fisher_g("env-gene vs trans-gene", gx_env, gx_tr, GU)
fisher_g("env-gene vs emqtl-gene", gx_env, gx_em, GU)
fisher_g("env-gene vs twas-gene", gx_env, gx_tw, GU)
fisher_g("cis-gene vs trans-gene", gx_cis, gx_tr, GU)
fisher_g("cis-gene vs emqtl-gene", gx_cis, gx_em, GU)
fisher_g("trans-gene vs twas-gene", gx_tr, gx_tw, GU)
fisher_g("emqtl-gene vs twas-gene", gx_em, gx_tw, GU)
gr = pd.DataFrame(rows)
gr["p_bonf"] = np.minimum(gr["p"] * len(gr), 1)
print(gr.to_string(index=False))

sec("T3. Environmental breadth: effect size, and regulatory connectivity")
breadth = ee.groupby("eQTL").agg(n_env=("Environment", "nunique"),
                                 max_absbeta=("Beta", lambda s: float(s.abs().max())))
breadth["icis"] = breadth.index.isin(cis_loci)
breadth["itrans"] = breadth.index.isin(tr_loci)
print("breadth (n environments) describe:")
print(breadth["n_env"].describe().to_string())
print("\nSpearman(n_env, max|beta|) =",
      round(float(stats.spearmanr(breadth["n_env"], breadth["max_absbeta"]).statistic), 4),
      " p =", f"{stats.spearmanr(breadth['n_env'], breadth['max_absbeta']).pvalue:.3e}")
print("\nmax|beta| by cis/trans membership (Mann-Whitney):")
for col in ("icis", "itrans"):
    A = breadth.loc[breadth[col], "max_absbeta"]
    B = breadth.loc[~breadth[col], "max_absbeta"]
    u = stats.mannwhitneyu(A, B)
    print(f"  {col}: n_with={len(A)} median={A.median():.3f} | "
          f"n_without={len(B)} median={B.median():.3f} | p={u.pvalue:.3e}")
print("\nbreadth by cis/trans class:")
cls = np.where(breadth["icis"] & breadth["itrans"], "cis+trans",
               np.where(breadth["icis"], "cis only",
                        np.where(breadth["itrans"], "trans only", "neither")))
breadth["cls"] = cls
print(breadth.groupby("cls")["n_env"].agg(["size", "median", "mean"]).to_string())
print("\nKruskal-Wallis across classes:",
      stats.kruskal(*[v["n_env"].values for _, v in breadth.groupby("cls")]))

sec("T3b. Which environments have the strongest / most pleiotropic associations")
envstat = ee.groupby("Environment").agg(
    n_loci=("eQTL", "nunique"), mean_absbeta=("Beta", lambda s: float(s.abs().mean())),
    max_absbeta=("Beta", lambda s: float(s.abs().max())),
    med_p=("p-value", "median"))
cis_frac = (ee.assign(is_cis=ee["eQTL"].isin(cis_loci))
            .groupby("Environment")["is_cis"].mean())
envstat["cis_frac"] = cis_frac
print(envstat.sort_values("n_loci", ascending=False).to_string())
envstat.to_csv(OUT / "environment_stats.csv")

sec("T4. Multi-layer convergence per gene")
# environment association per gene (SNP within the gene body/annotation)
env_by_gene = g[g["SNP"].isin(env_loci)].groupby("Genename").agg(
    n_env_snp=("SNP", "nunique"))
env_by_gene["env_beta_max"] = g[g["SNP"].isin(env_loci)].merge(
    ee.groupby("eQTL")["Beta"].apply(lambda s: float(s.abs().max())).rename("ab"),
    left_on="SNP", right_index=True, how="left").groupby("Genename")["ab"].max()

cis_by_gene = cis.groupby("Gene").agg(n_cis=("SNP", "nunique"),
                                      cis_t_max=("t-stat", lambda s: float(s.abs().max())))
tr_by_gene = tr.groupby("Gene").agg(n_trans=("SNP", "nunique"),
                                    trans_t_max=("t-stat", lambda s: float(s.abs().max())))
em_by_gene = gene_frames["emqtl"].groupby("Gene").agg(
    n_meth=("Methylation", "nunique"),
    meth_r_max=("Pearson", lambda s: float(s.abs().max())))
tw_by_gene = gene_frames["twas"].groupby("Gene").agg(
    n_trait=("Phenotype", "nunique"),
    twas_r_max=("Pearson", lambda s: float(s.abs().max())))

M = (cis_by_gene.join(tr_by_gene, how="outer")
     .join(em_by_gene, how="outer").join(tw_by_gene, how="outer")
     .join(env_by_gene, how="outer"))
for c in M.columns:
    M[c] = M[c].fillna(0)
M["n_layers"] = (M[["n_cis", "n_trans", "n_meth", "n_trait", "n_env_snp"]] > 0).sum(axis=1)
print("genes per number of layers present:")
print(M["n_layers"].value_counts().sort_index().to_string())
print("\nlayer co-occurrence (Fisher OR / p):")
lay = {"cis": M["n_cis"] > 0, "trans": M["n_trans"] > 0, "methylation": M["n_meth"] > 0,
       "trait": M["n_trait"] > 0, "environment": M["n_env_snp"] > 0}
import itertools
recs = []
for a, b in itertools.combinations(lay, 2):
    tb = pd.crosstab(lay[a], lay[b])
    orr, p = stats.fisher_exact(tb.values)
    recs.append(dict(layer_a=a, layer_b=b, OR=orr, p=p,
                     n_both=int(tb.values[1, 1]), n_a=int(tb.values[1].sum()),
                     n_b=int(tb.values[:, 1].sum())))
lr = pd.DataFrame(recs).sort_values("p")
lr["p_bonf"] = np.minimum(lr["p"] * len(lr), 1)
print(lr.to_string(index=False))

M.to_csv(OUT / "gene_multilayer_matrix.csv")
res.to_csv(OUT / "locus_level_fisher.csv", index=False)
gr.to_csv(OUT / "gene_level_fisher.csv", index=False)
lr.to_csv(OUT / "layer_cooccurrence.csv", index=False)
breadth.to_csv(OUT / "env_locus_breadth.csv")
json.dump({"n_env_loci": len(env_loci), "n_cis_loci": len(cis_loci),
           "n_trans_loci": len(tr_loci), "n_universe": len(universe)},
          open(OUT / "locus_counts.json", "w"), indent=2)
print("\nsaved results/*.csv")
