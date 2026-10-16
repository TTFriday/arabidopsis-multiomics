"""
13_convergence.py -- Analysis 3: cross-layer convergence.

Joins the layers that share identifiers or coordinates:
  * locus level  : environment-eQTL loci <-> cis/trans eQTL loci <-> GWAS loci
  * gene level   : cis/trans eGene symbols <-> EMQTL genes <-> TWAS genes,
                   plus positional gene assignment for environment-associated
                   loci using the recovered coordinate frame
  * trait level  : shared phenotype panels across EWAS / TWAS / pathway layers

Reports enrichment against chromosome-matched random locus sets and against
gene-set permutations, and assembles a per-gene multi-layer matrix.

Outputs: tables/conv_*.csv, results/convergence.json
"""
from __future__ import annotations

import itertools
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
RNG = np.random.default_rng(23)


def sec(t):
    print("\n" + "=" * 84 + f"\n{t}\n" + "=" * 84, flush=True)


env = A.load("environmenteqtl")
cis = A.load("ciseqtl")
trans = A.load("transeqtl")
em = A.load("emqtl")
twas = A.load("twas")
ewas_cols = ["Phenotype", "Phenotype ontology", "Methylation", "Pearson", "FDR_Pearson", "Study"]
ewas = A.load("ewas", usecols=ewas_cols + ["FDR_Spearman", "Spearman"])
gwe = A.load("gwaseqtl")
gwas = A.load("gwas")
GF = A.build_gene_frame()

cis_loci = set(cis["SNP"].astype(str))
trans_loci = set(trans["SNP"].astype(str))
env_loci = set(env["eQTL"].astype(str))
gwe_env_loci = set(gwe["eQTL loci"].astype(str))
gwe_gwas_loci = set(gwe["GWAS loci"].astype(str))
gwas_loci = set(gwas["SNP"].astype(str))
universe = sorted(env_loci | cis_loci | trans_loci)
pool = pd.DataFrame({"locus": universe})
pch, ppo = A.split_locus(pool["locus"])
pool["chrom"], pool["pos"] = pch, ppo
pool_by_chrom = {c: np.sort(g["pos"].values) for c, g in pool.dropna().groupby("chrom")}

sec("1. Identifier inventory")
print(f"environment-eQTL loci           : {len(env_loci):,}")
print(f"cis-eQTL loci / trans-eQTL loci : {len(cis_loci):,} / {len(trans_loci):,}")
print(f"GWAS-eQTL join loci (env side)  : {len(gwe_env_loci):,}  "
      f"(overlap with environment-eQTL loci: {len(gwe_env_loci & env_loci):,})")
print(f"GWAS-eQTL join loci (GWAS side) : {len(gwe_gwas_loci):,}  "
      f"(overlap with gwas.csv SNPs: {len(gwe_gwas_loci & gwas_loci):,})")
print(f"assayed locus universe          : {len(universe):,}")

sec("2. Environment-associated loci versus regulatory and association layers")
rows = []


def matched_test(name, target, partner, partner_name, n_perm=2000, seed=13):
    obs, null = A.matched_locus_null(target, partner, pool_by_chrom, n_perm=n_perm, seed=seed)
    z, p = A.z_empirical(obs, null)
    rec = dict(target=name, partner=partner_name, n_target=len(target), observed=obs,
               expected=float(null.mean()), sd=float(null.std()), z=z, p_emp=p,
               ratio=obs / null.mean() if null.mean() else np.nan)
    rows.append(rec)
    return rec


matched_test("environment", env_loci, cis_loci, "cis_eQTL")
matched_test("environment", env_loci, trans_loci, "trans_eQTL")
matched_test("environment", env_loci, gwe_gwas_loci, "GWAS_eQTL_partner")
matched_test("environment", env_loci, gwas_loci, "GWAS_association")
print(pd.DataFrame(rows).to_string(index=False))

sec("3. Are environment-associated loci a biased subset of cis loci? (continuous)")
locus = pd.DataFrame({"locus": sorted(env_loci)})
locus["is_cis"] = locus["locus"].isin(cis_loci)
locus["is_trans"] = locus["locus"].isin(trans_loci)
locus["n_env"] = locus["locus"].map(env.groupby("eQTL")["Environment"].nunique())
locus["max_abs_beta"] = locus["locus"].map(env.groupby("eQTL")["Beta"].apply(lambda s: s.abs().max()))
# compare environment-associated cis loci to cis loci that are not environment-associated
comp = {}
for partner, pset, pname in ((cis_loci, cis_loci, "cis"), (trans_loci, trans_loci, "trans")):
    inside = locus.loc[locus[f"is_{pname}"], "max_abs_beta"].dropna()
    outside_env = set(env_loci) - pset
    outside = locus.loc[locus["locus"].isin(outside_env), "max_abs_beta"].dropna()
    if len(outside) > 5:
        u = stats.mannwhitneyu(inside, outside)
        comp[pname] = dict(n_in=len(inside), n_out=len(outside),
                           med_in=float(inside.median()), med_out=float(outside.median()),
                           p=float(u.pvalue))
        print(f"max|beta| of environment loci inside vs outside the {pname} set: "
              f"median {inside.median():.3f} (n={len(inside)}) vs "
              f"{outside.median():.3f} (n={len(outside)}), p={u.pvalue:.3e}")

sec("4. Trait-level convergence across EWAS / TWAS / pathway layers")
pe = set(ewas["Phenotype"].unique())
pt = set(twas["Phenotype"].unique())
pp = set(A.load("phenotypepathway")["Phenotype"].unique())
core = sorted(pe & pt & pp)
print(f"EWAS {len(pe)} | TWAS {len(pt)} | pathway {len(pp)} -> three-layer panel {len(core)}")
trait_tbl = []
ppath = A.load("phenotypepathway")
for ph in core:
    e = ewas[ewas["Phenotype"] == ph]
    t = twas[twas["Phenotype"] == ph]
    p = ppath[ppath["Phenotype"] == ph]
    trait_tbl.append(dict(
        Phenotype=ph, trait_class=A.trait_class(ph),
        ewas_study=",".join(sorted(e["Study"].unique())),
        twas_study=",".join(sorted(t["Study"].unique())),
        n_meth_sig=len(e), n_gene_sig=len(t), n_pathway_sig=len(p),
        n_meth_sites=e["Methylation"].nunique(), n_genes=t["Gene"].nunique(),
        n_pathways=p["Pathway"].nunique(),
        med_meth_absr=float(e["Pearson"].abs().median()),
        med_gene_absr=float(t["Pearson"].abs().median()),
        med_path_absr=float(p["Pearson"].abs().median()),
        max_meth_absr=float(e["Pearson"].abs().max()),
        max_gene_absr=float(t["Pearson"].abs().max()),
        max_path_absr=float(p["Pearson"].abs().max()),
    ))
TT = pd.DataFrame(trait_tbl)
TT["layers_with_evidence"] = (TT[["n_meth_sig", "n_gene_sig", "n_pathway_sig"]] > 0).sum(axis=1)
print(TT.groupby("trait_class")[["n_meth_sig", "n_gene_sig", "n_pathway_sig"]].median().to_string())
print("\nSpearman across the three-layer panel:")
for a, b in (("n_meth_sig", "n_gene_sig"), ("n_meth_sig", "n_pathway_sig"),
             ("n_gene_sig", "n_pathway_sig"), ("med_meth_absr", "med_gene_absr"),
             ("med_meth_absr", "med_path_absr"), ("med_gene_absr", "med_path_absr")):
    r = stats.spearmanr(TT[a], TT[b])
    print(f"  {a} vs {b}: rho = {r.statistic:+.3f} (p={r.pvalue:.3e})")
TT.to_csv(A.TAB / "conv_trait_panel.csv", index=False)

sec("5. Gene-level layer membership")
# positional assignment of environment-associated loci to the nearest gene
gpos = GF.dropna(subset=["pos_med"])
idx = {}
for c, sub in gpos.groupby("gene_chr"):
    sub = sub.sort_values("pos_med")
    idx[c] = (sub["pos_med"].values, sub["Gene"].tolist())
trees = {c: cKDTree(v[0].reshape(-1, 1)) for c, v in idx.items()}
ech, epo = A.split_locus(env.drop_duplicates("eQTL")["eQTL"])
envpos = pd.DataFrame({"chrom": ech, "pos": epo}).dropna()
env_genes, env_dist = set(), []
for c, g in envpos.groupby("chrom"):
    if c not in trees:
        continue
    d, i = trees[c].query(g["pos"].values.reshape(-1, 1), k=1)
    env_dist.append(pd.Series(d))
    for k, dd in zip(i, d):
        if dd <= 10000:
            env_genes.add(idx[c][1][int(k)])
env_dist = pd.concat(env_dist, ignore_index=True) if env_dist else pd.Series(dtype=float)
print(f"environment-associated loci assigned to a gene within 10 kb: {len(env_genes):,} genes")
print(f"  (median distance to nearest gene {env_dist.median():,.0f} bp; "
      f"{(env_dist <= 10000).mean():.1%} of loci within 10 kb)")

sets = {
    "environment": env_genes,
    "cis_eGene": set(cis["Gene"].astype(str)),
    "trans_eGene": set(trans["Gene"].astype(str)),
    "methylation": set(em["Gene"].astype(str)),
    "trait_TWAS": set(twas["Gene"].astype(str)),
}
U = set(GF["Gene"]) | set().union(*sets.values())
print("\ngene set sizes (within the frame universe):")
for k, v in sets.items():
    print(f"  {k:12s} {len(v & U):6,}")
print(f"  universe      {len(U):6,}")

grows = []
for a, b in itertools.combinations(sets, 2):
    Aa, Bb = sets[a] & U, sets[b] & U
    n11 = len(Aa & Bb)
    n10 = len(Aa) - n11
    n01 = len(Bb) - n11
    n00 = len(U) - n11 - n10 - n01
    orr, p = stats.fisher_exact([[n11, n10], [n01, n00]])
    grows.append(dict(layer_a=a, layer_b=b, n_both=n11, n_a=len(Aa), n_b=len(Bb),
                      OR=orr, p=p))
gene_enr = pd.DataFrame(grows).sort_values("p")
gene_enr["p_bh"] = A.bh_fdr(gene_enr["p"].values)
print("\ngene-level pairwise enrichment:")
print(gene_enr.to_string(index=False))
gene_enr.to_csv(A.TAB / "conv_gene_enrichment.csv", index=False)

sec("6. Per-gene multi-layer matrix")
env_by_gene = envpos.copy()
nears = []
for c, g in envpos.groupby("chrom"):
    if c not in trees:
        nears.append(pd.Series([None] * len(g)))
        continue
    d, i = trees[c].query(g["pos"].values.reshape(-1, 1), k=1)
    nears.append(pd.Series([idx[c][1][int(k)] if dd <= 10000 else None for k, dd in zip(i, d)]))
envpos["near_gene"] = pd.concat(nears, ignore_index=True)
envbeta = env.groupby("eQTL")["Beta"].apply(lambda s: float(s.abs().max()))
envpos = envpos.merge(env.groupby("eQTL")["Environment"].nunique().rename("n_env"),
                      left_on="chrom", right_on="eQTL", how="left") if False else envpos
locus_meta = env.groupby("eQTL").agg(n_env=("Environment", "nunique"),
                                     max_abs_beta=("Beta", lambda s: float(s.abs().max())))
locus_meta.index.name = "locus"
env_locus_tbl = env.drop_duplicates("eQTL")[["eQTL"]].copy()
env_locus_tbl = env_locus_tbl.merge(locus_meta, left_on="eQTL", right_index=True, how="left")
env_locus_tbl = env_locus_tbl.merge(envpos[["chrom", "pos", "near_gene"]],
                                    left_on="eQTL", right_on="chrom", how="left") if False else env_locus_tbl
env_locus_tbl["chrom"] = A.split_locus(env_locus_tbl["eQTL"])[0]
env_locus_tbl["pos"] = A.split_locus(env_locus_tbl["eQTL"])[1]
env_locus_tbl = env_locus_tbl.merge(
    pd.DataFrame({"chrom": envpos["chrom"].values, "pos": envpos["pos"].values,
                  "near_gene": envpos["near_gene"].values}),
    on=["chrom", "pos"], how="left")
env_gene_agg = (env_locus_tbl.dropna(subset=["near_gene"])
                .groupby("near_gene").agg(n_env_locus=("eQTL", "nunique"),
                                          max_env_beta=("max_abs_beta", "max"),
                                          max_breadth=("n_env", "max")))

cis_by_gene = cis.groupby("Gene").agg(n_cis=("SNP", "nunique"),
                                      cis_absbeta_max=("Beta", lambda s: float(s.abs().max())),
                                      cis_t_max=("t_stat" if "t_stat" in cis.columns else "t-stat",
                                                 lambda s: float(s.abs().max())))
trans_by_gene = trans.groupby("Gene").agg(n_trans=("SNP", "nunique"),
                                          trans_absbeta_max=("Beta", lambda s: float(s.abs().max())),
                                          trans_t_max=("t_stat" if "t_stat" in trans.columns else "t-stat",
                                                       lambda s: float(s.abs().max())))
em_by_gene = em.groupby("Gene").agg(n_meth=("Methylation", "nunique"),
                                    meth_absr_max=("Pearson", lambda s: float(s.abs().max())),
                                    meth_absr_med=("Pearson", lambda s: float(s.abs().median())))
tw_by_gene = twas.groupby("Gene").agg(n_trait=("Phenotype", "nunique"),
                                      twas_absr_max=("Pearson", lambda s: float(s.abs().max())))
M = (env_gene_agg.join(cis_by_gene, how="outer").join(trans_by_gene, how="outer")
     .join(em_by_gene, how="outer").join(tw_by_gene, how="outer"))
counts = ["n_env_locus", "n_cis", "n_trans", "n_meth", "n_trait"]
M[counts] = M[counts].fillna(0)
M["n_layers"] = (M[counts] > 0).sum(axis=1)
print("genes by number of layers with evidence:")
print(M["n_layers"].value_counts().sort_index().to_string())
print("\nlayer co-occurrence across genes:")
lay = {"environment": M["n_env_locus"] > 0, "cis": M["n_cis"] > 0, "trans": M["n_trans"] > 0,
       "methylation": M["n_meth"] > 0, "trait": M["n_trait"] > 0}
recs = []
for a, b in itertools.combinations(lay, 2):
    tb = pd.crosstab(lay[a], lay[b])
    orr, p = stats.fisher_exact(tb.values)
    recs.append(dict(layer_a=a, layer_b=b, OR=orr, p=p, n_both=int(tb.values[1, 1]),
                     n_a=int(tb.values[1].sum()), n_b=int(tb.values[:, 1].sum())))
lc = pd.DataFrame(recs).sort_values("p")
lc["p_bh"] = A.bh_fdr(lc["p"].values)
print(lc.to_string(index=False))
M.to_csv(A.TAB / "conv_gene_multilayer.csv")
lc.to_csv(A.TAB / "conv_layer_cooccurrence.csv", index=False)

sec("7. Methylation load versus regulatory effect size")
j = em_by_gene.join(cis_by_gene, how="inner")
print(f"genes with both EMQTL and cis-eQTL data: {len(j)}")
print("Spearman(n_meth, n_cis)      =", round(float(stats.spearmanr(j['n_meth'], j['n_cis']).statistic), 4))
print("Spearman(meth|r|max, cis|t|) =", round(float(stats.spearmanr(j['meth_absr_max'], j['cis_t_max']).statistic), 4))
jt = em_by_gene.join(trans_by_gene, how="inner")
print(f"genes with both EMQTL and trans-eQTL data: {len(jt)}")
print("Spearman(n_meth, n_trans)    =", round(float(stats.spearmanr(jt['n_meth'], jt['n_trans']).statistic), 4))
jc = cis_by_gene.join(trans_by_gene, how="inner")
print(f"genes with both cis and trans eQTL data: {len(jc)}")
print("Fisher OR of cis/trans co-occurrence at gene level:",
      stats.fisher_exact(pd.crosstab(jc["n_cis"] > 0, jc["n_trans"] > 0).values))

sec("8. Environment-associated genes: do they carry stronger effects?")
for partner, col, lab in (("cis", "cis_t_max", "cis |t|"), ("trans", "trans_t_max", "trans |t|"),
                          ("methylation", "meth_absr_max", "max methylation |r|"),
                          ("trait", "n_trait", "number of associated traits")):
    sub = M[M[counts[{"cis": 1, "trans": 2, "methylation": 3, "trait": 4}[partner]]] > 0]
    inx = sub.loc[sub["n_env_locus"] > 0, col]
    out = sub.loc[sub["n_env_locus"] == 0, col]
    if len(inx) > 3 and len(out) > 3:
        u = stats.mannwhitneyu(inx, out)
        print(f"{lab}: with environment n={len(inx)} median={inx.median():.3f} | "
              f"without n={len(out)} median={out.median():.3f} | p={u.pvalue:.3e}")

sec("9. GWAS-eQTL join: trait-associated and environmentally variable loci")
ge = gwe.copy()
ge["same_chrom"] = (A.split_locus(ge["eQTL loci"])[0] == A.split_locus(ge["GWAS loci"])[0])
print(f"GWAS-eQTL pairs: {len(ge):,}")
print(f"  LD >= 0.5 in all pairs (minimum {ge['LD'].min():.3f})")
print("  Type x environment-association of the eQTL locus:")
ge["env_assoc"] = ge["eQTL loci"].isin(env_loci)
print(pd.crosstab(ge["Type"], ge["env_assoc"]).to_string())
print("\n  LD by type:")
print(ge.groupby("Type")["LD"].agg(["size", "median", "mean"]).to_string())
print("\n  LD: environment-associated vs not (Mann-Whitney):")
u = stats.mannwhitneyu(ge.loc[ge["env_assoc"], "LD"], ge.loc[~ge["env_assoc"], "LD"])
print(f"    median {ge.loc[ge['env_assoc'], 'LD'].median():.3f} vs "
      f"{ge.loc[~ge['env_assoc'], 'LD'].median():.3f}, p={u.pvalue:.3e}")
print("\n  GWAS phenotype classes among environment-associated pairs:")
print(ge.loc[ge["env_assoc"], "Phenotype ontology"].value_counts().head(10).to_string())

json.dump({
    "locus_matched_tests": rows,
    "n_env_loci": len(env_loci), "n_cis_loci": len(cis_loci),
    "n_trans_loci": len(trans_loci), "n_universe": len(universe),
    "n_env_genes_within10kb": len(env_genes),
    "median_env_locus_gene_distance": float(env_dist.median()),
    "frac_env_loci_within10kb": float((env_dist <= 10000).mean()),
    "gene_universe": len(U),
    "gene_set_sizes": {k: len(v & U) for k, v in sets.items()},
    "gene_enrichment": gene_enr.to_dict(orient="records"),
    "layer_cooccurrence": lc.to_dict(orient="records"),
    "n_genes_by_layer_count": {str(k): int(v) for k, v in M["n_layers"].value_counts().items()},
    "three_layer_trait_panel": int(len(core)),
    "methylation_vs_cis_gene": {
        "n": int(len(j)),
        "spearman_n_meth_n_cis": float(stats.spearmanr(j["n_meth"], j["n_cis"]).statistic),
        "spearman_effect": float(stats.spearmanr(j["meth_absr_max"], j["cis_t_max"]).statistic),
    },
}, open(A.RES / "convergence.json", "w"), indent=2, default=str)
print("\nsaved tables/conv_*.csv and results/convergence.json")
