"""
08_enrich.py -- Gene-level and locus-level enrichment, with a block-aware
permutation null.

Because environment-eQTL loci are tightly clustered (many adjacent SNPs within
1 kb), ordinary Fisher tests treat physically linked SNPs as independent and
inflate significance. Every locus-level test is therefore repeated on a
linkage-pruned locus set and against a block-resampling null.
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
LOC = re.compile(r"^(chr[1-5])_(\d+)$")
RNG = np.random.default_rng(20260529)


def sec(t):
    print("\n" + "=" * 80 + f"\n{t}\n" + "=" * 80, flush=True)


def rd(n, enc="utf-8", **kw):
    return pd.read_csv(DATA / n, encoding=enc, low_memory=False, **kw)


def split_locus(s):
    ex = s.astype(str).str.extract(LOC)
    return ex[0].values, pd.to_numeric(ex[1], errors="coerce").values


cis = pd.read_excel(DATA / "ciseqtl.xls")
tr = pd.read_excel(DATA / "transeqtl.xls")
ee = pd.read_excel(DATA / "environmenteqtl.xlsx")
em = rd("emqtl.csv")
twas = rd("twas.csv")
gw = rd("gwas.csv")
gene_frame = pd.read_csv(OUT / "gene_frame.csv")

env_loci = set(ee["eQTL"].astype(str))
cis_loci = set(cis["SNP"].astype(str))
tr_loci = set(tr["SNP"].astype(str))
universe = env_loci | cis_loci | tr_loci

sec("A. Re-check the cis/trans locus overlap arithmetic")
print("|env|=%d |cis|=%d |trans|=%d |union|=%d" %
      (len(env_loci), len(cis_loci), len(tr_loci), len(universe)))
print("env&cis=%d env&trans=%d cis&trans=%d env&cis&trans=%d" %
      (len(env_loci & cis_loci), len(env_loci & tr_loci), len(cis_loci & tr_loci),
       len(env_loci & cis_loci & tr_loci)))
tb = {
    "a_env&cis": len(env_loci & cis_loci),
    "b_env_only": len(env_loci - cis_loci),
    "c_cis_only": len(cis_loci - env_loci),
    "d_neither_in_union": len(universe - env_loci - cis_loci),
}
print("table floor:", tb, "-> sum =", sum(tb.values()), "(union size =", len(universe), ")")

sec("B. Linkage pruning of the environment-eQTL locus set")
ee2 = ee.assign(chrom=split_locus(ee["eQTL"])[0], pos=split_locus(ee["eQTL"])[1])
pruned = []
for c, sub in ee2.groupby("chrom"):
    p = np.sort(sub["pos"].dropna().unique())
    last = -np.inf
    for x in p:
        if x - last > 5000:  # 5 kb spacing
            pruned.append(f"{c}_{int(x)}")
            last = x
pruned = set(pruned)
print(f"environment-eQTL loci: {len(env_loci):,} -> pruned (5 kb): {len(pruned):,}")

sec("C. Locus-level enrichment: observed vs block-aware null")
def fisher_orr(a_set, b_set, U):
    a = len(a_set & b_set & U)
    b = len(a_set & U) - a
    c = len(b_set & U) - a
    d = len(U) - a - b - c
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return orr, p, a, b, c, d


rows = []
for label, A in (("environment", env_loci), ("environment-pruned", pruned)):
    for bname, B in (("cis", cis_loci), ("trans", tr_loci)):
        orr, p, a, b, c, d = fisher_orr(A, B, universe)
        rows.append(dict(locus_set=label, partner=bname, odds_ratio=orr, p=p,
                         a=a, b=b, c=c, d=d))
lr = pd.DataFrame(rows)
print(lr.to_string(index=False))

# block null: shift the environment locus set along each chromosome by a random
# offset, preserving clustering, and recount the overlap
print("\nBlock-shift permutation null (200 shifts per chromosome):")


def shifted_overlap(Aset, Bset, n_shift=200):
    A = pd.DataFrame({"locus": sorted(Aset)})
    ch, po = split_locus(A["locus"])
    A["chrom"], A["pos"] = ch, po
    B = pd.DataFrame({"locus": sorted(Bset)})
    bch, bpo = split_locus(B["locus"])
    B["chrom"], B["pos"] = bch, bpo
    sizes = {c: int(g["chrom"].max() if False else 0) for c, g in A.groupby("chrom")}
    span = ee2.groupby("chrom")["pos"].max().to_dict()
    obs = len(Aset & Bset)
    null = []
    for _ in range(n_shift):
        sh = {c: RNG.integers(0, max(int(span.get(c, 1)), 1)) for c in A["chrom"].unique()}
        cand = set()
        for c, g in A.groupby("chrom"):
            newp = (g["pos"].values + sh[c]) % max(int(span.get(c, 1)) + 1, 1)
            cand.update(f"{c}_{int(x)}" for x in newp)
        null.append(len(cand & Bset))
    null = np.array(null)
    z = (obs - null.mean()) / (null.std() + 1e-12)
    pemp = (np.sum(null >= obs) + 1) / (len(null) + 1)
    return obs, null.mean(), null.std(), z, pemp


for bname, Bset in (("cis", cis_loci), ("trans", tr_loci)):
    obs, mu, sd, z, pemp = shifted_overlap(env_loci, Bset)
    print(f"  env vs {bname}: observed={obs:,} null_mean={mu:,.0f} null_sd={sd:,.0f} "
          f"z={z:.1f} p_emp={pemp:.4f}")

sec("D. Gene-level environment association via the recovered coordinate frame")
# map every environment-associated SNP to the gene in whose body it sits
gw["Effect"] = gw["Effect"].astype(str).str.upper()
ann = gw[(gw["Effect"] != "INTERGENIC") & gw["Genename"].notna()][["SNP", "Genename"]]
ann = ann[ann["SNP"].isin(env_loci)]
env_genes = set(ann["Genename"].astype(str))
# every gene covered by the functional annotation = the fair universe
covered_genes = set(gw.loc[gw["Genename"].notna(), "Genename"].astype(str))
print(f"environment-associated SNPs with an intragenic annotation: {ann['SNP'].nunique():,}")
print(f"implied environment-associated genes: {len(env_genes):,} "
      f"(of {len(covered_genes):,} annotated genes)")

cis_gene_loci = cis[cis["SNP"].isin(env_loci)]
cis_env_genes = set(cis_gene_loci["Gene"].astype(str))
print(f"eGenes reached by environment-associated cis-eQTL SNPs: {len(cis_env_genes):,}")

all_assayed = (set(cis["Gene"].astype(str)) | set(tr["Gene"].astype(str))
               | set(em["Gene"].astype(str)) | set(twas["Gene"].astype(str)))
U = covered_genes & all_assayed
print(f"joint gene universe (annotated AND assayed): {len(U):,}")


def fisher_g(name, A, B, U):
    a = len(A & B & U)
    b = len(A & U) - a
    c = len(B & U) - a
    d = len(U) - a - b - c
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return dict(test=name, a=a, b=b, c=c, d=d, odds_ratio=orr, p=p)


sets = {
    "env_assoc": env_genes,
    "env_cis_eGene": cis_env_genes,
    "cis_eGene": set(cis["Gene"].astype(str)),
    "trans_eGene": set(tr["Gene"].astype(str)),
    "emqtl": set(em["Gene"].astype(str)),
    "twas": set(twas["Gene"].astype(str)),
}
grows = []
import itertools
for a, b in itertools.combinations(sets, 2):
    grows.append(fisher_g(f"{a} vs {b}", sets[a], sets[b], U))
gr = pd.DataFrame(grows)
gr["p_bonf"] = np.minimum(gr["p"] * len(gr), 1)
gr = gr.sort_values("p")
print(gr.to_string(index=False))

sec("E. Per-gene multi-layer matrix")
env_by_gene = ann.groupby("Genename").agg(n_env_snp=("SNP", "nunique"))
envbeta = ee.groupby("eQTL")["Beta"].apply(lambda s: float(s.abs().max()))
env_by_gene["env_absbeta_max"] = (ann.merge(envbeta.rename("ab"), left_on="SNP",
                                            right_index=True, how="left")
                                  .groupby("Genename")["ab"].max())
cis_by_gene = cis.groupby("Gene").agg(n_cis=("SNP", "nunique"),
                                      cis_t_max=("t-stat", lambda s: float(s.abs().max())),
                                      cis_beta_max=("Beta", lambda s: float(s.abs().max())))
tr_by_gene = tr.groupby("Gene").agg(n_trans=("SNP", "nunique"),
                                    trans_t_max=("t-stat", lambda s: float(s.abs().max())),
                                    trans_beta_max=("Beta", lambda s: float(s.abs().max())))
em_by_gene = em.groupby("Gene").agg(n_meth=("Methylation", "nunique"),
                                    meth_absr_max=("Pearson", lambda s: float(s.abs().max())),
                                    meth_absr_med=("Pearson", lambda s: float(s.abs().median())))
tw_by_gene = twas.groupby("Gene").agg(n_trait=("Phenotype", "nunique"),
                                      twas_absr_max=("Pearson", lambda s: float(s.abs().max())),
                                      twas_fdr_min=("FDR_Pearson", "min"))
M = (env_by_gene.join(cis_by_gene, how="outer").join(tr_by_gene, how="outer")
     .join(em_by_gene, how="outer").join(tw_by_gene, how="outer"))
count_cols = ["n_env_snp", "n_cis", "n_trans", "n_meth", "n_trait"]
for c in M.columns:
    if c not in count_cols:
        M[c] = M[c]
M[count_cols] = M[count_cols].fillna(0)
M["n_layers"] = (M[count_cols] > 0).sum(axis=1)
print("genes per number of layers present:")
print(M["n_layers"].value_counts().sort_index().to_string())
print("\nlayer pair co-occurrence (Fisher, gene universe = genes in >=1 layer):")
lay = {"environment": M["n_env_snp"] > 0, "cis": M["n_cis"] > 0, "trans": M["n_trans"] > 0,
       "methylation": M["n_meth"] > 0, "trait": M["n_trait"] > 0}
GM = M.index
recs = []
for a, b in itertools.combinations(lay, 2):
    tb = pd.crosstab(lay[a], lay[b])
    orr, p = stats.fisher_exact(tb.values)
    recs.append(dict(layer_a=a, layer_b=b, OR=orr, p=p,
                     n_both=int(tb.values[1, 1]),
                     n_a=int(tb.values[1].sum()), n_b=int(tb.values[:, 1].sum())))
lc = pd.DataFrame(recs).sort_values("p")
lc["p_bonf"] = np.minimum(lc["p"] * len(lc), 1)
print(lc.to_string(index=False))

sec("F. Do environment-associated genes carry stronger regulatory effects?")
sub = M[(M["n_cis"] > 0)]
A = sub.loc[sub["n_env_snp"] > 0, "cis_t_max"]
B = sub.loc[sub["n_env_snp"] == 0, "cis_t_max"]
print(f"cis-eGenes with vs without environmental association:")
print(f"  with   : n={len(A):4d} median |t| = {A.median():.3f}")
print(f"  without: n={len(B):4d} median |t| = {B.median():.3f}")
print("  Mann-Whitney p =", f"{stats.mannwhitneyu(A, B).pvalue:.3e}")
sub2 = M[M["n_meth"] > 0]
A2 = sub2.loc[sub2["n_env_snp"] > 0, "meth_absr_max"]
B2 = sub2.loc[sub2["n_env_snp"] == 0, "meth_absr_max"]
print(f"\nEMQTL genes with vs without environmental association:")
print(f"  with   : n={len(A2):4d} median max|r| = {A2.median():.3f}")
print(f"  without: n={len(B2):4d} median max|r| = {B2.median():.3f}")
print("  Mann-Whitney p =", f"{stats.mannwhitneyu(A2, B2).pvalue:.3e}")
sub3 = M[M["n_trait"] > 0]
A3 = sub3.loc[sub3["n_env_snp"] > 0, "n_trait"]
B3 = sub3.loc[sub3["n_env_snp"] == 0, "n_trait"]
print(f"\nTWAS genes, number of associated traits:")
print(f"  with   : n={len(A3):4d} median = {A3.median():.1f}")
print(f"  without: n={len(B3):4d} median = {B3.median():.1f}")
print("  Mann-Whitney p =", f"{stats.mannwhitneyu(A3, B3).pvalue:.3e}")

lr.to_csv(OUT / "locus_enrichment.csv", index=False)
gr.to_csv(OUT / "gene_enrichment.csv", index=False)
M.to_csv(OUT / "gene_multilayer_matrix.csv")
lc.to_csv(OUT / "layer_cooccurrence.csv", index=False)
pruned_df = pd.DataFrame({"locus": sorted(pruned)})
pruned_df.to_csv(OUT / "env_loci_pruned_5kb.csv", index=False)
json.dump({"n_env_loci": len(env_loci), "n_env_loci_pruned": len(pruned),
           "n_cis_loci": len(cis_loci), "n_trans_loci": len(tr_loci),
           "n_universe": len(universe), "n_env_genes": len(env_genes),
           "n_cis_env_egenes": len(cis_env_genes), "n_joint_gene_universe": len(U)},
          open(OUT / "enrich_counts.json", "w"), indent=2)
print("\nsaved results/locus_enrichment.csv, gene_enrichment.csv, "
      "gene_multilayer_matrix.csv, layer_cooccurrence.csv")
