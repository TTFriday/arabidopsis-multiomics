"""
09_validate.py -- Resolve the cis-label anomaly and rebuild the tests properly.

Issues addressed
  1. The released cis eQTL pairs have a median SNP-eGene distance of ~0.5 Mb,
     which is not a canonical cis window. Test whether the reported eGene is the
     physically nearest gene, i.e. whether the label is a distance problem or a
     label-semantics problem.
  2. The earlier block-shift permutation was degenerate. Replace it with a
     matched random-locus null drawn from the same coordinate frame.
  3. Recompute gene-level environment association using positional proximity in
     the recovered gene frame instead of the sparse functional annotation.
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
OUT = Path(__file__).resolve().parents[1] / "results"
LOC = re.compile(r"^(chr[1-5])_(\d+)$")
RNG = np.random.default_rng(7)


def sec(t):
    print("\n" + "=" * 80 + f"\n{t}\n" + "=" * 80, flush=True)


def rd(n, enc="utf-8", **kw):
    return pd.read_csv(DATA / n, encoding=enc, low_memory=False, **kw)


def split(s):
    ex = s.astype(str).str.extract(LOC)
    return ex[0].values, pd.to_numeric(ex[1], errors="coerce").values


cis = pd.read_excel(DATA / "ciseqtl.xls")
tr = pd.read_excel(DATA / "transeqtl.xls")
ee = pd.read_excel(DATA / "environmenteqtl.xlsx")
em = rd("emqtl.csv")
twas = rd("twas.csv")
gw = rd("gwas.csv")
GF = pd.read_csv(OUT / "gene_frame.csv")
print("gene frame:", GF.shape, "| observed:", int((GF['src'] == 'observed').sum()))

# ---------------------------------------------------------------- gene index
gene_index: dict[str, tuple[np.ndarray, np.ndarray, list[str]]] = {}
for c, sub in GF.dropna(subset=["pos_med"]).groupby("gene_chr"):
    sub = sub.sort_values("pos_med")
    gene_index[c] = (sub["pos_med"].values, sub["pos_med"].values,
                     sub["Gene"].tolist())
trees = {c: cKDTree(v[0].reshape(-1, 1)) for c, v in gene_index.items()}


def nearest_gene(chrom: str, pos: float):
    if chrom not in trees:
        return None, np.inf
    d, i = trees[chrom].query([[pos]], k=1)
    return gene_index[chrom][2][int(i[0])], float(d[0])


sec("1. Is the released eGene the physically nearest gene?")
rows = []
for nm, df, gcol in (("ciseqtl", cis, "Gene"), ("transeqtl", tr, "Gene")):
    ch, po = split(df["SNP"])
    rep_ch = ("chr" + df[gcol].astype(str).str.extract(r"^AT(\d)G")[0]).values
    rec = []
    for i in range(len(df)):
        if not np.isfinite(po[i]) or rep_ch[i] not in trees:
            rec.append((np.nan, np.nan, None, None))
            continue
        # distance to the reported gene
        gpos = GF.loc[GF["Gene"] == df[gcol].iat[i], "pos_med"]
        d_rep = abs(po[i] - float(gpos.iat[0])) if len(gpos) else np.nan
        ng, d_min = nearest_gene(rep_ch[i], po[i])
        rec.append((d_rep, d_min, ng, ng == df[gcol].iat[i]))
    r = pd.DataFrame(rec, columns=["d_reported", "d_nearest", "nearest_gene", "is_nearest"])
    print(f"\n{nm}: n={len(r):,}")
    print(f"  distance to reported eGene  : median={r['d_reported'].median():>10,.0f} "
          f"p90={r['d_reported'].quantile(.9):>12,.0f}")
    print(f"  distance to nearest gene    : median={r['d_nearest'].median():>10,.0f} "
          f"p90={r['d_nearest'].quantile(.9):>12,.0f}")
    print(f"  fraction where reported eGene IS the nearest gene: "
          f"{r['is_nearest'].mean():.3f}")
    print(f"  fraction within 100 kb of any gene: {(r['d_nearest'] <= 1e5).mean():.3f}")
    print(f"  fraction within 10 kb of any gene : {(r['d_nearest'] <= 1e4).mean():.3f}")
    print(f"  Spearman(d_reported, d_nearest) = "
          f"{stats.spearmanr(r['d_reported'], r['d_nearest']).statistic:.3f}")
    rows.append((nm, r))

sec("2. Cross-check with the genome-wide SNP->gene annotation in gwas.csv")
gw["Effect"] = gw["Effect"].astype(str).str.upper()
ann = gw[(gw["Effect"] != "INTERGENIC") & gw["Genename"].notna()][["SNP", "Genename"]]
annmap = ann.groupby("SNP")["Genename"].apply(list)
print("annotated SNP count:", len(annmap))
for nm, df, gcol in (("ciseqtl", cis, "Gene"), ("transeqtl", tr, "Gene")):
    hit = df["SNP"].isin(annmap.index)
    sub = df[hit]
    agree = [df[gcol].iat[i] in annmap[df["SNP"].iat[i]] for i in np.where(hit)[0]]
    print(f"{nm}: {hit.sum():,}/{len(df):,} SNPs annotated; "
          f"reported eGene agrees with annotation in {np.mean(agree) if agree else float('nan'):.3f}")

sec("3. Matched random-locus null for locus-level enrichment")
cis_loci = set(cis["SNP"].astype(str))
tr_loci = set(tr["SNP"].astype(str))
env_loci = set(ee["eQTL"].astype(str))
env_df = pd.DataFrame({"locus": sorted(env_loci)})
env_df["chrom"], env_df["pos"] = split(env_df["locus"])
# candidate pool: every locus present anywhere in the release (same frame/coverage)
cand = sorted(set(ee["eQTL"].astype(str)) | cis_loci | tr_loci)
cand_df = pd.DataFrame({"locus": cand})
cand_df["chrom"], cand_df["pos"] = split(cand_df["locus"])
by_chrom = {c: g["pos"].values for c, g in cand_df.dropna().groupby("chrom")}

N_PERM = 2000
null_cis, null_tr = [], []
env_pairs = list(zip(env_df["chrom"].values, env_df["pos"].values))
for _ in range(N_PERM):
    pick = set()
    for c, p in env_pairs:
        pool = by_chrom.get(c)
        if pool is None:
            continue
        pick.add(f"{c}_{int(pool[RNG.integers(len(pool))])}")
    null_cis.append(len(pick & cis_loci))
    null_tr.append(len(pick & tr_loci))
obs_cis, obs_tr = len(env_loci & cis_loci), len(env_loci & tr_loci)
for nm, obs, null in (("cis", obs_cis, np.array(null_cis)), ("trans", obs_tr, np.array(null_tr))):
    z = (obs - null.mean()) / (null.std() + 1e-12)
    p = (np.sum(null >= obs) + 1) / (N_PERM + 1)
    print(f"env vs {nm}: obs={obs:,} null={null.mean():,.0f}+-{null.std():,.0f} z={z:.1f} p={p:.4f}")

sec("4. Genomic distance of environment-associated loci to their nearest gene")
ch, po = split(ee["eQTL"])
nearest = [nearest_gene(ch[i], po[i]) for i in range(len(ee))]
ee["near_gene"] = [n[0] for n in nearest]
ee["near_dist"] = [n[1] for n in nearest]
print("distance from environment-associated SNP to nearest gene, quantiles (bp):")
print(ee["near_dist"].describe(percentiles=[.1, .25, .5, .75, .9]).round(0).to_string())
for w in (1000, 5000, 10000, 50000):
    print(f"  within {w:>6,} bp of a gene: {(ee['near_dist'] <= w).mean():.3f}")
env_near_genes = set(ee.loc[ee["near_dist"] <= 10000, "near_gene"].dropna())
print("distinct genes within 10 kb of an environment-associated locus:",
      len(env_near_genes))

sec("5. Same for cis/trans eQTL SNPs (comparison of genomic placement)")
for nm, df in (("ciseqtl", cis), ("transeqtl", tr)):
    ch2, po2 = split(df["SNP"])
    d = np.array([nearest_gene(ch2[i], po2[i])[1] for i in range(len(df))])
    print(f"{nm}: median nearest-gene distance = {np.median(d):,.0f} bp | "
          f"within 10 kb = {np.mean(d <= 1e4):.3f} | within 100 kb = {np.mean(d <= 1e5):.3f}")

sec("6. Gene-level enrichment using positional gene assignment")
all_assayed = (set(cis["Gene"]) | set(tr["Gene"]) | set(em["Gene"]) | set(twas["Gene"]))
GFg = set(GF["Gene"])
U = all_assayed & GFg
print("joint gene universe:", len(U))
sets = {
    "env_near10kb": env_near_genes & U,
    "cis_eGene": set(cis["Gene"]) & U,
    "trans_eGene": set(tr["Gene"]) & U,
    "emqtl": set(em["Gene"]) & U,
    "twas": set(twas["Gene"]) & U,
}
print({k: len(v) for k, v in sets.items()})
rows = []
import itertools
for a, b in itertools.combinations(sets, 2):
    A, B = sets[a], sets[b]
    aa = len(A & B)
    bb = len(A) - aa
    cc = len(B) - aa
    dd = len(U) - aa - bb - cc
    orr, p = stats.fisher_exact([[aa, bb], [cc, dd]])
    rows.append(dict(a=a, b=b, n_both=aa, n_a=len(A), n_b=len(B), OR=orr, p=p))
er = pd.DataFrame(rows).sort_values("p")
er["p_bonf"] = np.minimum(er["p"] * len(er), 1)
print(er.to_string(index=False))

er.to_csv(OUT / "gene_enrichment_positional.csv", index=False)
ee[["eQTL", "Environment", "Beta", "p-value", "FDR", "type", "near_gene", "near_dist"]] \
    .to_csv(OUT / "env_locus_nearest_gene.csv", index=False)
json.dump({"obs_cis": obs_cis, "obs_trans": obs_tr,
           "null_cis_mean": float(np.mean(null_cis)), "null_trans_mean": float(np.mean(null_tr)),
           "n_env_genes_near10kb": len(env_near_genes),
           "n_joint_universe": len(U)},
          open(OUT / "validate_summary.json", "w"), indent=2)
print("\nsaved results/validate_summary.json, gene_enrichment_positional.csv, env_locus_nearest_gene.csv")
