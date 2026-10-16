"""
01_overlap.py -- Cross-layer integration feasibility audit.

Quantifies every meaningful key overlap between the 10 released association
tables. The result determines which multi-omics research questions are
actually answerable with this data drop, and which are not.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

METH_RE = re.compile(r"^(\d+)_(\d+)_(CG|CHG|CHH)$")
LOCUS_RE = re.compile(r"^(chr[1-5])_(\d+)$")


def read_csv(name, enc="utf-8", **kw):
    return pd.read_csv(DATA / name, encoding=enc, low_memory=False, **kw)


def parse_locus(s: pd.Series) -> pd.DataFrame:
    ex = s.astype(str).str.extract(LOCUS_RE)
    return pd.DataFrame({"chrom": ex[0], "pos": pd.to_numeric(ex[1], errors="coerce")})


def jaccard(a: set, b: set) -> float:
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def jacc(a: pd.Series) -> set:
    return set(a.dropna().astype(str).unique())


report: dict = {}
lines: list[str] = ["# Cross-layer integration audit\n"]

# ---------------------------------------------------------------- methylation
ewas = read_csv("ewas.csv")
emqtl = read_csv("emqtl.csv")
epath = read_csv("epathway.csv")

ewas_meth = jacc(ewas["Methylation"])
emqtl_meth = jacc(emqtl["Methylation"])
epath_meth = jacc(epath["Methylation"])

# methylation feature space per study (EWAS has a Study column; subset it)
ewas_by_study = {k: jacc(v) for k, v in ewas.groupby("Study")["Methylation"]}


def meth_chrom_context(sites: set) -> dict:
    m = pd.Series(sorted(sites)).str.extract(METH_RE)
    return {"n": len(sites), "n_parseable": int(m[0].notna().sum()),
            "chrom_counts": m[0].value_counts().sort_index().to_dict(),
            "context_counts": m[2].value_counts().to_dict()}


report["methylation_sites"] = {
    "ewas_unique": len(ewas_meth),
    "emqtl_unique": len(emqtl_meth),
    "epathway_unique": len(epath_meth),
    "ewas_and_emqtl": len(ewas_meth & emqtl_meth),
    "ewas_and_epathway": len(ewas_meth & epath_meth),
    "emqtl_and_epathway": len(emqtl_meth & epath_meth),
    "all_three": len(ewas_meth & emqtl_meth & epath_meth),
    "jaccard_ewas_epathway": jaccard(ewas_meth, epath_meth),
}
report["methylation_site_annotation"] = meth_chrom_context(
    ewas_meth | emqtl_meth | epath_meth)

# ---------------------------------------------------------------- phenotype axes
twas = read_csv("twas.csv")
ppath = read_csv("phenotypepathway.csv")

twas_phen = jacc(twas["Phenotype"])
ewas_phen = jacc(ewas["Phenotype"])
ppath_phen = jacc(ppath["Phenotype"])
gwas_phen = jacc(read_csv("gwas.csv")["Phenotype"])
ge_phen = jacc(read_csv("gwaseqtl.csv", enc="gb18030")["GWAS phenotype"])
env_vars = jacc(pd.read_excel(DATA / "environmenteqtl.xlsx")["Environment"])

report["phenotype_overlap"] = {
    "n_twas": len(twas_phen), "n_ewas": len(ewas_phen), "n_ppathway": len(ppath_phen),
    "n_gwas": len(gwas_phen), "n_gwaseqtl": len(ge_phen),
    "twas_and_ewas": len(twas_phen & ewas_phen),
    "twas_and_ppathway": len(twas_phen & ppath_phen),
    "ewas_and_ppathway": len(ewas_phen & ppath_phen),
    "twas_ewas_ppathway": len(twas_phen & ewas_phen & ppath_phen),
    "twas_and_gwas": len(twas_phen & gwas_phen),
    "twas_and_gwaseqtl": len(twas_phen & ge_phen),
    "ppathway_and_gwas": len(ppath_phen & gwas_phen),
    "env_in_gwas": len(env_vars & gwas_phen),
    "env_in_twas": len(env_vars & twas_phen),
    "env_in_ewas": len(env_vars & ewas_phen),
    "env_in_gwaseqtl": len(env_vars & ge_phen),
    "env_vars": sorted(env_vars),
    "example_twas_ewas": sorted(twas_phen & ewas_phen)[:15],
    "example_twas_ewas_ppath": sorted(twas_phen & ewas_phen & ppath_phen)[:15],
}

# --------------------------------------------------------------- gene universe
gene_sets = {
    "emqtl": jacc(emqtl["Gene"]),
    "ciseqtl": jacc(pd.read_excel(DATA / "ciseqtl.xls")["Gene"]),
    "transeqtl": jacc(pd.read_excel(DATA / "transeqtl.xls")["Gene"]),
    "twas": jacc(twas["Gene"]),
    "gwas_annot": jacc(read_csv("gwas.csv")["Genename"]),
}
report["gene_overlap"] = {
    "sizes": {k: len(v) for k, v in gene_sets.items()},
    "emqtl_and_ciseqtl": len(gene_sets["emqtl"] & gene_sets["ciseqtl"]),
    "emqtl_and_twas": len(gene_sets["emqtl"] & gene_sets["twas"]),
    "twas_and_ciseqtl": len(gene_sets["twas"] & gene_sets["ciseqtl"]),
    "twas_and_transeqtl": len(gene_sets["twas"] & gene_sets["transeqtl"]),
    "emqtl_and_twas_and_ciseqtl": len(gene_sets["emqtl"] & gene_sets["twas"] & gene_sets["ciseqtl"]),
    "gwas_annot_and_twas": len(gene_sets["gwas_annot"] & gene_sets["twas"]),
    "gwas_annot_and_emqtl": len(gene_sets["gwas_annot"] & gene_sets["emqtl"]),
}

# ------------------------------------------------------------------- pathways
report["pathway_overlap"] = {
    "n_epathway": len(jacc(epath["Pathway"])),
    "n_ppathway": len(jacc(ppath["Pathway"])),
    "shared": len(jacc(epath["Pathway"]) & jacc(ppath["Pathway"])),
    "shared_examples": sorted(jacc(epath["Pathway"]) & jacc(ppath["Pathway"]))[:10],
}

# ------------------------------------------------------- genomic position joins
env_eqtl = pd.read_excel(DATA / "environmenteqtl.xlsx")
ciseqtl = pd.read_excel(DATA / "ciseqtl.xls")
transeqtl = pd.read_excel(DATA / "transeqtl.xls")
gwaseqtl = read_csv("gwaseqtl.csv", enc="gb18030")
gwas = read_csv("gwas.csv")

env_loc = parse_locus(env_eqtl["eQTL"])
cis_loc = parse_locus(ciseqtl["SNP"])
trans_loc = parse_locus(transeqtl["SNP"])
ge_env_loc = parse_locus(gwaseqtl["eQTL loci"])
ge_gwas_loc = parse_locus(gwaseqtl["GWAS loci"])
gwas_loc = parse_locus(gwas["SNP"])

report["locus_parse"] = {
    "environmenteqtl_eQTL": int(env_loc["pos"].notna().sum()),
    "ciseqtl_SNP": int(cis_loc["pos"].notna().sum()),
    "transeqtl_SNP": int(trans_loc["pos"].notna().sum()),
    "gwaseqtl_eQTL": int(ge_env_loc["pos"].notna().sum()),
    "gwaseqtl_GWAS": int(ge_gwas_loc["pos"].notna().sum()),
    "gwas_SNP": int(gwas_loc["pos"].notna().sum()),
    "gwas_SNP_total": int(len(gwas_loc)),
}

env_keys = set(env_eqtl["eQTL"].astype(str))
ge_env_keys = set(gwaseqtl["eQTL loci"].astype(str))
report["eqtl_key_overlap"] = {
    "environmenteqtl_x_gwaseqtl_eQTL": len(env_keys & ge_env_keys),
    "n_environmenteqtl": len(env_keys),
    "n_gwaseqtl_eqlt": len(ge_env_keys),
}

# how many GWAS loci in gwaseqtl also appear in the GWAS association table
report["gwaseqtl_x_gwas_loci"] = len(set(gwaseqtl["GWAS loci"].astype(str)) & set(gwas["SNP"].astype(str)))

# ---------------------------------------------- provenance: study identifiers
def study_dist(df, col="Study"):
    return {str(k): int(v) for k, v in df[col].astype(str).value_counts().items()} if col in df else {}


report["studies"] = {
    "ewas": study_dist(ewas),
    "twas": study_dist(twas),
    "ppathway": study_dist(ppath),
    "gwas_numeric_study_n": int(gwas["Study"].nunique()),
}

(OUT / "overlap_audit.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

# ------------------------------------------------------------------ markdown
L = lines
L.append("\n## Methylation feature space\n")
L.append("| set | n sites |")
L.append("|---|---|")
for k, v in report["methylation_sites"].items():
    L.append(f"| {k} | {v:,} |" if isinstance(v, int) else f"| {k} | {v} |")
L.append("\n### methylation site chrom/context distribution\n")
mc = report["methylation_site_annotation"]
L.append(f"- parseable as `chrom_pos_CONTEXT`: {mc['n_parseable']:,} / {mc['n']:,}")
L.append(f"- chrom counts: {mc['chrom_counts']}")
L.append(f"- context counts: {mc['context_counts']}")

L.append("\n## Phenotype axis overlap\n")
L.append("| key | value |")
L.append("|---|---|")
for k, v in report["phenotype_overlap"].items():
    L.append(f"| {k} | {v if not isinstance(v, list) else ', '.join(map(str, v[:15]))} |")

L.append("\n## Gene overlap\n")
L.append("| key | value |")
L.append("|---|---|")
for k, v in report["gene_overlap"].items():
    L.append(f"| {k} | {v} |")

L.append("\n## Pathway overlap\n")
L.append(f"```\n{json.dumps(report['pathway_overlap'], indent=2)}\n```")

L.append("\n## Locus parsing and positional joins\n")
L.append("| key | value |")
L.append("|---|---|")
for k, v in report["locus_parse"].items():
    L.append(f"| {k} | {v:,} |")
for k, v in report["eqtl_key_overlap"].items():
    L.append(f"| {k} | {v:,} |")
L.append(f"| gwaseqtl GWAS loci also in gwas.csv | {report['gwaseqtl_x_gwas_loci']:,} |")

L.append("\n## Study provenance\n")
L.append(f"```\n{json.dumps(report['studies'], indent=2)}\n```")

(OUT / "overlap_audit.md").write_text("\n".join(L), encoding="utf-8")
print("\n".join(L))
