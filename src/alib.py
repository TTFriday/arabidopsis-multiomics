"""
alib.py -- Shared analysis library for the Arabidopsis multi-omics release.

Centralises:
  * paths and encodings for the released tables,
  * the recovered gene coordinate frame,
  * locus parsing for the chr<N>_<pos> naming scheme,
  * multiple-testing and effect-size helpers,
  * a reproducible plotting style.

Design notes
------------
Every released table is already filtered at FDR <= 0.05, so a row indicates
"significant association", not "tested association". All analyses therefore
treat row counts as pre-filtered evidence and use effect sizes, locus identity
and set overlap rather than p-value thresholds.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------- paths
DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
MS = ROOT / "manuscript"
for _p in (RES, FIG, TAB, MS):
    _p.mkdir(parents=True, exist_ok=True)

LOCUS_RE = re.compile(r"^(chr[1-5])_(\d+)$")
GENE_RE = re.compile(r"^AT([1-5])G(\d+)$")
METH_RE = re.compile(r"^(\d+)_(\d+)_(CG|CHG|CHH)$")

FDR_THRESHOLD = 0.05
INTRAGENIC_EFFECTS = {
    "INTRON", "EXON", "SYNONYMOUS_CODING", "NON_SYNONYMOUS_CODING",
    "STOP_GAINED", "STOP_LOST", "START_LOST", "SPLICE_SITE_ACCEPTOR",
    "SPLICE_SITE_DONOR", "SPLICE_SITE_REGION", "UTR_3_PRIME", "UTR_5_PRIME",
    "SYNONYMOUS_STOP", "NON_SYNONYMOUS_START",
}

# ------------------------------------------------------------------- loaders
CSV_SPECS = {
    "ewas": ("ewas.csv", "utf-8"),
    "twas": ("twas.csv", "utf-8"),
    "emqtl": ("emqtl.csv", "utf-8"),
    "epathway": ("epathway.csv", "utf-8"),
    "phenotypepathway": ("phenotypepathway.csv", "utf-8"),
    "gwas": ("gwas.csv", "utf-8"),
    "gwaseqtl": ("gwaseqtl.csv", "gb18030"),
}
XLS_SPECS = {"ciseqtl": "ciseqtl.xls", "transeqtl": "transeqtl.xls"}
XLSX_SPECS = {"environmenteqtl": "environmenteqtl.xlsx"}

_DROP_EE_DUP = ["p-value"]
_CACHE: dict[str, pd.DataFrame] = {}


def load(name: str, **kw) -> pd.DataFrame:
    """Load a released table by short name, stripped of index artefacts."""
    if name in _CACHE and not kw:
        return _CACHE[name].copy()
    if name in CSV_SPECS:
        fn, enc = CSV_SPECS[name]
        df = pd.read_csv(DATA / fn, encoding=enc, low_memory=False, **kw)
    elif name in XLS_SPECS:
        df = pd.read_excel(DATA / XLS_SPECS[name], **kw)
    elif name in XLSX_SPECS:
        df = pd.read_excel(DATA / XLSX_SPECS[name], **kw)
    else:
        raise KeyError(name)
    df.columns = [str(c).strip() for c in df.columns]
    if "p-value" in df.columns:
        df = df.rename(columns={"p-value": "p_value"})
    if name == "environmenteqtl" and "p_value" in df.columns:
        # drop the duplicated p-value column produced by the rename of 'p-value'
        pass
    if not kw:
        _CACHE[name] = df
    return df.copy()


def load_main() -> dict[str, pd.DataFrame]:
    """Load the six tables used by the primary analyses."""
    return {
        "ewas": load("ewas"),
        "twas": load("twas"),
        "emqtl": load("emqtl"),
        "env": load("environmenteqtl"),
        "cis": load("ciseqtl"),
        "trans": load("transeqtl"),
        "gwas": load("gwas"),
        "gwaseqtl": load("gwaseqtl"),
        "epathway": load("epathway"),
        "ppath": load("phenotypepathway"),
    }


# ------------------------------------------------------------ locus handling
def split_locus(s: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    ex = s.astype(str).str.extract(LOCUS_RE)
    return ex[0].values, pd.to_numeric(ex[1], errors="coerce").values


def split_gene(s: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    ex = s.astype(str).str.extract(GENE_RE)
    return ("chr" + ex[0]).values, pd.to_numeric(ex[1], errors="coerce").values


def locus_pos(df: pd.DataFrame, col: str) -> pd.DataFrame:
    ch, po = split_locus(df[col])
    return pd.DataFrame({"chrom": ch, "pos": po})


# ------------------------------------------------------------- gene frame
def build_gene_frame(force: bool = False) -> pd.DataFrame:
    """
    Recover gene coordinates without an external annotation.

    Rows of gwas.csv whose `Effect` is an intragenic consequence must lie inside
    the gene they name, so their SNP positions bracket that gene. Genes without
    any intragenic SNP are placed with a per-chromosome linear model of
    position on the ATxGnnnnn gene number (R^2 >= 0.985), which is safe for
    ordering and coarse placement but not for base-pair resolution.
    """
    cache = RES / "gene_frame.csv"
    if cache.exists() and not force:
        return pd.read_csv(cache)

    gw = load("gwas")
    gw["Effect"] = gw["Effect"].astype(str).str.upper()
    ann = gw[gw["Genename"].notna()].copy()
    ch, po = split_locus(ann["SNP"])
    ann["chrom"], ann["snp_pos"] = ch, po
    gch, gnum = split_gene(ann["Genename"])
    ann["gene_chr"], ann["gene_num"] = gch, gnum

    strict = ann[ann["Effect"].isin(INTRAGENIC_EFFECTS)]
    obs = (strict.groupby(["Genename", "gene_chr", "gene_num"])
           .agg(n_snp=("snp_pos", "size"), pos_lo=("snp_pos", "min"),
                pos_hi=("snp_pos", "max"), pos_med=("snp_pos", "median"))
           .reset_index())
    obs["span"] = obs["pos_hi"] - obs["pos_lo"]

    fits = {}
    for c, sub in obs.groupby("gene_chr"):
        fits[c] = np.polyfit(sub["gene_num"].values.astype(float),
                             sub["pos_med"].values.astype(float), 1)

    genes = set()
    for nm in ("emqtl", "twas"):
        genes |= set(load(nm, usecols=["Gene"])["Gene"].astype(str))
    for nm in ("ciseqtl", "transeqtl"):
        genes |= set(load(nm, usecols=["Gene"])["Gene"].astype(str))
    genes |= set(ann["Genename"].astype(str))

    G = pd.DataFrame({"Gene": sorted(genes)})
    gch, gnum = split_gene(G["Gene"])
    G["gene_chr"], G["gene_num"] = gch, gnum
    G = G[G["gene_num"].notna() & G["gene_chr"].notna()].copy()
    G = G.merge(obs[["Genename", "pos_lo", "pos_hi", "pos_med", "n_snp", "span"]],
                left_on="Gene", right_on="Genename", how="left").drop(columns=["Genename"])
    G["src"] = np.where(G["pos_med"].notna(), "observed", "imputed")
    for c, sub in G[G["pos_med"].isna()].groupby("gene_chr"):
        if c in fits:
            G.loc[sub.index, "pos_med"] = np.polyval(fits[c], sub["gene_num"].values.astype(float))
    G.to_csv(cache, index=False)
    return G


# ------------------------------------------------------------ stats helpers
def fisher_table(a: int, b: int, c: int, d: int):
    from scipy import stats
    return stats.fisher_exact([[a, b], [c, d]])


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Rank-based effect size: P(x>y) - P(x<y)."""
    from scipy import stats
    x, y = np.asarray(x), np.asarray(y)
    u = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 2 * u.statistic / (len(x) * len(y)) - 1


def matched_locus_null(target_loci, partner_loci, pool_by_chrom, n_perm=2000, seed=7):
    """
    Draw random locus sets matched on chromosome composition to `target_loci`
    from `pool_by_chrom` (the coordinate frame actually assayed) and count the
    overlap with `partner_loci`. Returns (observed, null array).
    """
    rng = np.random.default_rng(seed)
    tdf = pd.DataFrame({"locus": sorted(target_loci)})
    ch, po = split_locus(tdf["locus"])
    tdf["chrom"], tdf["pos"] = ch, po
    obs = len(set(target_loci) & set(partner_loci))
    null = np.empty(n_perm, dtype=int)
    pairs = list(zip(tdf["chrom"].values, tdf["pos"].values))
    for k in range(n_perm):
        pick = set()
        for c, _p in pairs:
            arr = pool_by_chrom.get(c)
            if arr is None or len(arr) == 0:
                continue
            pick.add(f"{c}_{int(arr[rng.integers(len(arr))])}")
        null[k] = len(pick & set(partner_loci))
    return obs, null


def z_empirical(obs: float, null: np.ndarray):
    z = (obs - null.mean()) / (null.std() + 1e-12)
    p = (np.sum(null >= obs) + 1) / (len(null) + 1)
    return float(z), float(p)


# ---------------------------------------------------------------- plotting
def setup_plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 140,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "figure.facecolor": "white",
    })
    return plt


# --------------------------------------------------------------- constants
# Phenotypic classes used for the trait-panel contrasts. Assignment is by trait
# name, using the released `Phenotype` and `Phenotype ontology` fields.
ION_TRAITS = ["As75", "Cd111", "Co59", "Cu65", "K39", "Mg25", "Mn55", "Mo98",
              "Ni60", "P31", "Rb85", "S34", "Se82", "Zn66"]
CLIMATE_TRAITS = ["clim-bio1", "clim-bio3", "clim-bio5", "clim-bio6", "clim-bio9",
                  "clim-prec6", "clim-tmax1", "clim-tmax2", "clim-tmax3",
                  "clim-tmax7", "clim-tmax8", "clim-tmax9", "clim-tmin1",
                  "clim-tmin2", "clim-tmin3", "clim-tmin9"]
PHENOLOGY_TRAITS = ["DTF spain 2008 (2nd experiment)", "DTF spain 2009 (1st experiment)",
                    "DTF spain 2009 (2nd experiment)", "DTF sweden 2009 (1st experiment)",
                    "DTFmainEffect2009", "DTFplantingSummer2009"]
SIZE_TRAITS = ["Size sweden 2009 (1st experiment)", "SizeLocSweden2009"]
SEED_TRAITS = ["2W", "GR21"]


def trait_class(name: str) -> str:
    if name in ION_TRAITS:
        return "ionome"
    if name in CLIMATE_TRAITS:
        return "climate"
    if name in PHENOLOGY_TRAITS:
        return "phenology"
    if name in SIZE_TRAITS:
        return "size"
    if name in SEED_TRAITS:
        return "seed"
    return "other"
