"""
00_profile.py -- Structural profiling of the multi-omics data drop.

Reads every dataset in the data directory, reports shape, columns, dtypes,
missingness, uniqueness of key columns, and numeric ranges. Writes a machine
readable profile to results/profile.json and a human readable summary to
results/profile.md.
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

DATA = Path(r"C:\Users\TU\AppData\Local\Programs\DeepSeek Harness\data")
OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

CSV_FILES = {
    "emqtl": "emqtl.csv",
    "epathway": "epathway.csv",
    "gwas": "gwas.csv",
    "gwaseqtl": ("gwaseqtl.csv", "gb18030"),
    "phenotypepathway": "phenotypepathway.csv",
    "twas": "twas.csv",
    "ewas": "ewas.csv",
}
XLS_FILES = {"ciseqtl": "ciseqtl.xls", "transeqtl": "transeqtl.xls"}
XLSX_FILES = {"environmenteqtl": "environmenteqtl.xlsx"}


def series_profile(s: pd.Series) -> dict:
    d = {
        "dtype": str(s.dtype),
        "n_missing": int(s.isna().sum()),
        "pct_missing": round(float(s.isna().mean() * 100), 4),
        "n_unique": int(s.nunique(dropna=True)),
    }
    if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
        desc = s.describe()
        d["stats"] = {k: (None if pd.isna(v) else float(v)) for k, v in desc.items()}
    else:
        vc = s.astype(str).value_counts().head(5)
        d["top_values"] = {str(k): int(v) for k, v in vc.items()}
        d["examples"] = [str(x) for x in s.dropna().unique()[:5]]
    return d


def profile_df(name: str, df: pd.DataFrame, source: str) -> dict:
    return {
        "name": name,
        "source_file": source,
        "n_rows": int(len(df)),
        "n_cols": int(df.shape[1]),
        "columns": [str(c) for c in df.columns],
        "memory_mb": round(float(df.memory_usage(deep=True).sum()) / 1e6, 2),
        "fields": {str(c): series_profile(df[c]) for c in df.columns},
    }


def load_all() -> dict:
    frames: dict[str, pd.DataFrame] = {}
    for name, spec in CSV_FILES.items():
        if isinstance(spec, tuple):
            fn, enc = spec
        else:
            fn, enc = spec, "utf-8"
        p = DATA / fn
        print(f"[load] {fn} (enc={enc}) ...", flush=True)
        df = pd.read_csv(p, low_memory=False, encoding=enc)
        df.columns = [str(c).strip() for c in df.columns]
        frames[name] = df
        print(f"       -> {df.shape}", flush=True)
    for name, fn in XLS_FILES.items():
        p = DATA / fn
        print(f"[load] {fn} ...", flush=True)
        df = pd.read_excel(p)
        df.columns = [str(c).strip() for c in df.columns]
        frames[name] = df
        print(f"       -> {df.shape}", flush=True)
    for name, fn in XLSX_FILES.items():
        p = DATA / fn
        print(f"[load] {fn} ...", flush=True)
        xl = pd.ExcelFile(p)
        print(f"       sheets: {xl.sheet_names}", flush=True)
        df = xl.parse(xl.sheet_names[0])
        df.columns = [str(c).strip() for c in df.columns]
        frames[name] = df
        print(f"       -> {df.shape}", flush=True)
    return frames


def main() -> None:
    frames = load_all()
    profiles = {n: profile_df(n, df, (CSV_FILES.get(n)[0] if isinstance(CSV_FILES.get(n), tuple)
                                      else CSV_FILES.get(n)) or XLS_FILES.get(n) or XLSX_FILES.get(n))
                for n, df in frames.items()}

    with open(OUT / "profile.json", "w", encoding="utf-8") as fh:
        json.dump(profiles, fh, indent=2, ensure_ascii=False)

    lines = ["# Data profile\n"]
    for n, p in profiles.items():
        lines.append(f"\n## {n}  (`{p['source_file']}`)\n")
        lines.append(f"- rows: **{p['n_rows']:,}**  cols: {p['n_cols']}  (~{p['memory_mb']} MB in memory)")
        lines.append("\n| column | dtype | missing % | unique | notes |")
        lines.append("|---|---|---|---|---|")
        for c, s in p["fields"].items():
            if "stats" in s:
                st = s["stats"]
                notes = f"min={st.get('min'):.4g} med={st.get('50%'):.4g} max={st.get('max'):.4g}"
            else:
                notes = "; ".join(f"{k} ({v})" for k, v in list(s.get("top_values", {}).items())[:2])
                if len(notes) > 90:
                    notes = notes[:87] + "..."
            lines.append(f"| {c} | {s['dtype']} | {s['pct_missing']} | {s['n_unique']:,} | {notes} |")
    (OUT / "profile.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nWrote", OUT / "profile.json", "and", OUT / "profile.md")


if __name__ == "__main__":
    main()
