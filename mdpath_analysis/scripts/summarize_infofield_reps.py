#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Summarize 3 replicate "Mut-WT difference infofield" CSVs.
Outputs:
  - <outprefix>.residue_summary.csv      (per-residue mean/sd/CI across replicates)
  - <outprefix>.global_per_rep.csv       (per-rep global mean across residues)
  - <outprefix>.global_summary.csv       (global mean/sd/CI across reps)
  - (optional) <outprefix>.motif_per_rep.csv
  - (optional) <outprefix>.motif_summary.csv

Also prints Top-N residues ranked by one or more summary columns to terminal.

Example:
  python3 summarize_infofield_reps.py \
    --csv rep1.csv rep2.csv rep3.csv \
    --outprefix 6FYV \
    --rank_by mean_abs_delta_div_R3.0_mean,mean_dFmag_R3.0_mean \
    --topn 10
"""

import argparse
import json
import sys
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


METRICS = [
    "delta_div_at_CA",
    "mean_abs_delta_div_R3.0",
    "dFmag_at_CA",
    "mean_dFmag_R3.0",
]


# ---- stats helpers ----

def tcrit_975(df: int) -> float:
    """Two-sided 95% CI t critical value (approx) for small df; fallback to 1.96."""
    table = {
        1: 12.706,  # n=2
        2: 4.303,   # n=3
        3: 3.182,   # n=4
        4: 2.776,   # n=5
        5: 2.571,   # n=6
        6: 2.447,   # n=7
        7: 2.365,   # n=8
        8: 2.306,   # n=9
        9: 2.262,   # n=10
    }
    return table.get(df, 1.96)


def t_ci95(mean: float, sd: float, n: int) -> Tuple[float, float]:
    """Mean 95% CI using t distribution; returns (l, u)."""
    if n <= 1 or not np.isfinite(sd):
        return (np.nan, np.nan)
    df = n - 1
    tc = tcrit_975(df)
    half = tc * sd / np.sqrt(n)
    return (mean - half, mean + half)


def bootstrap_ci(vals: np.ndarray, n_boot: int = 5000, seed: int = 0) -> Tuple[float, float]:
    """Bootstrap CI for mean; returns (l, u)."""
    vals = np.array([v for v in vals if np.isfinite(v)], dtype=float)
    if len(vals) == 0:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    boots = rng.choice(vals, size=(n_boot, len(vals)), replace=True).mean(axis=1)
    return (float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975)))


# ---- IO helpers ----

def read_one(path: str, rep_id: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "res" not in df.columns:
        raise ValueError(f"[{path}] missing required column 'res'. Columns={list(df.columns)}")

    missing = [c for c in METRICS if c not in df.columns]
    if missing:
        raise ValueError(f"[{path}] missing required columns: {missing}. Columns={list(df.columns)}")

    keep = ["res"] + METRICS
    df = df[keep].copy()
    df["res"] = df["res"].astype(int)
    df["rep"] = rep_id
    return df


# ---- summaries ----

def make_residue_summary(df_all: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for res, g in df_all.groupby("res", sort=True):
        out: Dict[str, object] = {"res": int(res), "n_rep": int(g["rep"].nunique())}
        for m in METRICS:
            vals = g[m].astype(float).values
            finite = np.isfinite(vals)
            n = int(finite.sum())
            mean = float(np.nanmean(vals)) if n > 0 else np.nan
            sd = float(np.nanstd(vals, ddof=1)) if n >= 2 else np.nan

            t_l, t_u = t_ci95(mean, sd, n)
            b_l, b_u = bootstrap_ci(vals, n_boot=5000, seed=42)

            out[f"{m}_mean"] = mean
            out[f"{m}_sd"] = sd
            out[f"{m}_tci95_l"] = t_l
            out[f"{m}_tci95_u"] = t_u
            out[f"{m}_bci95_l"] = b_l
            out[f"{m}_bci95_u"] = b_u
        rows.append(out)

    return pd.DataFrame(rows).sort_values("res").reset_index(drop=True)


def make_global_summary(df_all: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # per replicate: global mean across residues
    per_rep_rows = []
    for rep, g in df_all.groupby("rep", sort=True):
        out = {"rep": rep, "n_res": int(g["res"].nunique())}
        for m in METRICS:
            out[m] = float(np.nanmean(g[m].astype(float).values))
        per_rep_rows.append(out)
    per_rep = pd.DataFrame(per_rep_rows).sort_values("rep").reset_index(drop=True)

    # across replicates: mean/sd/CI
    sum_rows = []
    n_rep = per_rep.shape[0]
    for m in METRICS:
        vals = per_rep[m].astype(float).values
        finite = np.isfinite(vals)
        n = int(finite.sum())
        mean = float(np.nanmean(vals)) if n > 0 else np.nan
        sd = float(np.nanstd(vals, ddof=1)) if n >= 2 else np.nan
        t_l, t_u = t_ci95(mean, sd, n)
        b_l, b_u = bootstrap_ci(vals, n_boot=5000, seed=7)
        sum_rows.append({
            "metric": m, "n_rep": n_rep,
            "mean": mean, "sd": sd,
            "tci95_l": t_l, "tci95_u": t_u,
            "bci95_l": b_l, "bci95_u": b_u,
        })
    global_sum = pd.DataFrame(sum_rows)
    return per_rep, global_sum


def make_motif_summary(df_all: pd.DataFrame, motif_json_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    motif_json format:
    {
      "agg": "max" | "mean",
      "motifs": {
        "β3–VAIK(K)": [41],
        "Hinge": [9,94]
      }
    }
    """
    with open(motif_json_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    motifs: Dict[str, List[int]] = cfg["motifs"]
    agg = cfg.get("agg", "max")
    if agg not in ("max", "mean"):
        raise ValueError(f"motif_json agg must be 'max' or 'mean', got: {agg}")

    per_rep_rows = []
    for rep, g in df_all.groupby("rep", sort=True):
        res_to_row = g.set_index("res")
        for motif, res_list in motifs.items():
            out = {"rep": rep, "motif": motif, "residues": ",".join(map(str, res_list))}
            for m in METRICS:
                vals = []
                for r in res_list:
                    r = int(r)
                    if r in res_to_row.index:
                        vals.append(float(res_to_row.loc[r, m]))
                    else:
                        vals.append(np.nan)
                vals = np.array(vals, dtype=float)
                if agg == "mean":
                    out[m] = float(np.nanmean(vals))
                else:
                    out[m] = float(np.nanmax(vals))
            per_rep_rows.append(out)
    motif_per_rep = pd.DataFrame(per_rep_rows).sort_values(["motif", "rep"]).reset_index(drop=True)

    # summarize across reps for each motif
    sum_rows = []
    for motif, gg in motif_per_rep.groupby("motif", sort=True):
        for m in METRICS:
            vals = gg[m].astype(float).values
            finite = np.isfinite(vals)
            n = int(finite.sum())
            mean = float(np.nanmean(vals)) if n > 0 else np.nan
            sd = float(np.nanstd(vals, ddof=1)) if n >= 2 else np.nan
            t_l, t_u = t_ci95(mean, sd, n)
            b_l, b_u = bootstrap_ci(vals, n_boot=5000, seed=13)
            sum_rows.append({
                "motif": motif, "metric": m,
                "mean": mean, "sd": sd,
                "tci95_l": t_l, "tci95_u": t_u,
                "bci95_l": b_l, "bci95_u": b_u,
            })
    motif_sum = pd.DataFrame(sum_rows).sort_values(["motif", "metric"]).reset_index(drop=True)
    return motif_per_rep, motif_sum


# ---- printing ----

def print_top_residues(res_sum: pd.DataFrame, rank_by_list: List[str], topn: int, title: str) -> None:
    print("\n" + title)
    for rank_by in rank_by_list:
        rank_by = rank_by.strip()
        if not rank_by:
            continue
        if rank_by not in res_sum.columns:
            print(f"[WARN] rank_by column not found: {rank_by}")
            # show hints
            hints = [c for c in res_sum.columns if c.endswith("_mean")]
            print("       Available *_mean columns include:", ", ".join(hints))
            continue

        tmp = res_sum[["res", rank_by]].copy()
        tmp = tmp.replace([np.inf, -np.inf], np.nan).dropna()
        tmp = tmp.sort_values(rank_by, ascending=False).head(topn)

        print(f"\n[Top {topn}] ranked by {rank_by}")
        for i, (res, val) in enumerate(tmp.values.tolist(), start=1):
            print(f"{i:>2d}. res {int(res):>4d}   {rank_by} = {float(val):.6g}")


# ---- main ----

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs=3, required=True,
                    help="3 replicate CSV paths (Mut-WT difference infofield). Each must contain columns: res,"
                         " delta_div_at_CA, mean_abs_delta_div_R3.0, dFmag_at_CA, mean_dFmag_R3.0")
    ap.add_argument("--outprefix", required=True, help="output prefix, e.g. 6FYV")
    ap.add_argument("--motif_json", default=None, help="optional motifs json, e.g. motifs_6FYV.json")
    ap.add_argument("--rank_by", default="mean_abs_delta_div_R3.0_mean",
                    help="comma-separated columns to rank by for terminal output. "
                         "Default: mean_abs_delta_div_R3.0_mean")
    ap.add_argument("--topn", type=int, default=10, help="how many residues to print (default: 10)")
    args = ap.parse_args()

    # read
    try:
        dfs = [read_one(p, rep_id=f"rep{i}") for i, p in enumerate(args.csv, start=1)]
    except Exception as e:
        print(f"[ERROR] failed to read input CSVs: {e}", file=sys.stderr)
        sys.exit(1)

    df_all = pd.concat(dfs, ignore_index=True)

    # residue summary
    res_sum = make_residue_summary(df_all)
    res_sum_path = f"{args.outprefix}.residue_summary.csv"
    res_sum.to_csv(res_sum_path, index=False)

    # global summary
    global_per_rep, global_sum = make_global_summary(df_all)
    global_per_rep_path = f"{args.outprefix}.global_per_rep.csv"
    global_sum_path = f"{args.outprefix}.global_summary.csv"
    global_per_rep.to_csv(global_per_rep_path, index=False)
    global_sum.to_csv(global_sum_path, index=False)

    # motif summary (optional)
    motif_per_rep_path = motif_sum_path = None
    if args.motif_json:
        try:
            motif_per_rep, motif_sum = make_motif_summary(df_all, args.motif_json)
            motif_per_rep_path = f"{args.outprefix}.motif_per_rep.csv"
            motif_sum_path = f"{args.outprefix}.motif_summary.csv"
            motif_per_rep.to_csv(motif_per_rep_path, index=False)
            motif_sum.to_csv(motif_sum_path, index=False)
        except Exception as e:
            print(f"[ERROR] failed motif summary: {e}", file=sys.stderr)
            sys.exit(2)

    # terminal Top-N
    rank_by_list = [s.strip() for s in args.rank_by.split(",") if s.strip()]
    print_top_residues(
        res_sum,
        rank_by_list=rank_by_list,
        topn=args.topn,
        title=f"=== {args.outprefix}: residues with largest changes (across 3 reps) ==="
    )

    # done
    print("\n[OK] written files:")
    print(f" - {res_sum_path}")
    print(f" - {global_per_rep_path}")
    print(f" - {global_sum_path}")
    if motif_per_rep_path and motif_sum_path:
        print(f" - {motif_per_rep_path}")
        print(f" - {motif_sum_path}")


if __name__ == "__main__":
    main()
