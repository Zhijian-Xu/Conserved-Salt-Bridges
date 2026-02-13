#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import ast
import json
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    from scipy.ndimage import gaussian_filter
except Exception:
    gaussian_filter = None


# ----------------------------
# Helpers
# ----------------------------
def undirected_key(i, j):
    return (i, j) if i < j else (j, i)


def parse_res_label(s: str) -> int:
    m = re.search(r"(\d+)", str(s))
    if not m:
        raise ValueError(f"Cannot parse residue index from: {s}")
    return int(m.group(1))


def load_nmi_map(nmi_csv: Path, shift: int = 0) -> dict:
    """
    Read nmi_df.csv. Column might be named 'MI Difference' but stores NMI in your code.
    Returns dict {(min(i,j), max(i,j)): nmi}.
    """
    df = pd.read_csv(nmi_csv)
    if "Residue Pair" not in df.columns:
        raise ValueError(f"{nmi_csv} missing column 'Residue Pair'")

    # choose value column
    val_col = None
    for c in df.columns:
        if c.lower() in ["mi difference", "nmi", "normalized_mutual_info", "normalized mutual info"]:
            val_col = c
            break
    if val_col is None:
        val_col = [c for c in df.columns if c != "Residue Pair"][0]

    pairs = df["Residue Pair"].apply(ast.literal_eval)
    ii = pairs.apply(lambda x: parse_res_label(x[0]) + shift)
    jj = pairs.apply(lambda x: parse_res_label(x[1]) + shift)
    ww = df[val_col].astype(float).values

    m = {}
    for a, b, w in zip(ii, jj, ww):
        a, b = int(a), int(b)
        k = undirected_key(a, b)
        if k in m:
            m[k] = 0.5 * (m[k] + float(w))
        else:
            m[k] = float(w)
    return m


def iter_paths(output_txt: Path, shift: int = 0):
    """
    Yield nodes list from output.txt lines:
      Path: [np.int64(8), 88, ...], Total Weight: ...
    """
    pat_npint = re.compile(r"np\.int64\((\d+)\)")
    pat = re.compile(r"^Path:\s*\[(.*)\],\s*Total Weight:\s*([0-9eE\.\+\-]+)\s*$")

    with output_txt.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            m = pat.match(line)
            if not m:
                continue
            raw = pat_npint.sub(r"\1", m.group(1))
            nodes = [int(x.strip()) + shift for x in raw.split(",") if x.strip()]
            if len(nodes) >= 2:
                yield nodes


def load_ca_coords_from_pdb(pdb_path: Path, index_mode: str = "sequential", shift: int = 0) -> dict:
    """
    Extract CA coords from PDB.
    index_mode:
      - sequential: CA residues encountered are labeled 1..N (recommended for mdpath 'Res i')
      - resseq: use PDB resSeq numbers (only if mdpath indices are actual resSeq)
    """
    coords = {}
    seen = set()
    seq_idx = 0

    with pdb_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            resseq = int(line[22:26].strip())
            icode = line[26].strip()
            key = (chain, resseq, icode)
            if key in seen:
                continue
            seen.add(key)

            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])

            if index_mode == "sequential":
                seq_idx += 1
                rid = seq_idx + shift
            elif index_mode == "resseq":
                rid = resseq + shift
            else:
                raise ValueError("index_mode must be 'sequential' or 'resseq'")

            coords[int(rid)] = np.array([x, y, z], dtype=float)

    if not coords:
        raise ValueError(f"No CA found in {pdb_path}")
    return coords


def kabsch_align(P: np.ndarray, Q: np.ndarray):
    """
    Find rotation R and translations to map Q -> P in least squares sense.
    """
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    C = Q0.T @ P0
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    return R, Pc, Qc


def apply_alignment(coords: dict, R: np.ndarray, Pc: np.ndarray, Qc: np.ndarray) -> dict:
    out = {}
    for k, v in coords.items():
        out[k] = (v - Qc) @ R + Pc
    return out


# ----------------------------
# Field build
# ----------------------------
def make_grid_from_points(points: np.ndarray, h: float, padding: float):
    mn = points.min(axis=0) - padding
    mx = points.max(axis=0) + padding
    size = mx - mn
    nx, ny, nz = (np.ceil(size / h).astype(int) + 1).tolist()
    return mn, (nx, ny, nz)


def world_to_idx(xyz: np.ndarray, mn: np.ndarray, h: float, shape):
    ijk = np.floor((xyz - mn) / h).astype(int)
    ijk = np.clip(ijk, [0, 0, 0], [shape[0] - 1, shape[1] - 1, shape[2] - 1])
    return tuple(ijk.tolist())


def build_vector_field_from_paths(coords: dict, nmi: dict, output_txt: Path,
                                  mn: np.ndarray, shape, h: float,
                                  vector_mode: str = "unit",
                                  shift_paths: int = 0):
    """
    Deposit each directed edge vector (u->v) at its midpoint voxel.
    vector_mode:
      - unit: f = NMI * unit(rv-ru)
      - raw : f = NMI * (rv-ru)
    """
    Fx = np.zeros(shape, dtype=float)
    Fy = np.zeros(shape, dtype=float)
    Fz = np.zeros(shape, dtype=float)

    used = 0
    skipped = 0

    for nodes in iter_paths(output_txt, shift=shift_paths):
        for u, v in zip(nodes[:-1], nodes[1:]):
            if u not in coords or v not in coords:
                skipped += 1
                continue
            w = nmi.get(undirected_key(u, v), None)
            if w is None:
                skipped += 1
                continue

            ru = coords[u]
            rv = coords[v]
            d = rv - ru
            norm = float(np.linalg.norm(d))
            if norm <= 1e-12:
                skipped += 1
                continue

            if vector_mode == "unit":
                vec = d / norm
            elif vector_mode == "raw":
                vec = d
            else:
                raise ValueError("vector_mode must be 'unit' or 'raw'")

            f = float(w) * vec
            mid = 0.5 * (ru + rv)
            ix, iy, iz = world_to_idx(mid, mn, h, shape)

            Fx[ix, iy, iz] += f[0]
            Fy[ix, iy, iz] += f[1]
            Fz[ix, iy, iz] += f[2]
            used += 1

    return Fx, Fy, Fz, {"used_edges": used, "skipped_edges": skipped}


def smooth_field(Fx, Fy, Fz, sigma: float):
    if sigma <= 0:
        return Fx, Fy, Fz
    if gaussian_filter is None:
        raise RuntimeError("scipy not available; install scipy or set --smooth_sigma 0")
    Fx = gaussian_filter(Fx, sigma=sigma, mode="nearest")
    Fy = gaussian_filter(Fy, sigma=sigma, mode="nearest")
    Fz = gaussian_filter(Fz, sigma=sigma, mode="nearest")
    return Fx, Fy, Fz


def divergence(Fx, Fy, Fz, h: float):
    dFx_dx = np.gradient(Fx, h, axis=0)
    dFy_dy = np.gradient(Fy, h, axis=1)
    dFz_dz = np.gradient(Fz, h, axis=2)
    return dFx_dx + dFy_dy + dFz_dz


# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Align WT/MUT PDBs, build 3D vector fields from mdpath paths with NMI-weighted edge vectors, compute divergence, compare."
    )
    ap.add_argument("--wt_dir", required=True, help="WT keyres folder containing nmi_df.csv and output.txt")
    ap.add_argument("--mut_dir", required=True, help="MUT keyres folder containing nmi_df.csv and output.txt")
    ap.add_argument("--wt_pdb", required=True, help="WT protein PDB (for CA coords)")
    ap.add_argument("--mut_pdb", required=True, help="MUT protein PDB (for CA coords)")

    ap.add_argument("--pdb_index_mode", default="sequential", choices=["sequential", "resseq"],
                    help="Map CA residues to indices. Use 'sequential' for mdpath 'Res i' most of the time.")
    ap.add_argument("--shift_mut", type=int, default=0,
                    help="If needed: WT_index = MUT_index + shift_mut. Applies to MUT coords/NMI/paths.")
    ap.add_argument("--grid_h", type=float, default=1.0, help="Voxel size in Å")
    ap.add_argument("--padding", type=float, default=4.0, help="Grid padding around WT coords in Å")
    ap.add_argument("--smooth_sigma", type=float, default=1.0,
                    help="Gaussian smoothing sigma in voxel units (0 disables).")
    ap.add_argument("--vector_mode", default="unit", choices=["unit", "raw"],
                    help="unit: NMI*unit_direction (recommended). raw: NMI*displacement")
    ap.add_argument("--out_dir", default=None, help="Output directory (default: mut_dir/DIVERGENCE_COMPARE)")

    args = ap.parse_args()

    wt_dir = Path(args.wt_dir)
    mut_dir = Path(args.mut_dir)
    out_dir = Path(args.out_dir) if args.out_dir else (mut_dir / "DIVERGENCE_COMPARE")
    out_dir.mkdir(parents=True, exist_ok=True)

    # files
    wt_nmi_csv = wt_dir / "nmi_df.csv"
    mut_nmi_csv = mut_dir / "nmi_df.csv"
    wt_out = wt_dir / "output.txt"
    mut_out = mut_dir / "output.txt"
    for p in [wt_nmi_csv, mut_nmi_csv, wt_out, mut_out, Path(args.wt_pdb), Path(args.mut_pdb)]:
        if not p.exists():
            raise FileNotFoundError(p)

    # 1) Load CA coords
    wt_coords = load_ca_coords_from_pdb(Path(args.wt_pdb), index_mode=args.pdb_index_mode, shift=0)
    mut_coords = load_ca_coords_from_pdb(Path(args.mut_pdb), index_mode=args.pdb_index_mode, shift=args.shift_mut)

    # 2) Align MUT -> WT using common residue indices
    common = sorted(set(wt_coords.keys()) & set(mut_coords.keys()))
    if len(common) < 20:
        raise RuntimeError(f"Too few common CA residues for alignment: {len(common)}. "
                           f"Check pdb_index_mode or shift_mut.")
    P = np.stack([wt_coords[i] for i in common], axis=0)
    Q = np.stack([mut_coords[i] for i in common], axis=0)
    R, Pc, Qc = kabsch_align(P, Q)
    mut_coords_aligned = apply_alignment(mut_coords, R, Pc, Qc)

    # 3) Load NMI maps
    wt_nmi = load_nmi_map(wt_nmi_csv, shift=0)
    mut_nmi = load_nmi_map(mut_nmi_csv, shift=args.shift_mut)

    # 4) Define a single common grid (based on WT coords after padding)
    wt_points = np.stack(list(wt_coords.values()), axis=0)
    mn, shape = make_grid_from_points(wt_points, h=args.grid_h, padding=args.padding)

    # 5) Build vector fields
    Fx_wt, Fy_wt, Fz_wt, info_wt = build_vector_field_from_paths(
        wt_coords, wt_nmi, wt_out, mn=mn, shape=shape, h=args.grid_h,
        vector_mode=args.vector_mode, shift_paths=0
    )
    Fx_mut, Fy_mut, Fz_mut, info_mut = build_vector_field_from_paths(
        mut_coords_aligned, mut_nmi, mut_out, mn=mn, shape=shape, h=args.grid_h,
        vector_mode=args.vector_mode, shift_paths=args.shift_mut
    )

    # 6) Smooth (optional)
    Fx_wt, Fy_wt, Fz_wt = smooth_field(Fx_wt, Fy_wt, Fz_wt, sigma=args.smooth_sigma)
    Fx_mut, Fy_mut, Fz_mut = smooth_field(Fx_mut, Fy_mut, Fz_mut, sigma=args.smooth_sigma)

    # 7) Divergence
    div_wt = divergence(Fx_wt, Fy_wt, Fz_wt, h=args.grid_h)
    div_mut = divergence(Fx_mut, Fy_mut, Fz_mut, h=args.grid_h)
    div_delta = div_mut - div_wt

    # 8) Save arrays
    np.save(out_dir / "wt_Fx.npy", Fx_wt); np.save(out_dir / "wt_Fy.npy", Fy_wt); np.save(out_dir / "wt_Fz.npy", Fz_wt)
    np.save(out_dir / "mut_Fx.npy", Fx_mut); np.save(out_dir / "mut_Fy.npy", Fy_mut); np.save(out_dir / "mut_Fz.npy", Fz_mut)
    np.save(out_dir / "wt_div.npy", div_wt)
    np.save(out_dir / "mut_div.npy", div_mut)
    np.save(out_dir / "delta_div_mut_minus_wt.npy", div_delta)

    summary = {
        "params": {
            "grid_h": args.grid_h,
            "padding": args.padding,
            "smooth_sigma": args.smooth_sigma,
            "vector_mode": args.vector_mode,
            "pdb_index_mode": args.pdb_index_mode,
            "shift_mut": args.shift_mut,
        },
        "grid": {
            "origin_mn": mn.tolist(),
            "shape": list(shape),
        },
        "wt": info_wt,
        "mut": info_mut,
        "div_stats": {
            "wt": {"min": float(div_wt.min()), "max": float(div_wt.max()), "mean": float(div_wt.mean()), "std": float(div_wt.std())},
            "mut": {"min": float(div_mut.min()), "max": float(div_mut.max()), "mean": float(div_mut.mean()), "std": float(div_mut.std())},
            "delta": {"min": float(div_delta.min()), "max": float(div_delta.max()), "mean": float(div_delta.mean()), "std": float(div_delta.std())},
        },
        "alignment": {
            "n_common_ca": len(common),
        }
    }
    with (out_dir / "SUMMARY.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[OK] Saved to: {out_dir}")
    print("[WT] ", info_wt)
    print("[MUT]", info_mut)
    print("Divergence WT  :", summary["div_stats"]["wt"])
    print("Divergence MUT :", summary["div_stats"]["mut"])
    print("Divergence Δ   :", summary["div_stats"]["delta"])


if __name__ == "__main__":
    main()
