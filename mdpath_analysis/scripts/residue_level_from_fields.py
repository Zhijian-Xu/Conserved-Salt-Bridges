#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def load_ca_coords(pdb_path: Path, index_mode="sequential"):
    coords = {}
    seen = set()
    seq = 0
    with pdb_path.open("r", errors="ignore") as f:
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
                seq += 1
                rid = seq
            else:
                rid = resseq
            coords[int(rid)] = np.array([x, y, z], float)
    if not coords:
        raise ValueError("No CA found in PDB.")
    return coords


def world_to_frac(origin, h, xyz):
    # voxel centers are at origin + (i+0.5)*h, so subtract 0.5
    return (xyz - origin) / h - 0.5


def trilinear(vol, origin, h, xyz):
    fi = world_to_frac(origin, h, xyz)
    i0 = np.floor(fi).astype(int)
    t = fi - i0
    nx, ny, nz = vol.shape
    i0 = np.clip(i0, [0, 0, 0], [nx - 2, ny - 2, nz - 2])
    x0, y0, z0 = i0
    tx, ty, tz = t

    c000 = vol[x0,   y0,   z0]
    c100 = vol[x0+1, y0,   z0]
    c010 = vol[x0,   y0+1, z0]
    c110 = vol[x0+1, y0+1, z0]
    c001 = vol[x0,   y0,   z0+1]
    c101 = vol[x0+1, y0,   z0+1]
    c011 = vol[x0,   y0+1, z0+1]
    c111 = vol[x0+1, y0+1, z0+1]

    c00 = c000*(1-tx) + c100*tx
    c10 = c010*(1-tx) + c110*tx
    c01 = c001*(1-tx) + c101*tx
    c11 = c011*(1-tx) + c111*tx

    c0 = c00*(1-ty) + c10*ty
    c1 = c01*(1-ty) + c11*ty

    return float(c0*(1-tz) + c1*tz)


def mean_abs_in_sphere(vol, origin, h, xyz, R):
    nx, ny, nz = vol.shape
    fc = world_to_frac(origin, h, xyz)
    ic = np.round(fc).astype(int)
    rad = int(np.ceil(R / h))

    i_min = max(0, ic[0] - rad); i_max = min(nx - 1, ic[0] + rad)
    j_min = max(0, ic[1] - rad); j_max = min(ny - 1, ic[1] + rad)
    k_min = max(0, ic[2] - rad); k_max = min(nz - 1, ic[2] + rad)

    acc = 0.0
    cnt = 0
    for i in range(i_min, i_max + 1):
        for j in range(j_min, j_max + 1):
            for k in range(k_min, k_max + 1):
                p = origin + h * (np.array([i, j, k], float) + 0.5)
                if np.linalg.norm(p - xyz) <= R:
                    acc += abs(float(vol[i, j, k]))
                    cnt += 1
    return acc / cnt if cnt else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--div_dir", required=True, help="DIVERGENCE_COMPARE directory")
    ap.add_argument("--wt_pdb", required=True, help="WT PDB (same index mode as mdpath)")
    ap.add_argument("--pdb_index_mode", default="sequential", choices=["sequential", "resseq"])
    ap.add_argument("--R", type=float, default=3.0, help="radius (Å) for neighborhood mean(|.|)")
    ap.add_argument("--top", type=int, default=50)
    args = ap.parse_args()

    d = Path(args.div_dir)
    summary = json.loads((d / "SUMMARY.json").read_text())
    origin = np.array(summary["grid"]["origin_mn"], float)
    h = float(summary["params"]["grid_h"])

    delta_div = np.load(d / "delta_div_mut_minus_wt.npy")
    wt_Fx = np.load(d / "wt_Fx.npy"); wt_Fy = np.load(d / "wt_Fy.npy"); wt_Fz = np.load(d / "wt_Fz.npy")
    mut_Fx = np.load(d / "mut_Fx.npy"); mut_Fy = np.load(d / "mut_Fy.npy"); mut_Fz = np.load(d / "mut_Fz.npy")

    # ΔF magnitude grid
    dFx = mut_Fx - wt_Fx
    dFy = mut_Fy - wt_Fy
    dFz = mut_Fz - wt_Fz
    dF_mag = np.sqrt(dFx*dFx + dFy*dFy + dFz*dFz)

    coords = load_ca_coords(Path(args.wt_pdb), index_mode=args.pdb_index_mode)

    rows = []
    for res, xyz in coords.items():
        delta_div_at_ca = trilinear(delta_div, origin, h, xyz)
        mean_abs_delta_div = mean_abs_in_sphere(delta_div, origin, h, xyz, R=args.R)

        dF_at_ca = trilinear(dF_mag, origin, h, xyz)
        mean_dF = mean_abs_in_sphere(dF_mag, origin, h, xyz, R=args.R)  # mean(|dF|) in sphere

        rows.append({
            "res": int(res),
            "x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2]),
            "delta_div_at_CA": float(delta_div_at_ca),
            f"mean_abs_delta_div_R{args.R:.1f}": float(mean_abs_delta_div),
            "dFmag_at_CA": float(dF_at_ca),
            f"mean_dFmag_R{args.R:.1f}": float(mean_dF),
        })

    df = pd.DataFrame(rows)

    out_csv = d / f"residue_level_changes_R{args.R:.1f}.csv"
    df.to_csv(out_csv, index=False)

    col = f"mean_abs_delta_div_R{args.R:.1f}"
    top_df = df.sort_values(col, ascending=False).head(args.top)
    out_top = d / f"top{args.top}_residues_by_{col}.csv"
    top_df.to_csv(out_top, index=False)

    print("[OK] wrote:")
    print(" ", out_csv)
    print(" ", out_top)
    print("\nTop residues:")
    print(top_df[["res", "delta_div_at_CA", col, "dFmag_at_CA", f"mean_dFmag_R{args.R:.1f}"]].to_string(index=False))


if __name__ == "__main__":
    main()
