#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import pandas as pd
from collections import defaultdict


def parse_variant(variant: str):
    # e.g. "M1A" -> ("M", 1, "A")
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", str(variant).strip())
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), m.group(3)


def main():
    ap = argparse.ArgumentParser(description="Extract AlphaMissense X->A scores for salt-bridge UniProt positions.")
    ap.add_argument("--mapping", required=True, help="mapping CSV from Step01")
    ap.add_argument("--am-tsv", required=True, help="AlphaMissense_aa_substitutions.tsv")
    ap.add_argument("--output", required=True, help="Output CSV (matched)")
    ap.add_argument("--not-found", default=None, help="Output CSV for not found (optional)")
    ap.add_argument("--non-saltbridge", default=None, help="Output CSV for non-saltbridge residues (optional)")
    ap.add_argument("--chunksize", type=int, default=2_000_000, help="TSV chunk size for streaming")
    args = ap.parse_args()

    map_df = pd.read_csv(args.mapping)
    required_cols = {"uniprot_id", "uniprot_pos"}
    if not required_cols.issubset(set(map_df.columns)):
        raise SystemExit(f"ERROR: mapping must contain columns: {sorted(required_cols)}")

    # positions to query
    required = set((str(r.uniprot_id), int(r.uniprot_pos)) for r in map_df.itertuples(index=False))
    print(f"[Step02] unique UniProt positions to query: {len(required)}")

    # Load AM TSV in chunks
    # Expected AM TSV columns (common):
    #   uniprot_id, protein_variant, am_pathogenicity, am_class
    keep_cols = ["uniprot_id", "protein_variant", "am_pathogenicity", "am_class"]
    found = defaultdict(list)

    print("[Step02] scanning AlphaMissense TSV (only X->A)...")
    reader = pd.read_csv(args.am_tsv, sep="\t", usecols=keep_cols, chunksize=args.chunksize, dtype={"uniprot_id": str})
    total = 0
    matched_lines = 0

    for chunk in reader:
        total += len(chunk)
        # parse variant -> new_aa == "A"
        # cheap filter first:
        chunk = chunk[chunk["protein_variant"].astype(str).str.endswith("A")]
        if chunk.empty:
            continue

        # parse positions
        for row in chunk.itertuples(index=False):
            unip = str(row.uniprot_id)
            oaa, pos, naa = parse_variant(row.protein_variant)
            if pos is None or naa != "A":
                continue
            key = (unip, pos)
            if key in required:
                found[key].append((oaa, naa, row.am_pathogenicity, row.am_class))
                matched_lines += 1

    print(f"[Step02] AM TSV scanned lines={total}, matched X->A records={matched_lines}, matched positions={len(found)}")

    saltbridge_res = {"K", "D", "E", "R"}  # your definition
    matched_rows = []
    non_sb_rows = []
    not_found_rows = []

    # dedupe by (unip, pos, oaa, naa)
    seen = set()
    for r in map_df.itertuples(index=False):
        unip = str(r.uniprot_id)
        pos = int(r.uniprot_pos)
        key = (unip, pos)

        if key not in found:
            not_found_rows.append(
                {
                    "uniprot_id": unip,
                    "uniprot_pos": pos,
                    "reason": "not_found_in_alphamissense",
                }
            )
            continue

        for (oaa, naa, score, am_class) in found[key]:
            dkey = (unip, pos, oaa, naa)
            if dkey in seen:
                continue
            seen.add(dkey)

            base = {
                "uniprot_id": unip,
                "uniprot_pos": pos,
                "wt_aa": oaa,
                "mut_aa": naa,
                "am_pathogenicity": score,
                "am_class": am_class,
            }
            # carry mapping meta if exists
            for col in map_df.columns:
                if col in base:
                    continue
                base[col] = getattr(r, col)
            if oaa in saltbridge_res:
                matched_rows.append(base)
            else:
                base2 = base.copy()
                base2["reason"] = f"wt_aa_not_in_saltbridge_set({sorted(saltbridge_res)})"
                non_sb_rows.append(base2)

    out = pd.DataFrame(matched_rows)
    out_path = args.output
    pd.DataFrame(matched_rows).to_csv(out_path, index=False, encoding="utf-8")
    print(f"[Step02] wrote matched: {out_path} rows={len(out)}")

    nf_path = args.not_found or (out_path.replace(".csv", "_not_found.csv"))
    pd.DataFrame(not_found_rows).drop_duplicates().to_csv(nf_path, index=False, encoding="utf-8")
    print(f"[Step02] wrote not_found: {nf_path} rows={len(not_found_rows)}")

    nsb_path = args.non_saltbridge or (out_path.replace(".csv", "_non_saltbridge.csv"))
    pd.DataFrame(non_sb_rows).to_csv(nsb_path, index=False, encoding="utf-8")
    print(f"[Step02] wrote non_saltbridge: {nsb_path} rows={len(non_sb_rows)}")


if __name__ == "__main__":
    main()
