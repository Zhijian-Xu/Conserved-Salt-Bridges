#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path
import pandas as pd


def parse_range(r: str):
    """Parse 'start-end' -> (start, end) ints."""
    if not r or "-" not in r:
        return None, None
    try:
        a, b = r.split("-", 1)
        return int(a), int(b)
    except Exception:
        return None, None


def convert_pdb_to_uniprot(pdb_pos: int, pdb_start: int, pdb_end: int, uniprot_start: int, uniprot_end: int):
    """
    Linear mapping by offset:
      pdb_start <-> uniprot_start
      pdb_end   <-> uniprot_end
    """
    if pdb_pos < pdb_start or pdb_pos > pdb_end:
        return None, "position_out_of_pdb_range"
    rel = pdb_pos - pdb_start
    uniprot_pos = uniprot_start + rel
    if uniprot_pos < uniprot_start or uniprot_pos > uniprot_end:
        return None, "mapped_out_of_uniprot_range"
    return uniprot_pos, "ok"


def parse_scop_file(scop_file_path: Path):
    """
    Build mapping:
      key: (PDBID, 'A:665-730')
      val: (UniProtID, '1-65', family_id)
    Accepts lines like (space-delimited):
      <scop_id> <pdbid> <chain:range> <uniprot> <uniprot_range> <family=xxxx>
    """
    mapping = {}
    with scop_file_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            pdb_id = parts[1]
            chain_info = parts[2]  # e.g. A:665-730
            uniprot_id = parts[3]
            uniprot_range = parts[4]
            fam_raw = parts[5]
            family_id = fam_raw.split("=", 1)[1] if "=" in fam_raw else fam_raw

            if ":" in chain_info and "-" in chain_info:
                mapping[(pdb_id, chain_info)] = (uniprot_id, uniprot_range, family_id)
    return mapping


def parse_saltbridge_file(file_path: Path):
    """
    Your original saltbridge txt format: each line has 10 fields
    and first field like: 5BP1-A665_730
    and contains fields like 'pi_resi:###' and 'cation_resi:###'
    """
    saltbridges = []
    seen = set()
    txt = file_path.read_text(encoding="utf-8", errors="ignore").strip().splitlines()
    for line in txt:
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) != 10:
            continue
        if line in seen:
            continue
        seen.add(line)

        first = fields[0]
        if "-" not in first:
            continue
        pdb_id, chain_seq = first.split("-", 1)

        m = re.match(r"([A-Z])(\d+)_(\d+)", chain_seq)
        if not m:
            continue
        chain = m.group(1)
        start_pos = int(m.group(2))
        end_pos = int(m.group(3))
        chain_range = f"{chain}:{start_pos}-{end_pos}"

        pi_resi = None
        cation_resi = None
        for field in fields:
            if field.startswith("pi_resi:"):
                pi_resi = int(field.split(":", 1)[1])
            elif field.startswith("cation_resi:"):
                cation_resi = int(field.split(":", 1)[1])

        if pi_resi is None or cation_resi is None:
            continue

        saltbridges.append(
            {
                "pdb_id": pdb_id,
                "chain_range": chain_range,
                "pi_resi": pi_resi,
                "cation_resi": cation_resi,
                "raw_line": line,
            }
        )
    return saltbridges


def infer_family_id_from_filename(p: Path) -> str:
    # default: use stem; you can customize here
    return p.stem


def main():
    ap = argparse.ArgumentParser(description="Map salt-bridge PDB residues to UniProt positions using SCOP mapping.")
    ap.add_argument("--saltbridge-dir", type=str, help="Directory containing *.txt saltbridge files")
    ap.add_argument("--saltbridge-file", type=str, help="Single saltbridge txt file")
    ap.add_argument("--scop-file", type=str, required=True, help="SCOP mapping file (e.g., scop-cla-latest.txt)")
    ap.add_argument("--output", type=str, required=True, help="Output mapping CSV path")
    ap.add_argument("--error-output", type=str, default=None, help="Output errors CSV path (optional)")
    args = ap.parse_args()

    if not args.saltbridge_dir and not args.saltbridge_file:
        raise SystemExit("ERROR: provide --saltbridge-dir or --saltbridge-file")

    scop_file = Path(args.scop_file)
    if not scop_file.exists():
        raise SystemExit(f"ERROR: SCOP file not found: {scop_file}")

    print("[Step01] Parsing SCOP mapping...")
    scop_mapping = parse_scop_file(scop_file)
    print(f"[Step01] SCOP mappings loaded: {len(scop_mapping)}")

    inputs = []
    if args.saltbridge_file:
        inputs = [Path(args.saltbridge_file)]
    else:
        d = Path(args.saltbridge_dir)
        if not d.exists():
            raise SystemExit(f"ERROR: saltbridge dir not found: {d}")
        inputs = sorted(d.glob("*.txt"))

    results = []
    errors = []

    for f in inputs:
        family_id = infer_family_id_from_filename(f)
        sbs = parse_saltbridge_file(f)
        print(f"[Step01] {f.name}: unique saltbridges={len(sbs)}")

        for sb in sbs:
            pdb_id = sb["pdb_id"]
            chain_range = sb["chain_range"]
            pi_resi = sb["pi_resi"]
            cation_resi = sb["cation_resi"]

            key = (pdb_id, chain_range)
            if key not in scop_mapping:
                errors.append(
                    {
                        "family_id": family_id,
                        "pdb_id": pdb_id,
                        "chain_range": chain_range,
                        "error_type": "scop_mapping_not_found",
                        "detail": f"{pdb_id} {chain_range}",
                    }
                )
                continue

            uniprot_id, uniprot_range, scop_family_id = scop_mapping[key]

            chain_part, pdb_range_str = chain_range.split(":", 1)
            pdb_start, pdb_end = parse_range(pdb_range_str)
            uniprot_start, uniprot_end = parse_range(uniprot_range)

            if pdb_start is None or uniprot_start is None:
                errors.append(
                    {
                        "family_id": family_id,
                        "pdb_id": pdb_id,
                        "chain_range": chain_range,
                        "error_type": "range_parse_failed",
                        "detail": f"pdb={pdb_range_str}, uniprot={uniprot_range}",
                    }
                )
                continue

            for site_type, pdb_pos in [("pi_resi", pi_resi), ("cation_resi", cation_resi)]:
                uni_pos, status = convert_pdb_to_uniprot(pdb_pos, pdb_start, pdb_end, uniprot_start, uniprot_end)
                if uni_pos is None:
                    errors.append(
                        {
                            "family_id": family_id,
                            "pdb_id": pdb_id,
                            "chain_range": chain_range,
                            "uniprot_id": uniprot_id,
                            "uniprot_range": uniprot_range,
                            "pdb_pos": pdb_pos,
                            "site_type": site_type,
                            "error_type": "convert_failed",
                            "detail": status,
                        }
                    )
                    continue

                results.append(
                    {
                        "family_id": family_id,
                        "pdb_id": pdb_id,
                        "chain_range": chain_range,
                        "uniprot_id": uniprot_id,
                        "uniprot_range": uniprot_range,
                        "pdb_pos": pdb_pos,
                        "uniprot_pos": uni_pos,
                        "site_type": site_type,
                        "scop_family_id": scop_family_id,
                    }
                )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(out, index=False, encoding="utf-8")
    print(f"[Step01] Wrote mapping: {out} (rows={len(results)})")

    err_out = Path(args.error_output) if args.error_output else out.with_name(out.stem + "_errors.csv")
    pd.DataFrame(errors).to_csv(err_out, index=False, encoding="utf-8")
    print(f"[Step01] Wrote errors:  {err_out} (rows={len(errors)})")


if __name__ == "__main__":
    main()
