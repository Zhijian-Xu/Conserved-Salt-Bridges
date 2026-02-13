#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import pandas as pd


def load_fasta_map(fasta_path: Path):
    """
    Load a FASTA and return dict {uniprot_id: sequence}
    Accepts header like:
      >sp|P12345|...
      >P12345 ...
    """
    seqs = {}
    cur_id = None
    cur = []
    with fasta_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush
                if cur_id and cur:
                    seqs[cur_id] = "".join(cur)
                cur = []
                header = line[1:]
                # parse id
                if "|" in header:
                    parts = header.split("|")
                    # sp|P12345|...
                    cur_id = parts[1] if len(parts) >= 2 else parts[0].split()[0]
                else:
                    cur_id = header.split()[0]
            else:
                cur.append(line)
        if cur_id and cur:
            seqs[cur_id] = "".join(cur)
    return seqs


def main():
    ap = argparse.ArgumentParser(description="Prepare ESM DMS input CSV from mapping + UniProt FASTA.")
    ap.add_argument("--mapping", required=True, help="mapping.csv from Step01")
    ap.add_argument("--fasta", required=True, help="UniProt FASTA (contains sequences for needed UniProt IDs)")
    ap.add_argument("--output", required=True, help="Output ESM input CSV")
    ap.add_argument("--mutation-aa", default="A", help="Mutate to this amino acid (default A)")
    args = ap.parse_args()

    df = pd.read_csv(args.mapping)
    if not {"uniprot_id", "uniprot_pos"}.issubset(df.columns):
        raise SystemExit("ERROR: mapping must contain columns: uniprot_id, uniprot_pos")

    seq_map = load_fasta_map(Path(args.fasta))
    print(f"[Step04] loaded FASTA sequences: {len(seq_map)}")

    rows = []
    missing_seq = 0
    out_of_range = 0

    for r in df.itertuples(index=False):
        uid = str(r.uniprot_id)
        pos = int(r.uniprot_pos)
        seq = seq_map.get(uid)
        if not seq:
            missing_seq += 1
            continue
        if pos < 1 or pos > len(seq):
            out_of_range += 1
            continue
        wt = seq[pos - 1]
        mutant = f"{wt}{pos}{args.mutation_aa}"
        rows.append(
            {
                "uniprot_id": uid,
                "uniprot_pos": pos,
                "sequence": seq,
                "mutant": mutant,
            }
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).drop_duplicates().to_csv(out, index=False, encoding="utf-8")
    print(f"[Step04] wrote: {out} rows={len(rows)}")
    print(f"[Step04] missing_seq={missing_seq}, out_of_range={out_of_range}")


if __name__ == "__main__":
    main()
