#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import pymol.cmd as cmd


def split_one_family(family_dir: str, list_txt: str, out_family_dir: str):
    """
    family_dir:    某个家族的原始结构目录（里面有 *.cif.gz）
    list_txt:      该家族的 domain 列表（每行：PDBID CHAIN RANGE...）
    out_family_dir:输出目录（该家族的 domain cif 输出到这里）
    """
    os.makedirs(out_family_dir, exist_ok=True)

    if not os.path.exists(list_txt):
        print(f"[WARN] list file not found: {list_txt}")
        return

    with open(list_txt, "r", encoding="utf-8") as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue

            # 支持 3/4/5 列
            if len(fields) == 3:
                pdbid, chain, res_range = fields
                selection_name = f"{pdbid}-{chain}{res_range.replace('-','_')}"
                sele = f"model {pdbid} and chain {chain} and resi {res_range}"

            elif len(fields) == 4:
                pdbid, chain, r1, r2 = fields
                selection_name = f"{pdbid}-{chain}{r1.replace('-','_')}-{r2.replace('-','_')}"
                sele = f"model {pdbid} and chain {chain} and (resi {r1} or resi {r2})"

            elif len(fields) == 5:
                pdbid, chain1, r1, chain2, r2 = fields
                selection_name = f"{pdbid}-{chain1}{r1.replace('-','_')}-{chain2}{r2.replace('-','_')}"
                sele = f"model {pdbid} and ((chain {chain1} and resi {r1}) or (chain {chain2} and resi {r2}))"

            else:
                print(f"[WARN] invalid line in {list_txt}: {line.strip()}")
                continue

            # 找到对应的 cif.gz（忽略大小写）
            pattern = re.compile(fr"^{re.escape(pdbid)}\.cif\.gz$", re.IGNORECASE)
            cif_gz_path = None
            for fn in os.listdir(family_dir):
                if pattern.match(fn):
                    cif_gz_path = os.path.join(family_dir, fn)
                    break

            if cif_gz_path is None:
                print(f"[WARN] {pdbid}.cif.gz not found under {family_dir}")
                continue

            try:
                cmd.delete("all")
                cmd.load(cif_gz_path, pdbid)
                cmd.select(selection_name, sele)

                out_path = os.path.join(out_family_dir, f"{selection_name}.cif")
                cmd.save(out_path, selection_name)
                print(f"[OK] saved: {out_path}")
            except Exception as e:
                print(f"[ERROR] failed on {cif_gz_path} line={line.strip()} err={e}")
            finally:
                cmd.delete("all")


def main():
    parser = argparse.ArgumentParser(
        description="Step01: Split SCOP domains into individual domain structure files (.cif)"
    )
    parser.add_argument("--input", required=True,
                        help="Input root dir of raw family structures. Each subdir is one family and contains *.cif.gz")
    parser.add_argument("--list-dir", required=True,
                        help="Directory containing <family_id>.txt (domain ranges per family)")
    parser.add_argument("--output", required=True,
                        help="Output root dir. Domains will be saved to <output>/<family_id>/*.cif")
    args = parser.parse_args()

    in_root = os.path.abspath(args.input)
    list_dir = os.path.abspath(args.list_dir)
    out_root = os.path.abspath(args.output)

    os.makedirs(out_root, exist_ok=True)

    # 遍历所有 family 子目录
    for fam in sorted(os.listdir(in_root)):
        fam_dir = os.path.join(in_root, fam)
        if not os.path.isdir(fam_dir):
            continue

        list_txt = os.path.join(list_dir, f"{fam}.txt")
        out_fam_dir = os.path.join(out_root, fam)

        print(f"\n=== Family {fam} ===")
        split_one_family(fam_dir, list_txt, out_fam_dir)


if __name__ == "__main__":
    main()
