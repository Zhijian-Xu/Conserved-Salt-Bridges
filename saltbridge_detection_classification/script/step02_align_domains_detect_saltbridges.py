#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import math
import argparse
import logging
from math import sqrt
from pymol import cmd, stored

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DEG_TO_RAD = math.pi / 180.0
RAD_TO_DEG = 180.0 / math.pi


def get_coord(selection_tmp):
    stored.xyz = []
    cmd.iterate_state(1, selection_tmp, "stored.xyz.append([x,y,z])")
    return stored.xyz[0]


def analyze_protein_saltbridge(protein_obj: str, out_fh, com_mod):
    """
    对单个 protein_obj 检测盐桥：
    Asp/Glu (CG/CD “center”) 与 Lys/Arg (NZ/CZ) 在 5.0 Å 内
    并输出 bridge_com（两点中点）
    """
    sele_Asp_CG_Glu_CD = f"((r. Asp and name CG) or (r. Glu and name CD)) and polymer and {protein_obj}"
    cmd.select("Asp_CG_Glu_CD", sele_Asp_CG_Glu_CD)

    if cmd.count_atoms("Asp_CG_Glu_CD") <= 0:
        out_fh.write(f"No atoms selected with Asp_CG_Glu_CD for protein {protein_obj}\n")
        return

    com_mod.com("Asp_CG_Glu_CD", mass=1, object="center_of_ring")
    CG_CD = cmd.get_model("Asp_CG_Glu_CD")

    for at in CG_CD.atom:
        if at.name == "CG":
            sele_ring = f"{protein_obj} and ((chain {at.chain} and resi {at.resi}) and name CG+OD1+OD2)"
        elif at.name == "CD":
            sele_ring = f"{protein_obj} and ((chain {at.chain} and resi {at.resi}) and name CD+OE1+OE2)"
        else:
            continue

        com_mod.com(sele_ring, mass=1, object="center_of_ring")
        center_xyz = get_coord("center_of_ring")

        sele_cations = f"({protein_obj} within 5.0 of center_of_ring) and ((r. Lys and name NZ) or (r. Arg and name CZ))"
        cmd.select("cations", sele_cations)
        if cmd.count_atoms("cations") <= 0:
            continue

        NZ_CZ = cmd.get_model("cations")
        for at_cation in NZ_CZ.atom:
            if at_cation.alt == "":
                sele_cation = f"{protein_obj} and (cations and chain {at_cation.chain} and resi {at_cation.resi} and name {at_cation.name})"
            else:
                sele_cation = f"{protein_obj} and (cations and chain {at_cation.chain} and resi {at_cation.resi} and name {at_cation.name} and alt {at_cation.alt})"

            cmd.select("cation", sele_cation)
            cation_xyz = get_coord("cation")

            dist = sqrt((cation_xyz[0] - center_xyz[0])**2 +
                        (cation_xyz[1] - center_xyz[1])**2 +
                        (cation_xyz[2] - center_xyz[2])**2)

            bridge_com = [
                (cation_xyz[0] + center_xyz[0]) / 2,
                (cation_xyz[1] + center_xyz[1]) / 2,
                (cation_xyz[2] + center_xyz[2]) / 2,
            ]

            out_fh.write(
                f"{protein_obj} pi_chain:{at.chain} pi_resn:{at.resn} pi_resi:{at.resi} "
                f"cation_chain:{at_cation.chain} cation_resn:{at_cation.resn} cation_resi:{at_cation.resi} "
                f"distance:{dist:.2f} bridge_com: [{bridge_com[0]:.2f},{bridge_com[1]:.2f},{bridge_com[2]:.2f}]\n"
            )


def parse_domain_names(list_txt: str):
    names = []
    with open(list_txt, "r", encoding="utf-8") as f:
        for line in f:
            fields = line.split()
            if not fields:
                continue

            if len(fields) == 3:
                pdbid, chain, r = fields
                names.append(f"{pdbid}-{chain}{r.replace('-','_')}")
            elif len(fields) == 4:
                pdbid, chain, r1, r2 = fields
                names.append(f"{pdbid}-{chain}{r1.replace('-','_')}-{r2.replace('-','_')}")
            elif len(fields) == 5:
                pdbid, c1, r1, c2, r2 = fields
                names.append(f"{pdbid}-{c1}{r1.replace('-','_')}-{c2}{r2.replace('-','_')}")
            else:
                logger.warning(f"[WARN] invalid line: {line.strip()}")
    return names


def main():
    parser = argparse.ArgumentParser(
        description="Step02: Align domains within each family and detect salt bridges"
    )
    parser.add_argument("--list-dir", required=True,
                        help="Directory containing <family_id>.txt (domain list)")
    parser.add_argument("--domain-dir", required=True,
                        help="Directory containing domains: <domain-dir>/<family_id>/<domain_name>.cif")
    parser.add_argument("--output", required=True,
                        help="Output directory to write <family_id>.txt salt bridge results")
    parser.add_argument("--utils-dir", required=True,
                        help="Directory containing com.py and com2.py (required)")
    args = parser.parse_args()

    list_dir = os.path.abspath(args.list_dir)
    domain_dir = os.path.abspath(args.domain_dir)
    out_dir = os.path.abspath(args.output)
    utils_dir = os.path.abspath(args.utils_dir)

    os.makedirs(out_dir, exist_ok=True)

    # 导入 com/com2
    sys.path.insert(0, utils_dir)
    import com  # noqa

    cmd.set("ignore_case", 1)

    for fn in sorted(os.listdir(list_dir)):
        if not fn.endswith(".txt"):
            continue

        family_id = os.path.splitext(fn)[0]
        list_txt = os.path.join(list_dir, fn)
        fam_domain_dir = os.path.join(domain_dir, family_id)
        if not os.path.isdir(fam_domain_dir):
            logger.warning(f"[WARN] domain folder not found: {fam_domain_dir}")
            continue

        domain_names = parse_domain_names(list_txt)
        if not domain_names:
            continue

        out_path = os.path.join(out_dir, f"{family_id}.txt")
        logger.info(f"=== Family {family_id} -> {out_path}")

        cmd.reinitialize()
        cmd.feedback("disable", "all", "everything")

        with open(out_path, "w", encoding="utf-8") as out_fh:
            ref_obj = None

            for i, name in enumerate(domain_names):
                cif_path = os.path.join(fam_domain_dir, f"{name}.cif")
                if not os.path.exists(cif_path):
                    logger.warning(f"[WARN] missing domain file: {cif_path}")
                    continue

                if i == 0:
                    ref_obj = name
                    cmd.load(cif_path, ref_obj)
                    analyze_protein_saltbridge(ref_obj, out_fh, com)
                else:
                    cmd.load(cif_path, name)
                    rmsd = cmd.super(name, ref_obj)[0]
                    out_fh.write(f"{name} aligned_to {ref_obj} RMSD:{rmsd:.2f}\n")
                    analyze_protein_saltbridge(name, out_fh, com)
                    cmd.delete(name)

        cmd.delete("all")


if __name__ == "__main__":
    main()
