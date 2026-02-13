# MdPath Information-Flow Analysis (WT vs MUT)

This submodule runs **MdPath** on WT and mutant MD trajectories, then quantifies **information-field reorganization** using:
1) divergence-field comparison (Mut − WT)  
2) residue-level attribution around Cα atoms  
3) replicate summarization (mean/SD/CI across 3 runs)

> Commands below are written **without hard-coded absolute paths**.  
> Use environment variables or relative paths so the workflow is portable.

---

## 0) Prerequisites

### Software
- `mdpath` executable available in `$PATH` (or provide full path to binary)
- Python ≥ 3.9 with:
  - numpy, pandas
  - scipy (recommended; required if you use smoothing in divergence-field step)

Note: `compare_divergence_field.py` applies Gaussian smoothing via `scipy.ndimage.gaussian_filter` when `--smooth_sigma > 0` fileciteturn5file0L17-L23.  
If you do not want scipy, set `--smooth_sigma 0`.

### Required MdPath outputs per run
For each run directory (WT or MUT), the scripts assume:
- `nmi_df.csv`
- `output.txt`

`compare_divergence_field.py` explicitly looks for these two files in both `--wt_dir` and `--mut_dir` fileciteturn5file0L313-L321.

---

## 1) Recommended directory layout

```text
mdpath_pipeline/
├── wt/
│   ├── rep1/              # mdpath output folder (contains nmi_df.csv + output.txt)
│   ├── rep2/
│   └── rep3/
├── mut/
│   ├── rep1/
│   ├── rep2/
│   └── rep3/
├── pdb/
│   ├── wt_first_frame.pdb
│   └── mut_first_frame.pdb
├── traj/
│   ├── wt_align.xtc
│   └── mut_align.xtc
├── scripts/
│   ├── compare_divergence_field.py
│   ├── residue_level_from_fields.py
│   └── summarize_infofield_reps.py
└── results/
    ├── rep1/
    ├── rep2/
    ├── rep3/
    └── summary/
```

---

## 2) Step01 — Run MdPath (mutant trajectories)

You provided two mutant systems. Example `nohup` submission commands:

### System A: 6FYV-E79A

```bash
nohup mdpath   -top 6FYV-E79A.pdb   -traj ../MD1-200ns/6FYV-E79A1_align.xtc   -cpu 20   -fardist 40   -lig 8 79 92 94 41 56 138 175   -numpath 500   -graphdist 8   -closedist 4   > mdpath.out 2>&1 &
```

### System B: 4IAN-R79A (trajectory file name uses R80A3_align.xtc in your example)

```bash
nohup mdpath   -top 4IAN-R79A.pdb   -traj ../MD3-200ns/4IAN-R80A3_align.xtc   -cpu 20   -fardist 40   -lig 8 79 41 56 158 157 93 31 33   -numpath 500   -graphdist 8   -closedist 4   > mdpath.out 2>&1 &
```

**Outputs (per run folder)**: `nmi_df.csv` + `output.txt` (used by downstream scripts).

---

## 3) Step02 — Divergence-field comparison (Mut − WT)

This step:
1) extracts Cα coordinates from WT and MUT PDBs
2) aligns MUT→WT by Kabsch alignment on common residues
3) builds 3D vector fields from MdPath edges weighted by NMI
4) computes divergence fields for WT and MUT
5) saves Δdivergence = div(MUT) − div(WT)

The script stores all arrays and a `SUMMARY.json` under the output directory fileciteturn5file0L404-L437.

### Command template (portable)

Set your paths as variables:

```bash
WT_DIR="wt/rep1"
MUT_DIR="mut/rep1"
WT_PDB="pdb/wt_first_frame.pdb"
MUT_PDB="pdb/mut_first_frame.pdb"
OUT_DIR="results/rep1/DIVERGENCE_COMPARE"
```

Run:

```bash
python3 scripts/compare_divergence_field.py   --wt_dir  "${WT_DIR}"   --mut_dir "${MUT_DIR}"   --wt_pdb  "${WT_PDB}"   --mut_pdb "${MUT_PDB}"   --pdb_index_mode sequential   --shift_mut 0   --grid_h 1.0   --smooth_sigma 1.0   --vector_mode unit   --out_dir "${OUT_DIR}"
```

### Key options you used
- `--pdb_index_mode sequential` (recommended for mdpath “Res i”) fileciteturn5file0L292-L299  
- `--shift_mut 0` (if indices already match) fileciteturn5file0L297-L300  
- `--grid_h 1.0` voxel size in Å fileciteturn5file0L299-L301  
- `--smooth_sigma 1.0` smoothing sigma (set 0 to disable) fileciteturn5file0L302-L305  
- `--vector_mode unit` uses NMI * unit(direction) fileciteturn5file0L305-L308  

### Outputs (inside `${OUT_DIR}`)
- `wt_Fx.npy, wt_Fy.npy, wt_Fz.npy`
- `mut_Fx.npy, mut_Fy.npy, mut_Fz.npy`
- `wt_div.npy`, `mut_div.npy`
- `delta_div_mut_minus_wt.npy`
- `SUMMARY.json`

Saved exactly by the script fileciteturn5file0L404-L437.

---

## 4) Step03 — Residue-level attribution from fields

This step samples (via trilinear interpolation) the divergence change and |ΔF| near each residue Cα, and reports:
- `delta_div_at_CA`
- `mean_abs_delta_div_R{R}`
- `dFmag_at_CA`
- `mean_dFmag_R{R}`

It writes:
- `residue_level_changes_R{R}.csv`
- `top{N}_residues_by_mean_abs_delta_div_R{R}.csv` fileciteturn5file1L133-L146

### Command template

```bash
DIV_DIR="results/rep1/DIVERGENCE_COMPARE"
WT_PDB="pdb/wt_first_frame.pdb"

python3 scripts/residue_level_from_fields.py   --div_dir "${DIV_DIR}"   --wt_pdb "${WT_PDB}"   --pdb_index_mode sequential   --R 3.0   --top 50
```

---

## 5) Step04 — Summarize 3 replicates (mean across reps)

You run three independent MUT replicates (rep1/rep2/rep3). Each replicate should produce:

```
results/repX/DIVERGENCE_COMPARE/residue_level_changes_R3.0.csv
```

Then summarize them:

```bash
python3 scripts/summarize_infofield_reps.py   --csv     results/rep1/DIVERGENCE_COMPARE/residue_level_changes_R3.0.csv     results/rep2/DIVERGENCE_COMPARE/residue_level_changes_R3.0.csv     results/rep3/DIVERGENCE_COMPARE/residue_level_changes_R3.0.csv   --outprefix 6FYV   --rank_by mean_abs_delta_div_R3.0_mean,mean_dFmag_R3.0_mean   --topn 10   --motif_json motifs_6FYV.json
```

This script outputs:
- `<outprefix>.residue_summary.csv` (per-residue mean/sd/CI across reps) fileciteturn5file2L1-L11  
- `<outprefix>.global_per_rep.csv` and `<outprefix>.global_summary.csv` fileciteturn5file2L1-L11  
- optional motif summaries when `--motif_json` is provided fileciteturn5file2L1-L11  

It also prints the **Top-N residues** ranked by one or more summary columns to terminal fileciteturn5file2L14-L21.

### Motif JSON format (example)

```json
{
  "agg": "max",
  "motifs": {
    "β3–VAIK": [41],
    "Hinge": [94, 95],
    "DFG": [184, 185, 186]
  }
}
```

- `agg="max"` means take the maximum metric among residues in the motif.
- `agg="mean"` means take the average.

(See docstring in script header) fileciteturn5file2L68-L95

---

## 6) Practical notes (to avoid common pitfalls)

### A) Avoid absolute paths
Your original commands had absolute paths (e.g., `/home/databank/...`).  
For portability, replace them with variables like `${WT_DIR}` and relative paths as shown above.

### B) Index consistency matters
If residue indices between WT and MUT are off by a constant, use `--shift_mut` to fix:
> WT_index = MUT_index + shift_mut fileciteturn5file0L297-L300

If indices are from PDB `resSeq`, switch:
- `--pdb_index_mode resseq`
Otherwise keep:
- `--pdb_index_mode sequential` fileciteturn5file0L292-L299

### C) Alignment requires enough common residues
The script requires ≥20 common Cα residues for alignment, otherwise it errors (check shift/index mode) fileciteturn5file0L337-L343.

### D) Runtime
- `compare_divergence_field.py` cost scales with grid size (`--grid_h`, `--padding`) and number of edges deposited.
- `residue_level_from_fields.py` cost scales with number of residues and radius `--R`.

---

## 7) What to report in the paper

Typical reporting items:
- Δdivergence field (Mut − WT) visualization (hotspots around salt-bridge neighborhood)
- Top residues by `mean_abs_delta_div_R3.0_mean` and `mean_dFmag_R3.0_mean`
- Motif-level summaries (β3–VAIK / αC / HRD / DFG / hinge, etc.)
- Replicate-averaged statistics and confidence intervals

---

## License / data note

Large MdPath outputs and MD trajectories are not committed to GitHub.  
Only scripts and small metadata are tracked.
