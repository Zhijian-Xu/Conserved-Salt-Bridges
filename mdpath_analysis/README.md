# MdPath Information-Flow Analysis (WT vs MUT)

This module performs **information-flow analysis** using MdPath on wild-type (WT) and mutant (MUT) molecular dynamics trajectories, followed by divergence-field comparison and residue-level attribution.

---

# Environment Setup (MdPath)

MdPath must be installed before running this pipeline.

Official repository:

https://github.com/wolberlab/mdpath

Follow the installation instructions provided in the MdPath repository to compile and configure the executable.

After installation, ensure the binary is accessible:

```bash
mdpath --help
```

If not, add it to your PATH:

```bash
export PATH=/path/to/mdpath:$PATH
```

---

# Data Policy

Due to storage constraints, **MD trajectories are NOT distributed in this repository.**

Only the following are tracked:

- original protein structure files (PDB)
- analysis scripts
- lightweight metadata

Users must generate trajectories independently before running MdPath.

---

# Pipeline Overview

1. Run MdPath on WT and mutant trajectories  
2. Compare divergence fields (Mut − WT)  
3. Compute residue-level information changes  
4. Average results across replicate simulations  

---

# Recommended Directory Layout

```
mdpath_pipeline/
├── wt/
│   ├── rep1/
│   ├── rep2/
│   └── rep3/
├── mut/
│   ├── rep1/
│   ├── rep2/
│   └── rep3/
├── pdb/
│   ├── wt_first_frame.pdb
│   └── mut_first_frame.pdb
├── scripts/
│   ├── compare_divergence_field.py
│   ├── residue_level_from_fields.py
│   └── summarize_infofield_reps.py
└── results/
```

---

# Step01 — Run MdPath

Example submissions.

## System A — 6FYV-E79A

```bash
nohup mdpath   -top 6FYV-E79A.pdb   -traj ../MD1-200ns/6FYV-E79A1_align.xtc   -cpu 20   -fardist 40   -lig 8 79 92 94 41 56 138 175   -numpath 500   -graphdist 8   -closedist 4   > mdpath.out 2>&1 &
```

## System B — 4IAN-R79A

```bash
nohup mdpath   -top 4IAN-R79A.pdb   -traj ../MD3-200ns/4IAN-R80A3_align.xtc   -cpu 20   -fardist 40   -lig 8 79 41 56 158 157 93 31 33   -numpath 500   -graphdist 8   -closedist 4   > mdpath.out 2>&1 &
```

Each run produces:

- `nmi_df.csv`
- `output.txt`

These files are required for downstream analysis.

---

# Step02 — Divergence Field Comparison

```bash
python3 scripts/compare_divergence_field.py   --wt_dir wt/rep1   --mut_dir mut/rep1   --wt_pdb pdb/wt_first_frame.pdb   --mut_pdb pdb/mut_first_frame.pdb   --pdb_index_mode sequential   --shift_mut 0   --grid_h 1.0   --smooth_sigma 1.0   --vector_mode unit
```

Output directory:

```
DIVERGENCE_COMPARE/
```

Key outputs:

- divergence fields (WT and MUT)
- delta divergence (Mut − WT)
- SUMMARY.json

---

# Step03 — Residue-Level Information Change

```bash
python3 scripts/residue_level_from_fields.py   --div_dir DIVERGENCE_COMPARE   --wt_pdb pdb/wt_first_frame.pdb   --pdb_index_mode sequential   --R 3.0   --top 50
```

Output:

```
residue_level_changes_R3.0.csv
```

---

# Step04 — Average Across Replicates

Run after generating residue tables for three independent simulations.

```bash
python3 scripts/summarize_infofield_reps.py   --csv     results/rep1/DIVERGENCE_COMPARE/residue_level_changes_R3.0.csv     results/rep2/DIVERGENCE_COMPARE/residue_level_changes_R3.0.csv     results/rep3/DIVERGENCE_COMPARE/residue_level_changes_R3.0.csv   --outprefix 6FYV   --rank_by mean_abs_delta_div_R3.0_mean,mean_dFmag_R3.0_mean   --topn 10   --motif_json motifs_6FYV.json
```

Outputs include:

- residue_summary.csv  
- global_summary.csv  

These represent replicate-averaged information-field changes.

---

# Reproducibility

All analyses can be reproduced by:

1. generating MD trajectories  
2. running MdPath  
3. executing the scripts in this directory  

Large trajectory files are intentionally excluded from version control.
