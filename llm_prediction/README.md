# ML Prediction Pipeline  
## Deep-Learning Assessment of Spatially Conserved Salt Bridges

This module evaluates the mutational sensitivity of spatially conserved salt-bridge residues using state-of-the-art protein variant predictors:

- **AlphaMissense** — pathogenicity prediction  
- **ESM-1v** — zero-shot mutational effect estimation  

The pipeline maps salt-bridge residues from **PDB coordinates to UniProt positions**, retrieves mutation scores, and validates evolutionary constraints through large-scale language models.

---

# 🔬 Scientific Rationale

Spatially conserved salt bridges are hypothesized to represent structurally critical interactions under strong evolutionary constraint.

To test this hypothesis, we quantify mutational tolerance at these residues using:

- genome-scale pathogenicity predictions (AlphaMissense)
- protein language model likelihood shifts (ESM-1v)

Lower mutational tolerance supports functional importance.

---

# Pipeline Overview

```mermaid
flowchart LR
  A["Salt bridge residues (PDB)"]
  --> B["Step01 Map PDB → UniProt (SCOP alignment)"]
  --> C["Step02 Extract AlphaMissense scores"]
  --> D["Step03 Visualize conservation signal"]
  --> E["Step04 Prepare ESM mutation table"]
  --> F["Step05 Run ESM-1v ensemble"]
  --> G["Step06 Merge predictions"]
  --> H["Step07 Plot mutational constraint"]
```

---

# ⚡ Quick Start (Minimal Reproduction)

After downloading required datasets:

```bash
conda create -n saltbridge-ml python=3.10 -y
conda activate saltbridge-ml

pip install torch fair-esm pandas numpy matplotlib tqdm
```

Run the core pipeline:

```bash
python scripts/step01_map_pdb_to_uniprot.py ...
python scripts/step02_extract_alphamissense.py ...
python scripts/step05_run_esm_ensemble.py ...
```

Results will reproduce the mutation constraint analyses reported in the manuscript.

---

# Repository Structure

```text
ml_prediction_pipeline/
│
├── scripts/
│   ├── step01_map_pdb_to_uniprot.py
│   ├── step02_extract_alphamissense.py
│   ├── step04_prepare_esm_input.py
│   ├── step05_run_esm_ensemble.py
│   └── step06_merge_results.py
│
├── notebooks/
│   ├── step03_plot_alphamissense.ipynb
│   └── step07_plot_esm_results.ipynb
│
├── data/                 # NOT tracked in git
│   ├── mapping/
│   ├── alphamissense/
│   └── esm_input/
│
└── results/              # reproducible outputs
```

---

# Data Availability

This pipeline depends exclusively on **public datasets**.  
Due to licensing and file-size constraints, these datasets are **not redistributed**.

---

## SCOP Structural Classification

Used for mapping PDB residues to UniProt coordinates.

- Source:  
https://www.ebi.ac.uk/pdbe/scop/files/scop-cla-latest.txt

Recommended location:

```
data/scop/
```

⚠️ Always record the release version for reproducibility.

---

## AlphaMissense Variant Predictions

- Publication:  
**Cheng et al., Nature (2023)**

- Download:  
https://alphamissense.hegelab.org/

Required file:

```
AlphaMissense_aa_substitutions.tsv
```

⚠️ File size is large (tens of GB).  
High-speed storage recommended.

---

## ESM-1v Protein Language Model

- Official repository:  
https://github.com/facebookresearch/esm

Models used:

```
esm1v_t33_650M_UR90S_1–5
```

Install:

```bash
pip install fair-esm
```

Models download automatically.

⚠️ GPU strongly recommended.

---

# Reproducibility Statement

All results can be reproduced using:

1. Public SCOP classification  
2. AlphaMissense dataset  
3. Official ESM models  

No proprietary datasets are required.

Intermediate files are intentionally excluded from version control.

---

# Step-by-Step Usage

---

## Step01 — Map Salt-Bridge Residues to UniProt

Uses curated SCOP alignment rather than heuristic sequence matching.

```bash
python scripts/step01_map_pdb_to_uniprot.py \
  --saltbridge-file data/saltbridges/cla.txt \
  --scop-file data/scop/scop-cla-latest.txt \
  --output data/mapping/cla_mapping.csv
```

---

## Step02 — Extract AlphaMissense Scores

Only **alanine substitutions (X→A)** are analyzed for consistency across residue types.

```bash
python scripts/step02_extract_alphamissense.py \
  --mapping data/mapping/cla_mapping.csv \
  --am-tsv data/alphamissense/AlphaMissense_aa_substitutions.tsv \
  --output results/alphamissense/cla_scores.csv
```

---

## Step03 — Visual Inspection

Notebook:

```
notebooks/step03_plot_alphamissense.ipynb
```

Used to confirm mutational constraint patterns.

---

## Step04 — Prepare ESM Mutation Table

```bash
python scripts/step04_prepare_esm_input.py \
  --mapping data/mapping/cla_mapping.csv \
  --output data/esm_input/cla_esm.csv
```

---

## Step05 — Run ESM-1v Ensemble

Five independent models are averaged to reduce stochastic variance.

```bash
python scripts/step05_run_esm_ensemble.py \
  --model-location esm1v_t33_650M_UR90S_1 \
                   esm1v_t33_650M_UR90S_2 \
                   esm1v_t33_650M_UR90S_3 \
                   esm1v_t33_650M_UR90S_4 \
                   esm1v_t33_650M_UR90S_5 \
  --dms-input data/esm_input/cla_esm.csv \
  --dms-output results/esm/cla_scores.csv
```

---

## Step06 — Merge Predictions

```bash
python scripts/step06_merge_results.py \
  --mapping data/mapping/cla_mapping.csv \
  --esm-results results/esm/cla_scores.csv \
  --output results/final/cla_mutational_constraint.csv
```

---

# Key Methodological Choices

- Alanine scanning ensures uniform physicochemical perturbation.
- Ensemble inference improves robustness.
- SCOP-based mapping avoids alignment artifacts.

---

# Compute Requirements

Recommended:

- GPU with ≥16GB VRAM (for ESM)
- ≥64GB RAM if processing large AlphaMissense tables
- SSD / high-speed network storage

---

# Citation

If you use this pipeline, please cite:

- AlphaMissense (Nature 2023)  
- ESM-1v  
- SCOP  
- The associated manuscript  

---

# License

Recommended: MIT
