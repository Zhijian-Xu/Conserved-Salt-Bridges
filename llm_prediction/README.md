# ML Prediction Pipeline  
## AlphaMissense + ESM-1v Analysis of Salt-Bridge Residues

This module evaluates mutational sensitivity of spatially conserved salt-bridge residues using:

- **AlphaMissense**
- **ESM-1v ensemble**

The pipeline maps salt-bridge residues from **PDB positions → UniProt positions**, retrieves mutation scores, and validates conservation signals using deep-learning predictors.

---

# Overview

```mermaid
flowchart LR
  A["Salt bridge residues (PDB)"]
  --> B["Step01 Map PDB → UniProt"]
  --> C["Step02 Extract AlphaMissense ALA scores"]
  --> D["Step03 Plot AlphaMissense results"]
  --> E["Step04 Prepare ESM mutation input"]
  --> F["Step05 Run ESM-1v ensemble"]
  --> G["Step06 Merge & visualize results"]
```

---

# Installation

Recommended environment:

```bash
conda create -n saltbridge-ml python=3.10 -y
conda activate saltbridge-ml
conda install pandas numpy matplotlib tqdm -y
pip install torch
pip install fair-esm
```

GPU strongly recommended for ESM inference.

---

# Required external data

## 1) SCOP mapping file

Download:

https://www.ebi.ac.uk/pdbe/scop/files/scop-cla-latest.txt

Used for PDB → UniProt residue alignment.

---

## 2) AlphaMissense dataset

Download from:

https://alphamissense.hegelab.org/

Example file:

```
AlphaMissense_aa_substitutions.tsv
```

---

# Directory structure

```text
ml_prediction_pipeline/
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
├── data/
│   ├── mapping/
│   ├── alphamissense/
│   └── esm_input/
│
└── results/
    ├── alphamissense/
    └── esm/
```

---

# Step-by-step usage

---

## Step01 — Map PDB residues to UniProt

```bash
python scripts/step01_map_pdb_to_uniprot.py \
  --saltbridge-file data/saltbridges/cla.txt \
  --scop-file data/scop/scop-cla-latest.txt \
  --output data/mapping/cla_mapping.csv
```

Output:

```
mapping_*.csv
```

---

## Step02 — Extract AlphaMissense ALA mutation scores

```bash
python scripts/step02_extract_alphamissense.py \
  --mapping data/mapping/cla_mapping.csv \
  --am-tsv data/alphamissense/AlphaMissense_aa_substitutions.tsv \
  --output results/alphamissense/cla_AM_scores.csv
```

Only ALA substitutions (X→A) are retained.

---

## Step03 — Plot AlphaMissense scores

Notebook:

```
notebooks/step03_plot_alphamissense.ipynb
```

---

## Step04 — Prepare ESM mutation input

```bash
python scripts/step04_prepare_esm_input.py \
  --mapping data/mapping/cla_mapping.csv \
  --output data/esm_input/cla_ESM_input.csv
```

---

## Step05 — Run ESM-1v ensemble prediction

```bash
python scripts/step05_run_esm_ensemble.py \
  --model-location esm1v_t33_650M_UR90S_1 \
                   esm1v_t33_650M_UR90S_2 \
                   esm1v_t33_650M_UR90S_3 \
                   esm1v_t33_650M_UR90S_4 \
                   esm1v_t33_650M_UR90S_5 \
  --dms-input data/esm_input/cla_ESM_input.csv \
  --dms-output results/esm/cla_esm_scores.csv
```

Five models are averaged for robustness.

---

## Step06 — Merge predictions

```bash
python scripts/step06_merge_results.py \
  --mapping data/mapping/cla_mapping.csv \
  --esm-results results/esm/cla_esm_scores.csv \
  --output results/esm/cla_final.csv
```

---

# Key Design Choices

- Only **ALA substitutions** are analyzed for consistency.
- Ensemble ESM prediction improves stability.
- Mapping relies on curated SCOP alignment.

---

# Citation

Please cite:

- AlphaMissense  
- ESM-1v  
- SCOP  
- Associated manuscript

