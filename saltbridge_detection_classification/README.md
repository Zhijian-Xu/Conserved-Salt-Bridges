# SCOP Salt Bridge Identification & Classification

A reproducible pipeline to **split SCOP domains**, **align domains within a family**, **detect salt bridges**, and classify them into **conserved / unconserved**, and further into **classical conserved / non-classical conserved** salt bridges.

---

## Pipeline overview

```mermaid
flowchart LR
  A[SCOP family list (*.txt)] --> B[Step01 Split domains]
  B --> C[Step02 Align domains & detect salt bridges]
  C --> D[Step03 Spatial clustering: conserved vs unconserved]
  D --> E[Step04 CA distance + geometry (cla/nocla/other)]
  E --> F[Step05 Distance statistics & plots]
  F --> G[Step06 Classical vs non-classical classification]
  G --> H[Step07 Finalize true conserved: recurrence >= 5]
Key thresholds (as used in this project)
Salt bridge detection: Asp/Glu ↔ Lys/Arg within 5.0 Å (see Step02 output)
Spatial conservation (Step03): bridge center distance ≤ 2.0 Å
True conservation (Step07): same-position salt bridge observed ≥ 5 times


Repository structure (recommended)
scripts/
  step01_split_scop_domains.py
  step02_align_domains_detect_saltbridges.py
  step03_cluster_conserved_vs_unconserved.ipynb
  step04_extract_ca_coordinates.py
  step05_distance_statistics_visualization.ipynb 
  step06_classify_classical_vs_nonclassical.ipynb
  step07_finalize_true_conserved.py


data/
  raw_families/                 # input raw structures grouped by family
  family_lists/                 # <family_id>.txt domain lists + reference.txt
  domains/                      # Step01 output: <family_id>/*.cif


results/
  step02_saltbridges_raw/       # Step02 output: per-family txt
  step03_clusters/              # Step03 output (user-defined)
  step04_ca_distances/          # Step04 output
  figures/                      # Step05 output
  final_tables/                 # Step06 output
  step07_true_conserved_k5/     # Step07 output: final filtered txt


Installation
This project relies on PyMOL Python API.
Recommended (conda-forge):
conda create -n saltbridge python=3.10 -y
conda activate saltbridge
conda install -c conda-forge pymol-open-source numpy pandas matplotlib jupyterlab -y


Input formats
1) Family domain list: data/family_lists/<family_id>.txt
Each line defines a domain to extract. Supported formats:
3 columns: PDBID CHAIN RES_RANGE
1abc A 10-120
4 columns (discontinuous segments): PDBID CHAIN RES_RANGE1 RES_RANGE2
1abc A 10-50 80-120
5 columns (two-chain domain): PDBID CHAIN1 RES_RANGE1 CHAIN2 RES_RANGE2
1abc A 10-50 B 1-60

2) Domain structures layout
Raw structures (input to Step01):
data/raw_families/<family_id>/*.cif.gz
Extracted domains (output of Step01, input to Step02/04):
data/domains/<family_id>/<domain_name>.cif

3) Reference list (for Step04)
data/family_lists/reference.txt format:
<family_id>  <reference_domain_name>
Example:
4007548  3xyz-A10_120


Step-by-step usage
Step01 — Split SCOP domains into domain structures
Script: scripts/step01_split_scop_domains.py
python scripts/step01_split_scop_domains.py \
  --input data/raw_families \
  --list-dir data/family_lists \
  --output data/domains
Output:data/domains/<family_id>/*.cif

Step02 — Align domains within each family & detect salt bridges
Script: scripts/step02_align_domains_detect_saltbridges.py
Note: this step requires com.py/com2.py utilities (place them in scripts/utils/).
python scripts/step02_align_domains_detect_saltbridges.py \
  --list-dir data/family_lists \
  --domain-dir data/domains \
  --output results/step02_saltbridges_raw \
  --utils-dir scripts/utils
Output:results/step02_saltbridges_raw/<family_id>.txt
Each salt bridge line contains fields like:
pi_chain/pi_resn/pi_resi
cation_chain/cation_resn/cation_resi
distance
bridge_com: [x,y,z]

Step03 — Cluster conserved vs unconserved (spatial)
Notebook: notebooks/step03_cluster_conserved_vs_unconserved.ipynb
Definition used:bridge center distance ≤ 2 Å → conserved
Input (recommended):results/step02_saltbridges_raw/<family_id>.txt
Output:user-defined (recommended): results/step03_clusters/<family_id>.txt

Step04 — Extract Cα distances & geometry category (cla/nocla/other/error)
Script: scripts/step04_extract_ca_coordinates.py
python scripts/step04_extract_ca_coordinates.py \
  -f <family_id> \
  -i results/step03_clusters \
  -o results/step04_ca_distances \
  -r data/family_lists/reference.txt \
  -p data/domains
Optional: generate PyMOL sessions for visualization:
python scripts/step04_extract_ca_coordinates.py \
  -f <family_id> \
  -i results/step03_clusters \
  -o results/step04_ca_distances \
  -r data/family_lists/reference.txt \
  -p data/domains \
  --pse --pse-dir results/step04_ca_distances/pse
Output:results/step04_ca_distances/{all,cla,nocla,other,error}/<family_id>.txt

Step05 — Distance statistics & plots
Notebook: notebooks/step05_distance_statistics_visualization.ipynb
Input:results/step04_ca_distances/*/<family_id>.txt
Output (recommended):results/figures/

Step06 — Classify classical vs non-classical conserved salt bridges
Notebook: notebooks/step06_classify_classical_vs_nonclassical.ipynb
Output (recommended):
results/final_tables/classical conserved
results/final_tables/non-classical conserved

Step07 — Finalize “true conserved” salt bridges (recurrence ≥ 5)
Script: scripts/step07_finalize_true_conserved.py
This is the final post-processing step. A salt bridge is considered true conserved if it recurs ≥ 5 times in the family-level results.
python scripts/step07_finalize_true_conserved.py results/final_tables \
  --outdir results/step07_true_conserved_k5 \
  --min_occurrence 5 \
  --key-mode residue
--key-mode:
line: strict (entire line must match)
residue: residue-identity mode (chain + resname + resi), robust to numeric variations
Output:results/step07_true_conserved_k5/<family_id>.txt

