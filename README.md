# final_year_dissertation_code_files
All code used in project: Computational prediction of non-antibiotic molecules that bind penicillin binding proteins and beta lactamases. 
Code was written, edited and commented with the assistance of Claude Opus 4.7
## Repository contents

| Folder | Description |
|---|---|
| `library_preparation/` | Python script for fetching compounds from the ChEMBL API and Bash script for converting them to PDBQT format |
| `docking/` | SLURM array job script for running AutoDock Vina on the University of Nottingham Ada HPC cluster |
| `analysis/` | Python script for collecting and ranking docking results, R scripts for RMSD calculation and figure generation |
| `data/` | Final ranked results from both target screens |

## Software dependencies

- AutoDock Vina 1.2.5 (build 23d1252-mod, ADFRsuite 1.1dev)
- Open Babel 2.4.1
- Meeko 0.6.1
- PyMOL 3.1.0 (open source)
- Python 3.10.18
- R 4.5.1 with ggplot2 (v4.0.2)

## Pipeline overview

1. `library_preparation/chembl_fetch_v3.py` — query ChEMBL REST API for compounds in the molecular weight window 467–567 Da; output as SDF
2. `library_preparation/prep_ligands.sh` — split SDF, convert to PDBQT via Open Babel and Meeko
3. `docking/slurm_dock.sh` — submit one Vina docking task per ligand on the HPC cluster
4. `analysis/collect_results.py` — parse Vina logs, rank by best-pose ΔG, output ranked CSV
5. `analysis/*.R` — generate figures and compute RMSDs

## Reproducibility notes

- Receptor preparation, redocking validation, and PyMOL-based pose analysis steps are described in the dissertation methods section but were performed interactively rather than scripted
- PDB structures used: 4KQO, 5OIZ, 3Q07, 5DF9 (downloadable from rcsb.org)

## Author

Adit Kadkol, BSc Biotechnology, University of Nottingham, 2026
