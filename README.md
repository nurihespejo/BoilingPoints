# Molecular Dynamics Workflow for Boiling Point Prediction

This repository contains the full computational pipeline used to generate, minimize, and run molecular dynamics simulations in OpenFF/GROMACS workflow for boiling point prediction. The OPLS4 workflow reported in the paper was run using Schrödinger's Materials Science Suite (Desmond), which requires a commercial license. The same ML training pipeline and configuration apply; only the force field and simulation engine used to generate the input descriptors differ.

## Prerequisites
To reproduce this workflow, you need the following software installed:
- **Python 3** (with `pandas`, `rdkit`, `openff-toolkit`)
- **GROMACS** (version 2025.2 or compatible)
- **PackMol**

## Workflow Pipeline

The simulation setup is fully automated and divided into sequential scripts. 

### 1. System Setup (0_initial.py)
Reads the input Excel file (`TrainingDataset.xlsx`) containing SMILES and molecular weights. It uses the OpenFF 2.0.0 force field and Gasteiger charges to generate GROMACS topologies. It calculates the required box size for ~15,000 atoms, targeting an initial 1/8th box for PackMol.
### 2. Single Molecule Minimization (1_minimSingle.sh)
Iterates through all generated compounds, adds a 3.0 nm buffer box, and performs an energy minimization on the isolated molecule to prepare it for packing.
### 3. Box Packing & Minimization (2_packmolMin.sh)
Uses PackMol to pack the minimized molecule into the calculated 1/8th simulation box. Afterward, it performs a second energy minimization on the packed system to resolve any steric clashes.
### 4. Production Setup (3_production.sh)
Takes the minimized packed box and replicates it 2x2x2 using gmx genconf to reach the final ~15,000 atom target. It creates index groups splitting the molecules into two halves, and sets up individual production directories for multiple temperatures (300K, 400K, 500K), automatically adjusting the sim.mdp parameters.
### 5. Trajectory Analysis (4_analysis.sh)
Once the production simulations are finished, this script extracts the relevant physical properties (Coulomb interactions, Lennard-Jones interactions, Density, Volume, and Heat Capacity) from the GROMACS .edr files and compiles them into a single analysis_results.csv file, ready for Machine Learning training.

## Machine Learning Model

The `ML_Model/` directory contains the trained CatBoost regression model and the complete training pipeline built with the [MotherML](https://github.com/Bayer-Group/MotherML) framework.

We provide the OpenFF-based model so that the entire workflow—from MD simulation to boiling-point prediction—can be reproduced using open-source tools (GROMACS + OpenFF + MotherML). The OPLS4-based model, which yields the primary results reported in the paper, was trained with the same pipeline and configuration; the only difference is the force field used to generate the input descriptors.

### Contents

| File | Description |
|------|-------------|
| `inference_model.pkl.gz` | Trained model (OpenFF descriptors, 6 features selected by CatBoost) |
| `train.py` | Training script |
| `config/parameters.yaml` | Feature columns, target, CV and feature-selection settings |
| `config/model/model.yaml` | CatBoost hyperparameters |
| `config/optuna/optuna.yaml` | Bayesian hyperparameter optimization (50 trials, 4-fold outer / 3-fold inner CV) |

### Inference

To predict boiling points for new compounds:

1. Run the OpenFF simulation pipeline (scripts `0_initial.py` through `4_analysis.sh`) at 300, 400, and 500 K.
2. Collect the thermodynamic descriptors into a CSV with the column names defined in `config/parameters.yaml`.
3. Run inference:

```python
import mother.inference as infer

model = infer.load_model("ML_Model/inference_model.pkl.gz")
predictions = infer.predict(model, your_dataframe)
```

See the [MotherML documentation](https://github.com/Bayer-Group/MotherML) for details on the inference API.

### Retraining

To retrain the model with your own data:

```bash
pip install motherml
cd ML_Model
# Place your training data in data/compounds.csv
python train.py
```

### Note on OPLS4

The OPLS4 simulations reported in the paper were performed using the Schrödinger Materials Science Suite (Desmond), which requires a commercial license. The same training pipeline and configuration apply; only the input descriptors differ.


