# Molecular Dynamics Workflow for Boiling Point Prediction

This repository contains the full computational pipeline used to generate, minimize, and run molecular dynamics simulations in GROMACS for boiling point prediction.

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

## Machine Learning Model (OPLS4 MD-only)

The trained CatBoost regression model and the complete training pipeline are provided in the `ML_Model/` directory.

### Contents

- `inference_model.pkl.gz` — Trained model for boiling-point prediction from three MD-derived descriptors: ΔH_vap (300 K), Hildebrand solubility parameter (300 K), and C_P (500 K).
- `train.py` — Training script using the [MotherML](https://github.com/Bayer-Group/MotherML) framework.
- `config/` — Configuration files:
  - `parameters.yaml` — Feature columns, target column, cross-validation and feature-selection settings.
  - `model/model.yaml` — CatBoost model hyperparameters.
  - `optuna/optuna.yaml` — Bayesian hyperparameter optimization settings (50 trials, 4-fold outer / 3-fold inner CV).
- `data/compounds.csv` — Training dataset (1,280 compounds) with MD-derived descriptors at 300, 400, and 500 K from OPLS4 simulations.

### Inference on new compounds

To predict boiling points for new compounds, you need to:

1. Run NPT simulations at 300, 400, and 500 K for each compound (see OpenFF pipeline above or use OPLS4/Desmond).
2. Extract ΔH_vap (300 K), Hildebrand solubility parameter (300 K), and C_P (500 K).
3. Load the model and predict:

```python
import mother.inference as infer

model = infer.load_model("ML_Model/inference_model.pkl.gz")
predictions = infer.predict(model, your_dataframe)
```

See [MotherML documentation](https://github.com/Bayer-Group/MotherML) for details on the inference API.

### Retraining

```bash
pip install motherml
cd ML_Model
python train.py
```

### OPLS4 simulations

The OPLS4 simulations were performed using the Schrödinger Materials Science Suite (Desmond), which requires a commercial license. The OPLS4 input files and force-field parameters are proprietary and cannot be redistributed. The OpenFF pipeline provided in this repository produces equivalent descriptors using open-source tools.

