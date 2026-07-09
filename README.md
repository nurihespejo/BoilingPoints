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
