#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GROMACS Input File Generator for Boiling Point Simulations.

This script reads an Excel file containing SMILES and molecular weights,
and generates all necessary GROMACS input files (topologies, coordinates,
PackMol inputs) for each compound.

Force field: OpenFF 2.0.0 (Sage)
Charges: Gasteiger (via RDKit)
Target system size: ~15,000 atoms (packed in 1/8 box, then replicated 2x2x2)
"""

import os
import math
import pandas as pd
import traceback
import argparse
from pathlib import Path

from openff.toolkit.topology import Molecule
from openff.interchange.components.interchange import Interchange
from openff.toolkit.typing.engines.smirnoff import ForceField
from rdkit import Chem

def generate_gromacs_files(index, smiles, mw, output_dir):
    comp_name = f"COMP{index}"
    comp_dir = output_dir / comp_name

    try:
        print(f"--- Processing: {comp_name} | SMILES: {smiles} ---")
        comp_dir.mkdir(parents=True, exist_ok=True)

        print("   1. Creating base molecule and computing parameters...")
        rdkit_mol = Chem.MolFromSmiles(smiles)
        molecule = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        molecule.name = "MOL"
        molecule.generate_conformers(n_conformers=4)
        molecule.assign_partial_charges(partial_charge_method='gasteiger')
        forcefield = ForceField('openff-2.0.0.offxml', allow_cosmetic_attributes=True)

        target_atoms = 15000
        num_atoms_per_molecule = molecule.n_atoms
        total_molecules_needed = math.ceil(target_atoms / num_atoms_per_molecule)
        num_molecules_packmol = math.ceil(total_molecules_needed / 8)
        num_molecules_final = num_molecules_packmol * 8
        total_atoms_final = num_molecules_final * num_atoms_per_molecule
        
        print(f"      -> Atoms per molecule: {num_atoms_per_molecule}")
        print(f"      -> Target total molecules: ~{total_molecules_needed} (for ~15000 atoms)")
        print(f"      -> Molecules for PackMol: {num_molecules_packmol} (1/8 of system)")
        print(f"      -> Final molecules after 2x2x2 replication: {num_molecules_final}")
        print(f"      -> Final total atoms: {total_atoms_final}")

        rho_target = 1000  # kg/m^3
        Na = 6.022e23
        
        mass_packmol_kg = (num_molecules_packmol * mw) / (Na * 1000)
        volume_packmol_m3 = mass_packmol_kg / rho_target
        box_length_packmol_nm = (volume_packmol_m3)**(1/3) * 1e9
        box_length_final_nm = box_length_packmol_nm * 2

        print(f"      -> PackMol box length: {box_length_packmol_nm:.4f} nm")
        print(f"      -> Final box length: {box_length_final_nm:.4f} nm")

        interchange = Interchange.from_smirnoff(forcefield, molecule.to_topology(), 
                                                box=[6.0, 6.0, 6.0], 
                                                charge_from_molecules=[molecule])

        single_mol_top_path = comp_dir / "MOL.top"
        interchange.to_top(single_mol_top_path)
        lines = single_mol_top_path.read_text().split('\n')
        with open(single_mol_top_path, 'w') as f:
            for line in lines:
                if line.strip() in ['[ system ]', '[ molecules ]']:
                    break
                f.write(line + '\n')

        interchange.to_gro(comp_dir / "MOL.gro")

        with open(comp_dir / "topol_single.top", "w") as f:
            f.write(single_mol_top_path.read_text())
            f.write("[ system ]\nSingle Molecule\n[ molecules ]\nMOL \t 1\n")

        with open(comp_dir / "topol_packmol.top", "w") as f:
            f.write(single_mol_top_path.read_text())
            f.write(f"[ system ]\nPacked Box\n[ molecules ]\nMOL \t {num_molecules_packmol}\n")

        with open(comp_dir / "topol_final.top", "w") as f:
            f.write(single_mol_top_path.read_text())
            f.write(f"[ system ]\nFinal Production Box\n[ molecules ]\nMOL \t {num_molecules_final}\n")
            
        with open(comp_dir / "packing_params.txt", "w") as f:
            f.write(f"box_length_packmol {box_length_packmol_nm}\n")
            f.write(f"box_length_final {box_length_final_nm}\n")
            f.write(f"num_molecules_packmol {num_molecules_packmol}\n")
            f.write(f"num_molecules_final {num_molecules_final}\n")

        packmol_input = f"""tolerance 2.0
filetype pdb
output conf_packed.pdb

structure MOL_minimized.pdb
  number {num_molecules_packmol}
  inside box 0. 0. 0. {box_length_packmol_nm*10} {box_length_packmol_nm*10} {box_length_packmol_nm*10}
end structure
"""
        with open(comp_dir / "packmol.inp", "w") as f:
            f.write(packmol_input)

        print(f"   SUCCESS: Files for {comp_name} are ready.\n")
        return True

    except Exception as e:
        print(f"   ERROR processing {comp_name}: {e}")
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(description="GROMACS Base Files Generator")
    parser.add_argument("excel_file", type=Path, help="Input Excel file.")
    parser.add_argument("-o", "--output_dir", type=Path, default=Path("runsims"), 
                        help="Output directory.")
    parser.add_argument("--start", type=int, default=1,
                    help="Row to start from (1 = first data row after header).")
    args = parser.parse_args()
    
    df = pd.read_excel(args.excel_file)
    args.output_dir.mkdir(exist_ok=True)
    
    start_idx = args.start - 1
    for i, row in df.iloc[start_idx:].iterrows():
        if pd.isna(row.get('smiles')) or pd.isna(row.get('mw')):
            print(f"WARNING: Skipping row {i+1} due to missing data.")
            continue
        generate_gromacs_files(i + 1, row['smiles'], row['mw'], args.output_dir)

    print("\nDONE! Folders are ready in the output directory.")

if __name__ == "__main__":
    main()

