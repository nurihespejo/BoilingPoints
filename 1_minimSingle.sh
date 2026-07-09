#!/bin/bash
# =============================================================================
# Single-molecule energy minimization
# =============================================================================

set -e

echo "=== Step 1: Single Molecule Minimization ==="

for dir in runsims/COMP*/; do
    [ ! -d "$dir" ] && continue
    comp=$(basename "$dir")
    echo "--- Processing: ${comp} ---"

    rm -rf "${dir}/single_mol_minim"
    mkdir -p "${dir}/single_mol_minim"

    cp "${dir}/MOL.gro" "${dir}/single_mol_minim/conf.gro"
    cp "${dir}/topol_single.top" "${dir}/single_mol_minim/topol.top"

    cd "${dir}/single_mol_minim"

    # 1. Create a large box to avoid periodic interactions
    gmx editconf -f conf.gro -o conf.gro -d 3.0 -c

    # 2. Minimization
    gmx grompp -f ../../../minim.mdp -c conf.gro -p topol.top -o topol.tpr -maxwarn 5
    if [ ! -f "topol.tpr" ]; then 
        echo "ERROR: grompp failed for ${comp}"
        cd ../../..; continue
    fi
    
    gmx mdrun -deffnm topol

    # 3. Create final PDB for PackMol
    if [ -f "confout.gro" ]; then
        echo "0" | gmx trjconv -s topol.tpr -f confout.gro -o ../MOL_minimized.pdb -pbc none
        echo "SUCCESS: MOL_minimized.pdb created for ${comp}"
    else
        echo "ERROR: confout.gro not found for ${comp}"
    fi

    cd ../../..
done

