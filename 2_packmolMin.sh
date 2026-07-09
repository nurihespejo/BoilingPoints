#!/bin/bash
# =============================================================================
# Pack molecules into a simulation box and minimize the packed system
# =============================================================================

set -e

echo "=== Step 2: PackMol Packing & Box Minimization ==="

for dir in runsims/COMP*/; do
    [ ! -d "$dir" ] && continue
    comp=$(basename "$dir")
    echo "--- Processing: ${comp} ---"

    cd "$dir"

    if [ ! -f "MOL_minimized.pdb" ]; then
        echo "ERROR: MOL_minimized.pdb not found for ${comp}."
        cd ../..; continue
    fi

    # 1. Run PackMol
    packmol < packmol.inp
    if [ ! -f "conf_packed.pdb" ]; then 
        echo "ERROR: PackMol failed for ${comp}."
        cd ../..; continue
    fi

    # 2. Prepare minimization directory
    rm -rf packed_minim
    mkdir -p packed_minim
    cd packed_minim

    cp ../topol_packmol.top topol.top

    # 3. Convert PDB to GRO with correct box size
    box_len=$(grep "box_length_packmol" ../packing_params.txt | awk '{print $2}')
    gmx editconf -f ../conf_packed.pdb -o conf.gro -box "$box_len" "$box_len" "$box_len" -c

    # 4. Minimization
    gmx grompp -f ../../../minim.mdp -c conf.gro -p topol.top -o topol.tpr -maxwarn 5
    if [ ! -f "topol.tpr" ]; then 
        echo "ERROR: grompp failed for ${comp}"
        cd ../../..; continue
    fi

    gmx mdrun -deffnm topol

    if [ -f "confout.gro" ]; then
        echo "SUCCESS: Packed box minimized for ${comp}"
    else
        echo "ERROR: Minimization failed for ${comp}"
    fi

    cd ../../..
done

