#!/bin/bash
# =============================================================================
# Replicate the system 2x2x2, create index groups, and setup temperatures
# =============================================================================

T_INIT=300
T_STEP=100
T_FINAL=500

echo "=== Step 3: Production Setup ==="

for dir in runsims/COMP*/; do
    [ ! -d "$dir" ] && continue
    comp=$(basename "$dir")
    echo "--- Processing: ${comp} ---"
    
    cd "$dir"

    if [ ! -f "packed_minim/confout.gro" ]; then
        echo "ERROR: packed_minim/confout.gro not found. Skipping..."
        cd ../../; continue
    fi

    # --- Base Setup (2x2x2 Replication) ---
    if [ ! -d "prod_base" ]; then
        mkdir -p prod_base
        cd prod_base

        echo "   Replicating box 2x2x2..."
        gmx genconf -f ../packed_minim/confout.gro -o conf.gro -nbox 2 2 2

        cp ../topol_final.top topol.top

        num_mols=$(grep "MOL" topol.top | tail -n 1 | awk '{print $2}')
        half_mols=$((num_mols / 2))

        echo "   Creating index groups..."
        gmx make_ndx -f conf.gro << EOF
r 1-${half_mols}
r $((half_mols + 1))-${num_mols}
name 3 1st
name 4 2nd
q
EOF
        cd ../
    fi

    # --- Temperature Folders Setup ---
    for T in $(seq $T_INIT $T_STEP $T_FINAL); do
        echo "   -> Setting up production at ${T}K"
        
        rm -rf "prod_${T}K"
        cp -r prod_base "prod_${T}K"
        cd "prod_${T}K"

        # Copy and edit grompp.mdp
        cp ../../../../grompp.mdp grompp.mdp
        sed -i "s/\(ref-t\s*=\s*\)[0-9.]*/\1$T.00/" grompp.mdp
        sed -i "s/\(gen_temp\s*=\s*\)[0-9.]*/\1$T/" grompp.mdp

        # Generate .tpr for production
        gmx grompp -n index.ndx -maxwarn 10 -o topol.tpr > grompp.log 2>&1
        
        if [ ! -f "topol.tpr" ]; then 
            echo "      ERROR: grompp failed for ${T}K"
        else
            echo "      SUCCESS: topol.tpr generated for ${T}K"
        fi
        
        cd ../
    done
    cd ../../
done

echo "=== Production Setup Complete ==="

