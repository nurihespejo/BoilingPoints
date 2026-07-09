#!/bin/bash
# =============================================================================
# Extract physical properties from GROMACS energy files (.edr)
# =============================================================================

set -e

echo "=== Step 4: Trajectory Analysis ==="

START_TIME=15000   # Analysis window start time in ps
OUTCSV="analysis_results.csv"

# Initialize CSV file
echo "compound,temp,property,value" > "$OUTCSV"

# Loop through all production directories
for dir in runsims/COMP*/prod_*K/; do
    [ ! -d "$dir" ] && continue
    
    COMP=$(basename "$(dirname "$dir")")
    TEMP=$(basename "$dir" | sed 's/prod_//; s/K//')
    
    # Locate the energy file (can be ener.edr or prod.edr depending on run)
    if [ -f "$dir/ener.edr" ]; then
        EDR_FILE="$dir/ener.edr"
    elif [ -f "$dir/prod.edr" ]; then
        EDR_FILE="$dir/prod.edr"
    else
        echo "   WARNING: No .edr file found in $dir. Skipping..."
        continue
    fi

    echo "--- Analyzing $COMP at ${TEMP}K ---"

    # --- 1. BASIC PROPERTIES ---
    declare -a PROPS=("Coul-SR:1st-2nd" "LJ-SR:1st-2nd" "Density" "Volume")

    for PROP in "${PROPS[@]}"; do
        VAL=$(echo "$PROP" | gmx energy -f "$EDR_FILE" -b $START_TIME 2>/dev/null | awk -v p="$PROP" '$1==p {print $2}')
        [[ -z "$VAL" ]] && VAL="NaN"
        echo "$COMP,$TEMP,$PROP,$VAL" >> "$OUTCSV"
    done

    # --- 2. HEAT CAPACITY (Cp) EXTRACTION ---
    # Find the dynamic indices for Enthalpy and Temperature in the energy file
    ENTHALPY_NUM=$(gmx energy -f "$EDR_FILE" < /dev/null 2>&1 | awk '{for(i=1;i<=NF;i++) if($i=="Enthalpy") print $(i-1)}')
    TEMP_NUM=$(gmx energy -f "$EDR_FILE" < /dev/null 2>&1 | awk '{for(i=1;i<=NF;i++) if($i=="Temperature") print $(i-1)}')

    if [[ -z "$ENTHALPY_NUM" || -z "$TEMP_NUM" ]]; then
        CP_VAL="NaN"
    else
        CP_OUT=$(printf "%s\n%s\n0\n" "$ENTHALPY_NUM" "$TEMP_NUM" | gmx energy -f "$EDR_FILE" -b $START_TIME -fluct_props 2>&1)
        CP_VAL=$(echo "$CP_OUT" | awk '/Heat capacity at constant pressure Cp/ {print $8}')
        [[ -z "$CP_VAL" ]] && CP_VAL="NaN"
    fi

    echo "$COMP,$TEMP,Cp,$CP_VAL" >> "$OUTCSV"
done

echo "=== Analysis Complete. Results saved to $OUTCSV ==="

