#!/usr/bin/env python3
"""
Step 5: Compute derived thermodynamic descriptors from raw GROMACS analysis.

This script reads the raw energy terms extracted by 4_analysis.sh and computes
the derived descriptors used by the ML model: cohesive energy, heat of
vaporization, Hildebrand solubility parameter, and normalized heat capacity.

This script is part of the OpenFF/GROMACS workflow. The OPLS4 workflow uses
Schrödinger Desmond (commercial) and is not covered by these open-source scripts.

Usage:
    python 5_createExcel.py
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ===================== CONFIG =====================
SYSTEMS = ["G1_FULL", "G234"]
BASE = Path.cwd()
INPUT_CSV_NAME = "analysis_results.csv"
TEMPS = ["300K", "400K", "500K"]
PROPS_RAW = ["Coul-SR:1st-2nd", "LJ-SR:1st-2nd", "Density", "Volume", "Cp"]

R_J_MOL_K = 8.314       # Gas constant in J/(mol·K)
KJ_TO_KCAL = 1 / 4.184
J_TO_KCAL = 1 / 4184.0
NM3_TO_A3 = 1000.0


def read_nmols(topol_path):
    """Read number of molecules from a GROMACS topology file."""
    if not topol_path.exists():
        return np.nan
    with open(topol_path) as f:
        lines = f.readlines()
    in_block = False
    for line in lines:
        line = line.strip()
        if line.startswith("[ molecules ]"):
            in_block = True
            continue
        if in_block:
            if not line or line.startswith(";"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    return int(parts[1])
                except ValueError:
                    return np.nan
    return np.nan


def read_molar_mass(topol_path):
    """Read molar mass from the [ atoms ] section of a GROMACS topology file."""
    if not topol_path.exists():
        return np.nan
    with open(topol_path) as f:
        lines = f.readlines()
    in_atoms, total_mass = False, 0.0
    for line in lines:
        line = line.strip()
        if line.startswith("[ atoms ]"):
            in_atoms = True
            continue
        if in_atoms:
            if line.startswith("["):
                break
            if not line or line.startswith(";"):
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    total_mass += float(parts[7])
                except (ValueError, IndexError):
                    continue
    return total_mass if total_mass > 0 else np.nan


def process_system(system):
    print(f"\nProcessing {system}")
    csv_file = BASE / f"{system}_{INPUT_CSV_NAME}"
    if not csv_file.exists():
        print(f"  Not found: {csv_file}")
        return

    df = pd.read_csv(csv_file)
    if df.empty:
        print(f"  Empty file: {csv_file}")
        return

    df["comp_num"] = pd.to_numeric(
        df["compound"].str.extract(r"(\d+)", expand=False), errors="coerce"
    )
    df.dropna(subset=["comp_num"], inplace=True)
    df["comp_num"] = df["comp_num"].astype(int)
    df["Temp"] = df["temp"].str.replace("prod_", "", regex=False)

    existing = dict(zip(df["comp_num"], df["compound"]))
    if not existing:
        print(f"  No valid compounds in {csv_file}")
        return

    max_comp = df["comp_num"].max()
    print(f"  Compounds {df['comp_num'].min()} to {max_comp}")

    rows = []
    for comp_num in range(1, max_comp + 1):
        comp_name = f"COMP{comp_num:03d}"
        row = {"Compound": comp_name}

        if comp_num not in existing:
            row.update({"nmols": np.nan, "M_gmol": np.nan,
                        "Normalization_Factor": np.nan})
            for T in TEMPS:
                for prop in PROPS_RAW:
                    row[f"{prop}_{T}"] = np.nan
                for col in ["CohE", "HeatVap", "V_A3_per_molecule",
                            "SolubilityParam_cal_cm3_0.5", "Solubility",
                            "Cp_kcal_molK", "Cp"]:
                    row[f"{col}_{T}"] = np.nan
            rows.append(row)
            continue

        comp = existing[comp_num]
        topol = BASE / system / "runsims" / comp / "prod_300K" / "topol.top"
        nmols = read_nmols(topol)
        M = read_molar_mass(topol)
        row.update({"nmols": nmols, "M_gmol": M})

        # Normalization: accounts for the two-group energy decomposition
        norm_factor = (
            (2 * (nmols - 1)) / (nmols ** 2)
            if not np.isnan(nmols) and nmols > 1
            else np.nan
        )
        row["Normalization_Factor"] = norm_factor

        for T in TEMPS:
            T_K = int(T.replace("K", ""))

            # Extract raw values for this compound and temperature
            raw_vals = {}
            for prop in PROPS_RAW:
                subset = df[
                    (df["comp_num"] == comp_num)
                    & (df["Temp"] == T)
                    & (df["property"] == prop)
                ]
                raw_vals[prop] = (
                    subset["value"].values[0] if not subset.empty else np.nan
                )

            for prop in PROPS_RAW:
                row[f"{prop}_{T}"] = raw_vals[prop]

            # --- Cohesive energy (kJ/mol per molecule) ---
            e_coul = raw_vals.get("Coul-SR:1st-2nd", np.nan)
            e_lj = raw_vals.get("LJ-SR:1st-2nd", np.nan)
            e_inter = (
                e_coul + e_lj
                if not np.isnan(e_coul) and not np.isnan(e_lj)
                else np.nan
            )
            e_coh_kj = (
                abs(e_inter * norm_factor)
                if not np.isnan(e_inter) and not np.isnan(norm_factor)
                else np.nan
            )

            # --- Heat of vaporization: CohE + RT ---
            rt_kj = (R_J_MOL_K * T_K) / 1000.0
            heat_vap_kj = (
                e_coh_kj + rt_kj if not np.isnan(e_coh_kj) else np.nan
            )

            row[f"CohE_{T}"] = (
                e_coh_kj * KJ_TO_KCAL if not np.isnan(e_coh_kj) else np.nan
            )
            row[f"HeatVap_{T}"] = (
                heat_vap_kj * KJ_TO_KCAL if not np.isnan(heat_vap_kj) else np.nan
            )

            # --- Molar volume per molecule (A^3) ---
            volume_nm3 = raw_vals.get("Volume", np.nan)
            V_A3 = (
                (volume_nm3 * NM3_TO_A3) / nmols
                if not np.isnan(volume_nm3) and nmols > 0
                else np.nan
            )
            row[f"V_A3_per_molecule_{T}"] = V_A3

            # --- Hildebrand solubility parameter: delta = sqrt(CohE / Vm) ---
            density_g_cm3 = raw_vals.get("Density", np.nan) / 1000.0
            Vm_cm3_mol = (
                M / density_g_cm3
                if not np.isnan(density_g_cm3)
                and density_g_cm3 > 0
                and not np.isnan(M)
                else np.nan
            )
            coh_kcal = row[f"CohE_{T}"]
            if (
                not np.isnan(coh_kcal)
                and not np.isnan(Vm_cm3_mol)
                and Vm_cm3_mol > 0
                and coh_kcal > 0
            ):
                # delta in (cal/cm^3)^0.5
                delta_cal = np.sqrt((coh_kcal * 1000) / Vm_cm3_mol)
                # delta in (J/cm^3)^0.5 = (MPa)^0.5
                delta_mpa = delta_cal * np.sqrt(4.184)
            else:
                delta_cal = np.nan
                delta_mpa = np.nan

            row[f"SolubilityParam_cal_cm3_0.5_{T}"] = delta_cal
            row[f"Solubility{T}"] = delta_mpa

            # --- Heat capacity: normalize by nmols ---
            cp_raw = raw_vals.get("Cp", np.nan)
            if not np.isnan(cp_raw) and not np.isnan(nmols) and nmols > 0:
                cp_J_molK = cp_raw / nmols
            else:
                cp_J_molK = np.nan

            row[f"Cp_kcal_molK_{T}"] = (
                cp_J_molK * J_TO_KCAL if not np.isnan(cp_J_molK) else np.nan
            )
            row[f"Cp_{T}"] = cp_J_molK

        rows.append(row)

    # --- Build final DataFrame with ordered columns ---
    column_order = ["Compound", "nmols", "M_gmol", "Normalization_Factor"]
    for T in TEMPS:
        for prop in PROPS_RAW:
            column_order.append(f"{prop}_{T}")
        column_order.extend([
            f"CohE_{T}",
            f"HeatVap_{T}",
            f"V_A3_per_molecule_{T}",
            f"SolubilityParam_cal_cm3_0.5_{T}",
            f"Solubility{T}",
            f"Cp_kcal_molK_{T}",
            f"Cp_{T}",
        ])

    final = pd.DataFrame(rows).reindex(columns=column_order)

    out_csv = BASE / f"{system}_descriptors.csv"
    final.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"  Results written to {out_csv}")


def main():
    for system in SYSTEMS:
        process_system(system)


if __name__ == "__main__":
    main()
