import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px

# constants

R = 8.314  # universal gas constant, J/(mol·K)
T0 = 298.15  # reference temperature (25 °C) in Kelvin

# Gas properties: H0 (Henry's constant at T0), delta_H_sol (enthalpy of solution),
# k_s (salinity correction factor), molar_mass (g/mol)
gas_data = {
    'methane': {'H0': 1400.0, 'delta_H_sol': 8300.0, 'k_s': 0.12, 'molar_mass': 16.04},
    'co2':     {'H0': 29.4,   'delta_H_sol': 19000.0, 'k_s': 0.10, 'molar_mass': 44.01},
    'nitrogen':{'H0': 1600.0, 'delta_H_sol': 5800.0,  'k_s': 0.13, 'molar_mass': 28.02},
    'air':     {'H0': 770.0,  'delta_H_sol': 10000.0, 'k_s': 0.12, 'molar_mass': 28.97}
}

def fahrenheit_to_celsius(F):
    return (F - 32.0) * (5.0 / 9.0)

def barg_to_atm(barg):
    """
    Convert gauge pressure (bar gauge) to absolute pressure in atm.
    """
    # 1 bar gauge = (barg + 1) bar absolute; 1 bar = 0.986923 atm
    return (barg + 1.0) * 0.986923

def psi_to_atm(psi):
    """Convert pressure from psi to atm."""
    return psi / 14.696

# Salinity conversion (mg/L or ppt → mg/L)
def salinity_to_mg_L(value, unit):
    """Convert salinity input from mg/L or ppt to mg/L."""
    if unit == 'ppt':
        # 1 ppt ≈ 1 g/L = 1000 mg/L
        return value * 1000.0
    return value

def salinity_to_mol_L(salinity_mg_L):
    """Convert salinity from mg/L NaCl to mol/L."""
    molar_mass_NaCl = 58.44  # g/mol
    return (salinity_mg_L / 1000.0) / molar_mass_NaCl

def solubility_mol_L_to_g_m3(sol_mol_L, gas):
    """
    Convert solubility from mol/L to g/m³ for the given gas.
    1 mol/L = 1000 mol/m³, then × molar_mass (g/mol) → g/m³.
    """
    mm = gas_data[gas]['molar_mass']  # g/mol
    return sol_mol_L * 1000.0 * mm

# Gas solubility using van 't Hoff and salinity correction
def gas_solubility(gas, Temp_C, Pressure_atm, salinity_mg_L):
    """
    Calculate gas solubility in water (mol/L) for a given gas at
    temperature (°C), pressure (atm), and salinity (mg/L NaCl).
    """
    if gas.lower() not in gas_data:
        raise ValueError(f"Unsupported gas type: {gas}")

    props = gas_data[gas.lower()]
    Temp_K = Temp_C + 273.15

    # van 't Hoff relation for Henry's constant
    ln_H = np.log(props['H0']) + (props['delta_H_sol'] / R) * (1.0 / T0 - 1.0 / Temp_K)
    H_T = np.exp(ln_H)

    # salinity correction (Sechenov equation)
    salinity_mol = salinity_to_mol_L(salinity_mg_L)
    correction = 10.0 ** (-props['k_s'] * salinity_mol)

    # solubility (C) = P/H
    C0 = Pressure_atm / H_T
    return C0 * correction


def unesco_density(Temp_C, salinity_mg_L):
    """
    UNESCO (1980) equation of state for seawater density (kg/m^3).
    Assumes salinity in psu approximated by g/kg (~salinity_mg_L/1000).
    """
    # convert salinity to PSU (practical salinity units)
    S = salinity_mg_L / 1000.0
    # pure water density coefficients
    A0 = 999.842594
    A1 = 6.793952e-2
    A2 = -9.095290e-3
    A3 = 1.001685e-4
    A4 = -1.120083e-6
    A5 = 6.536332e-9

    rho_w = (
        A0
        + A1 * Temp_C
        + A2 * Temp_C**2
        + A3 * Temp_C**3
        + A4 * Temp_C**4
        + A5 * Temp_C**5
    )

    # seawater corrections
    B0 = 0.824493
    B1 = -0.0040899
    B2 = 7.6438e-5
    B3 = -8.2467e-7
    B4 = 5.3875e-9

    rho_sw = rho_w + S * (
        B0
        + B1 * Temp_C
        + B2 * Temp_C**2
        + B3 * Temp_C**3
        + B4 * Temp_C**4
    )
    return rho_sw

# Viscosity: dynamic (Pa·s) and kinematic (m²/s)
def calculate_dynamic_viscosity(Temp_C, salinity_mg_L):
    A, B, C = 2.414e-5, 247.8, 140.0
    T_K = Temp_C + 273.15
    mu_w = A * 10.0**(B / (T_K - C))
    S = salinity_mg_L / 1000.0
    return mu_w * (1.0 + 0.001 * S)

def calculate_kinematic_viscosity(mu_dyn, rho):
    return mu_dyn / rho

# Electrical conductivity (S/m)
def calculate_conductivity(Temp_C, salinity_mg_L):
    S = salinity_mg_L / 1000.0
    C_ref, a, b = 4.2914, 0.0207, -0.00004
    return C_ref * (S / 35.0) * (1 + a*(Temp_C - 15.0) + b*(Temp_C - 15.0)**2)

# Viscosity (cP)
def calculate_viscosity(Temp_C, salinity_mg_L):
    """
    Estimate dynamic viscosity of saline water (cP) using empirical relation.
    """
    S = salinity_mg_L / 1000.0
    A, B, C = 2.414e-5, 247.8, 140.0
    # viscosity of pure water (Pa·s) = A * 10^(B/(T- C))
    mu_water = A * 10.0 ** (B / (Temp_C + 273.15 - C))
    # convert Pa·s to cP (1 Pa·s = 1000 cP)
    mu_cP = mu_water * 1000.0
    # salinity correction
    return mu_cP * (1.0 + 0.001 * S)

# Density uncertainty between models
def density_uncertainty(Temp_C,salinity_mg_L) -> tuple[float, float, float]:
    """
    Compute both simple and UNESCO densities and their difference (kg/m^3).
    Returns (rho_simple, rho_unesco, abs_delta).
    """
    rho_s = simple_density(Temp_C, salinity_mg_L)
    rho_u = unesco_density(Temp_C, salinity_mg_L)
    return rho_s, rho_u, abs(rho_s - rho_u)

# Gas flow requirements
def gas_flow_rate_needed(solubility, water_flow_L_min):
    """Compute gas amount needed (mol/min) to achieve saturation.
    solubility in mol/L, flow in L/min."""
    return solubility * water_flow_L_min


def convert_mol_min_to_m3_hr(mol_per_min, Temp_C, Pressure_atm):
    """
    Convert a molar flow rate (mol/min) to volumetric flow (m^3/hr)
    at given temperature (°C) and pressure (atm) using ideal gas law.
    """
    Temp_K = Temp_C + 273.15
    Pressure_Pa = Pressure_atm * 101325.0
    mol_per_sec = mol_per_min / 60.0
    # volume per second (m^3/s) = n*R*T/P
    vol_m3_s = (mol_per_sec * R * Temp_K) / Pressure_Pa
    # convert to m3/hr
    return vol_m3_s * 3600.0

def convert_mol_min_to_normal_m3_hr(mol_per_min, T_norm_K=273.15, P_norm_Pa=101325.0):
    """
    Convert molar flow (mol/min) to Normal m³/hr,
    i.e. at 0 °C (273.15 K) and 1 atm (101 325 Pa).
    """
    mol_per_sec = mol_per_min / 60.0
    vol_m3_s = (mol_per_sec * R * T_norm_K) / P_norm_Pa
    return vol_m3_s * 3600.0

# --- App ---
st.set_page_config('Gas Solubility and Flow Calculator', page_icon="ↂ", layout='wide')
st.title("Gas Solubility and Flow Calculator")

col1, col2 = st.columns(2)

# Temp
with col1:
    temp_input = st.number_input("Temperatur", value=25.0)
    press_input = st.number_input("Pressure", value=1.0)
    t_s = st.number_input("Salinity value", value=35.0)
    gas = st.selectbox("Gas Type", options=list(gas_data.keys()), index=0)


with col2:
    temp_unit = st.selectbox("Unit", ["°C", "°F"], index=0)
    press_unit = st.selectbox("Pressure unit", ['atm','barg','psi'], index=0 )
    unit_s = st.selectbox("Salinity unit", ["mg/L","ppt"], index=1)
    flow_m3_hr = st.number_input("Water Flow Rate (m³/hr)", value=1.0)

# convert to min for functions   
flow_L_min = (flow_m3_hr * 1000.0) / 60.0
salinity_mg_L = salinity_to_mg_L(t_s, unit_s)

# Temp Convert    
Temp_C = fahrenheit_to_celsius(temp_input) if temp_unit == "°F" else temp_input

# Pressure
if press_unit == "barg":
    Pressure_atm = barg_to_atm(press_input)
elif press_unit == "psi":
    Pressure_atm = psi_to_atm(press_input)
else:
    Pressure_atm = press_input


if st.button("Calculate"):
    # solubility & flow
    sol = gas_solubility(gas, Temp_C, Pressure_atm, salinity_mg_L)
    mol_min = gas_flow_rate_needed(sol, flow_L_min)
    actual_m3_hr = convert_mol_min_to_m3_hr(mol_min, Temp_C, Pressure_atm)
    normal_m3_hr = convert_mol_min_to_normal_m3_hr(mol_min)

    rho = unesco_density(Temp_C, salinity_mg_L)
    mu_dyn = calculate_dynamic_viscosity(Temp_C, salinity_mg_L)
    nu = calculate_kinematic_viscosity(mu_dyn, rho)
    cond = calculate_conductivity(Temp_C, salinity_mg_L)

    st.markdown(" ## Results💡")
    st.write(f"**Solubility:** {sol:.6f} mol/L")
    st.write(f"**Gas flow needed:** {mol_min:.4f} mol/min")
    st.write(f"**Gas volume (actual):** {actual_m3_hr:.4f} m³/hr")
    st.write(f"**Gas volume (Normal 0 °C,1 atm):** {normal_m3_hr:.4f} Nm³/hr")
    st.write(f"**Density (UNESCO):** {rho:.2f} kg/m³")
    st.write(f"**Dynamic viscosity:** {mu_dyn:.6e} cP")
    st.write(f"**Kinematic viscosity:** {nu:.6e} m²/s")
    st.write(f"**Electrical conductivity:** {cond:.3f} S/m")

    # Graphs

    # ----- Solubility vs Temp --------

    st.markdown("## Graphs 📊")


    temps = np.linspace(Temp_C - 10.0, Temp_C + 10.0, 20)
    sol_vs_T_mol = [gas_solubility(gas, T, Pressure_atm, salinity_mg_L) for T in temps]
    sol_vs_T_gm3 = [solubility_mol_L_to_g_m3(sol, gas) for sol in sol_vs_T_mol]

    df_Temp = pd.DataFrame({
        "Temperature (°C)": temps,
        "Solubility (g/m³)" : sol_vs_T_gm3
        })
    
    fig_T = px.line(df_Temp, x="Temperature (°C)", y="Solubility (g/m³)",
                    title=f"Solubility vs Temperature @ {Pressure_atm:.2f} atm, {salinity_mg_L:.0f} mg/L for {gas.title()} 🥶",
                    markers=True)
    fig_T.update_traces(line_color = 'orange')
    st.plotly_chart(fig_T)

    # ------- Solubility vs Pressure ------
    pressures = np.linspace(Pressure_atm - 0.9, Pressure_atm + 0.9, 100)
    sol_vs_P_mol  = [gas_solubility(gas, Temp_C, P, salinity_mg_L ) for P in pressures]
    sol_vs_P_gm3 = [solubility_mol_L_to_g_m3(sol, gas) for sol in sol_vs_P_mol]

    df_P = pd.DataFrame({
        "Pressure (atm)"   : pressures,
        "Solubility (g/m³)": sol_vs_P_gm3
        })
    
    fig_P = px.line(df_P, x="Pressure (atm)", y="Solubility (g/m³)", 
                    title=f"Solubility vs Pressure @ {Temp_C:.1f} °C, {salinity_mg_L:.0f} mg/L for {gas.title()} 💎",
                    markers=True,
                    line_shape='linear')
    
    fig_P.update_traces(line_color = 'lightgreen')
    st.plotly_chart(fig_P)

    # ----- Combined Heat Map -----
    data = []

    for T in temps:
        for P in pressures:
            data.append({
            "Temperature (°C)": T,
            "Pressure (atm)"  : P,
            "Solubility"      : solubility_mol_L_to_g_m3(sol_mol_L=(gas_solubility(gas, T, P, salinity_mg_L)), gas=gas)
            })
    
    df_combined = pd.DataFrame(data)

    pivot_combo = df_combined.pivot(index="Temperature (°C)", columns="Pressure (atm)", values="Solubility")

    fig_combo = px.imshow(pivot_combo, aspect='auto', origin='lower',
                          labels={
                                    "x": "Pressure (atm)",
                                    "y": "Temperature (°C)",
                                    "color": "Solubility (g/m³)"
                          },
                          title=f"Solubility Heatmap @ {salinity_mg_L:.0f} mg/L, Gas = {gas.title()} 🥵",
                          color_continuous_scale='Greys')
    
    st.plotly_chart(fig_combo)
