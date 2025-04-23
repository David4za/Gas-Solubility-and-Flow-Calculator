import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import AbstractState, PT_INPUTS

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

map_names = {
    'methane': 'Methane',
    'co2':     'CO2',
    'nitrogen':'Nitrogen',
    'air':     'Air'
}

def real_gas_props(gas, Temp_C, Pressure_atm):
    """Get real‐gas Z, φ, ρ_mass, μ at T,P from CoolProp using AbstractState."""
    # Map your gas name to the CoolProp fluid string
    fluid = map_names[gas.lower()]    # e.g. "Methane"
    T_K = Temp_C + 273.15
    P_Pa = Pressure_atm * 101325.0

    # Build an AbstractState for the HEOS backend
    AS = AbstractState("HEOS", fluid)
    AS.update(PT_INPUTS, P_Pa, T_K)

    # Now pull out what you need:
    Z        = AS.compressibility_factor()       # dimensionless Z
    phi      = AS.fugacity_coefficient(0)        # fugacity coeff of component 0
    rho_mass = AS.rhomass()                      # kg/m³
    mu       = AS.viscosity()                    # Pa·s

    return Z, phi, rho_mass, mu

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

def gas_solubility(gas, Temp_C, Pressure_atm, salinity_mg_L):
    """Real‐gas Henry’s law + salinity correction."""
    # 1) get real‐gas correction
    Z, phi, _, _ = real_gas_props(gas, Temp_C, Pressure_atm)

    props = gas_data[gas.lower()]
    Temp_K = Temp_C + 273.15

    # 2) van ’t Hoff for H(T)
    ln_H = np.log(props['H0']) + (props['delta_H_sol']/R)*(1/T0 - 1/Temp_K)
    H_T  = np.exp(ln_H)

    # 3) salinity correction
    sal_mol = salinity_to_mol_L(salinity_mg_L)
    corr    = 10.0**(-props['k_s'] * sal_mol)

    # 4) real‐gas solubility: choose φ or Z correction
    C0 = (phi * Pressure_atm) / H_T        # fugacity‐based
    # C0 = (Z   * Pressure_atm) / H_T      # or compressibility‐based

    return C0 * corr

# Gas solubility using van 't Hoff and salinity correction
# def gas_solubility(gas, Temp_C, Pressure_atm, salinity_mg_L):
#     """
#     Calculate gas solubility in water (mol/L) for a given gas at
#     temperature (°C), pressure (atm), and salinity (mg/L NaCl).
#     """
#     if gas.lower() not in gas_data:
#         raise ValueError(f"Unsupported gas type: {gas}")

#     props = gas_data[gas.lower()]
#     Temp_K = Temp_C + 273.15

#     # van 't Hoff relation for Henry's constant
#     ln_H = np.log(props['H0']) + (props['delta_H_sol'] / R) * (1.0 / T0 - 1.0 / Temp_K)
#     H_T = np.exp(ln_H)

#     # salinity correction (Sechenov equation)
#     salinity_mol = salinity_to_mol_L(salinity_mg_L)
#     correction = 10.0 ** (-props['k_s'] * salinity_mol)

#     # solubility (C) = P/H
#     C0 = Pressure_atm / H_T
#     return C0 * correction


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
def calculate_dynamic_viscosity_cP(Temp_C, salinity_mg_L):
    """
    Dynamic viscosity of saline water, returned in centipoise (cP).
    """
    A, B, C = 2.414e-5, 247.8, 140.0
    T_K = Temp_C + 273.15
    mu_water_pa_s = A * 10.0**(B / (T_K - C))
    S = salinity_mg_L / 1000.0
    # include salinity correction (still Pa·s)
    mu_pa_s = mu_water_pa_s * (1.0 + 0.001 * S)
    # convert Pa·s → cP (1 Pa·s = 1 000 cP)
    mu_cP = mu_pa_s * 1000.0
    return mu_cP

# def simple_density(Temp_C, salinity_mg_L):
#      """Simple linear model for water density (kg/m^3)."""
#      S = salinity_mg_L / 1000.0  # convert mg/L -> g/L -> kg/kg approx
#      return 1000.0 + (0.668 * S) - (0.001 * S * Temp_C)

# Density uncertainty between models
def density_uncertainty(Temp_C,salinity_mg_L):
    """
    Compute both simple and UNESCO densities and their difference (kg/m^3).
    Returns (rho_simple, rho_unesco, abs_delta).
    """
    # rho_s = simple_density(Temp_C, salinity_mg_L)
    rho_u = unesco_density(Temp_C, salinity_mg_L)
    return rho_s, rho_u, abs(rho_s - rho_u)

# Gas flow requirements
def gas_flow_rate_needed(solubility, water_flow_L_min):
    """Compute gas amount needed (mol/min) to achieve saturation.
    solubility in mol/L, flow in L/min."""
    return solubility * water_flow_L_min

def convert_mol_min_to_m3_hr(gas, mol_per_min, Temp_C, Pressure_atm):
    """Convert flow [mol/min] → volume [m³/hr] using real‐gas Z from CoolProp."""
    # get real‐gas props (we only care about Z here)
    Z, _, _, _ = real_gas_props(gas, Temp_C, Pressure_atm)

    # convert units
    T_K    = Temp_C + 273.15
    P_Pa   = Pressure_atm * 101325.0
    mol_s  = mol_per_min / 60.0

    # V̇ = ṅ·R·T / (P·Z)
    vol_s  = (mol_s * R * T_K) / (P_Pa * Z)

    # back to m³/hr
    return vol_s * 3600.0

# def convert_mol_min_to_m3_hr(mol_per_min, Temp_C, Pressure_atm):
#     """
#     Convert a molar flow rate (mol/min) to volumetric flow (m^3/hr)
#     at given temperature (°C) and pressure (atm) using ideal gas law.
#     """
#     Temp_K = Temp_C + 273.15
#     Pressure_Pa = Pressure_atm * 101325.0
#     mol_per_sec = mol_per_min / 60.0
#     # volume per second (m^3/s) = n*R*T/P
#     vol_m3_s = (mol_per_sec * R * Temp_K) / Pressure_Pa
#     # convert to m3/hr
#     return vol_m3_s * 3600.0

def convert_mol_min_to_normal_m3_hr(mol_per_min, T_norm_K=273.15, P_norm_Pa=101325.0):
    """
    Convert molar flow (mol/min) to Normal m³/hr,
    i.e. at 0 °C (273.15 K) and 1 atm (101 325 Pa).
    """
    mol_per_sec = mol_per_min / 60.0
    vol_m3_s = (mol_per_sec * R * T_norm_K) / P_norm_Pa
    return vol_m3_s * 3600.0

# ------------------------------------ Streamlit Application ------------------------------------
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
    # Real-gas props for gas (optional display)
    Z, phi, rho_gas, mu_gas = real_gas_props(gas, Temp_C, Pressure_atm)

    # Solubility & flow
    sol         = gas_solubility(gas, Temp_C, Pressure_atm, salinity_mg_L)
    mol_min     = gas_flow_rate_needed(sol, flow_L_min)
    actual_m3_hr= convert_mol_min_to_m3_hr(gas, mol_min, Temp_C, Pressure_atm)  # ← pass gas
    normal_m3_hr= convert_mol_min_to_normal_m3_hr(mol_min)

    # Water properties (you can swap to PropsSI if desired)
    rho_sw      = unesco_density(Temp_C, salinity_mg_L)
    mu_dyn      = calculate_dynamic_viscosity_cP(Temp_C, salinity_mg_L)
    nu          = calculate_kinematic_viscosity(mu_dyn, rho_sw)
    cond        = calculate_conductivity(Temp_C, salinity_mg_L)

    st.markdown(" # Results💡")

    # Two-column layout
    col1, col2 = st.columns(2)

    with col1:
        st.metric("Solubility (mol/L)",          f"{sol:.6f}")
        st.metric("Gas Flow Needed (mol/min)",   f"{mol_min:.6f}")
        st.metric("Actual Gas Volume (m³/hr)",   f"{actual_m3_hr:.6f}")
        st.metric("Normal Gas Volume (Nm³/hr)",  f"{normal_m3_hr:.6f}")

    with col2:
        st.metric("Density (kg/m³)",             f"{rho_sw:62f}")
        st.metric("Dynamic Viscosity (cP)",      f"{mu_dyn:.6f}")
        st.metric("Kinematic Viscosity (m²/s)",  f"{nu:.6f}")
        st.metric("Conductivity (S/m)",          f"{cond:.6f}")

    # Graphs

    # ----- Solubility vs Temp --------

    st.markdown("# Graphs 📊")


    temps = np.linspace(0, Temp_C + 15.0, 20)
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

    # ------- Solubility vs Pressure -------
    pressures = np.linspace(0.1, Pressure_atm + 3, 20)
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

    # ------- Combined Heat Map -------
    data = []

    for T in temps:
        for P in pressures:
            data.append({
            "Temperature (°C)": T,
            "Pressure (atm)"  : P,
            "Solubility"      : solubility_mol_L_to_g_m3(sol_mol_L=(gas_solubility(gas, T, P, salinity_mg_L)), gas=gas)
            })
    
    df_combined = pd.DataFrame(data)

    pivot_combo = df_combined.pivot(index="Pressure (atm)", columns="Temperature (°C)", values="Solubility")

    fig_combo = px.imshow(pivot_combo, aspect='auto', origin='lower',
                          labels={
                                    "x": "Temperature (°C)",
                                    "y": "Pressure (atm)",
                                    "color": "Solubility (g/m³)"
                          },
                          title=f"Solubility Heatmap @ {salinity_mg_L:.0f} mg/L, Gas = {gas.title()} 🥵",
                          color_continuous_scale='Greys')
    
    st.plotly_chart(fig_combo)
