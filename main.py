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

# Salinity conversion
def salinity_to_mol_L(salinity_mg_L):
    """Convert salinity from mg/L (assumed NaCl) to mol/L."""
    molar_mass_NaCl = 58.44  # g/mol
    # mg/L -> g/L -> mol/L
    return (salinity_mg_L / 1000.0) / molar_mass_NaCl


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

# Density calculations
def simple_density(Temp_C, salinity_mg_L):
    """Simple linear model for water density (kg/m^3)."""
    S = salinity_mg_L / 1000.0  # convert mg/L -> g/L -> kg/kg approx
    return 1000.0 + (0.668 * S) - (0.001 * S * Temp_C)

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

# --- App ---

st.title("Gas Solubility and Flow Calculator")

col1, col2 = st.columns(2)

# Temp
with col1:
    temp_input = st.number_input("Temperatur", value=25.0)
    press_input = st.number_input("Pressure", value=1.0)
    salinity = st.number_input("Salinity (mg/L)", value=35000.0)

with col2:
    temp_unit = st.selectbox("Unit", ["°C", "°F"], index=0)
    press_unit = st.selectbox("Pressure unit", ['atm','barg','psi'], index=0 )
    flow_rate = st.number_input("Water Flow Rate (L/min)", value=100.0)

# Gas
gas = st.selectbox("Gas Type", options=list(gas_data.keys()), index=0)



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
    sol = gas_solubility(gas, Temp_C, Pressure_atm, salinity)
    mol_per_min = gas_flow_rate_needed(sol, flow_rate)
    m3_per_hr = convert_mol_min_to_m3_hr(mol_per_min, Temp_C, Pressure_atm)
    rho_s = simple_density(Temp_C, salinity)
    rho_u = unesco_density(Temp_C, salinity)
    rho_diff = abs(rho_s - rho_u)
    visc = calculate_viscosity(Temp_C, salinity)

    st.markdown(" ## Results ∰")
    st.write(f"**Solubility:** {sol:.6f} mol/L")
    st.write(f"**Gas flow needed:** {mol_per_min:.4f} mol/min")
    st.write(f"**Gas volume:** {m3_per_hr:.4f} m³/hr")
    st.write(f"**Simple density:** {rho_s:.2f} kg/m³")
    st.write(f"**UNESCO density:** {rho_u:.2f} kg/m³")
    st.write(f"**Density difference:** ±{rho_diff:.2f} kg/m³")
    st.write(f"**Viscosity:** {visc:.6f} cP")

    hours = list(range(0,25))
    buildup = [m3_per_hr * h for h in hours]

    df = pd.DataFrame({
        "Hours": hours,
        "Cumulative Gas Volume (m³)": buildup
    })

    fig = px.line(df, x="Hours", y="Cumulative Gas Volume (m³)", title="Gas Buildup Over Time",
                labels={"Hours": 'Time (Hours)', "Cumulative Gas Volume (m³)": 'Gas Volume (m³)'},
                markers=True)
    st.markdown("### Ploty Chart")
    st.plotly_chart(fig)
