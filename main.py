import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px

from gas_calcs.constants import GAS_DATA
from gas_calcs.gas_properties import real_gas_props, solubility_mol_L_to_g_m3
from gas_calcs.flow import gas_flow_rate_needed, convert_mol_min_to_m3_hr, convert_mol_min_to_normal_m3_hr
from gas_calcs.units import barg_to_atm, psi_to_atm, fahrenheit_to_celsius
from gas_calcs.solubility import gas_solubility
from gas_calcs.salinity import salinity_to_mg_L
from gas_calcs.water_properties import unesco_density, calculate_dynamic_viscosity_cP, calculate_kinematic_viscosity, calculate_conductivity


# ------------------------------------ Streamlit Application ------------------------------------
st.set_page_config('Gas Solubility and Flow Calculator', page_icon="ↂ", layout='wide')
st.title("Gas Solubility and Flow Calculator")

with st.expander("Input Parameters"):
    st.info("Default values")
    col1, col2 = st.columns(2)

    # value
    with col1:
        temp_input = st.number_input("Temperatur", value=25.0)
        press_input = st.number_input("Pressure", value=1.0)
        t_s = st.number_input("Salinity value", value=35.0)
        gas = st.selectbox("Gas Type", options=list(GAS_DATA.keys()), index=0)

    # unit of measure 
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
