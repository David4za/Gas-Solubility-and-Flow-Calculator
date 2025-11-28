import numpy as np
from CoolProp.CoolProp import AbstractState, PT_INPUTS
from .constants import MAP_NAMES, GAS_DATA, R, T0
from .salinity import salinity_to_mol_L

def real_gas_props(gas, Temp_C, Pressure_atm):
    """Get real‐gas Z, φ, ρ_mass, μ at T,P from CoolProp using AbstractState."""
    # Map your gas name to the CoolProp fluid string
    fluid = MAP_NAMES[gas.lower()]    # e.g. "Methane"
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

def solubility_mol_L_to_g_m3(sol_mol_L, gas):
    """
    Convert solubility from mol/L to g/m³ for the given gas.
    1 mol/L = 1000 mol/m³, then × molar_mass (g/mol) → g/m³.
    """
    mm = GAS_DATA[gas]['molar_mass']  # g/mol
    return sol_mol_L * 1000.0 * mm