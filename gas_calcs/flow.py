from .constants import R
from .gas_properties import real_gas_props

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

def convert_mol_min_to_normal_m3_hr(mol_per_min, T_norm_K=273.15, P_norm_Pa=101325.0):
    """
    Convert molar flow (mol/min) to Normal m³/hr,
    i.e. at 0 °C (273.15 K) and 1 atm (101 325 Pa).
    """
    mol_per_sec = mol_per_min / 60.0
    vol_m3_s = (mol_per_sec * R * T_norm_K) / P_norm_Pa
    return vol_m3_s * 3600.0