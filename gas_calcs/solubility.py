import imp
import numpy as np
from .constants import R, T0, GAS_DATA
from .gas_properties import real_gas_props
from .salinity import salinity_to_mol_L

def gas_solubility(gas, Temp_C, Pressure_atm, salinity_mg_L):
    """Real‐gas Henry’s law + salinity correction."""
    # 1) get real‐gas correction
    Z, phi, _, _ = real_gas_props(gas, Temp_C, Pressure_atm)

    props = GAS_DATA[gas.lower()]
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