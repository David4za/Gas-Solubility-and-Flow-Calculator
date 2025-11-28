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

# Density uncertainty between models
def density_uncertainty(Temp_C,salinity_mg_L):
    """
    Compute both simple and UNESCO densities and their difference (kg/m^3).
    Returns (rho_simple, rho_unesco, abs_delta).
    """
    # rho_s = simple_density(Temp_C, salinity_mg_L)
    rho_u = unesco_density(Temp_C, salinity_mg_L)
    return rho_s, rho_u, abs(rho_s - rho_u)