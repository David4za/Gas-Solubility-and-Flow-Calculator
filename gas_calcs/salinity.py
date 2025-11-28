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