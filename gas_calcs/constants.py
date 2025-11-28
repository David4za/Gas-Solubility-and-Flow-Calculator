R = 8.314  # universal gas constant, J/(mol·K)
T0 = 298.15  # reference temperature (25 °C) in Kelvin

# Gas properties: H0 (Henry's constant at T0), delta_H_sol (enthalpy of solution),
# k_s (salinity correction factor), molar_mass (g/mol)
GAS_DATA = {
    'methane': {'H0': 1400.0, 'delta_H_sol': 8300.0, 'k_s': 0.12, 'molar_mass': 16.04},
    'co2':     {'H0': 29.4,   'delta_H_sol': 19000.0, 'k_s': 0.10, 'molar_mass': 44.01},
    'nitrogen':{'H0': 1600.0, 'delta_H_sol': 5800.0,  'k_s': 0.13, 'molar_mass': 28.02},
    'air':     {'H0': 770.0,  'delta_H_sol': 10000.0, 'k_s': 0.12, 'molar_mass': 28.97}
}

MAP_NAMES = {
    'methane': 'Methane',
    'co2':     'CO2',
    'nitrogen':'Nitrogen',
    'air':     'Air'
}