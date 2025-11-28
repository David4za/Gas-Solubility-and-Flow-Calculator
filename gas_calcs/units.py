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