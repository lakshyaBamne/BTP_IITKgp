"""
    Module to store the global constants used in the program
"""

GAMMA = 1.4 # specific heat
THETA = 1.3 # tuning parameter
EPSILON = 1.0E-12 # epsilon used in Flux calculation
CFL = 0.475 # CFL Number for calculating time steps
SMOOTHP = 0.00625 # Smooting parameter used in Kelvin Helmholtz Instability

def show_constants():
    """
        Function to show the constants being used in the problem
    """
    print(f"GAMMA = {GAMMA}")
    print(f"TUNING PARAMETER = {THETA}")
    print(f"EPSILON = {EPSILON}")
    print(f"CFL NUMBER = {CFL}")
    print(f"SMOOTHING PARAMETER = {SMOOTHP}")

