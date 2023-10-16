"""
    Module to calculate the Flux vectors from the Eulers System of Gas Dynamics
"""

import numpy as np

from CentralUpwindScheme.primitive_variables import *

def get_fflux(rho, mx, my, E) -> tuple:
    """
        Function to calculate F(U) given U
    """

    u, v, p = get_primitive_variables(rho, mx, my, E)

    f1 = rho*u
    f2 = rho*(u**2) + p
    f3 = rho*u*v
    f4 = u*(E + p)

    return (f1,f2,f3,f4)

def get_gflux(rho, mx, my, E) -> tuple:
    """
        Function to calculate G(U) given U
    """

    u, v, p = get_primitive_variables(rho, mx, my, E)

    g1 = rho*v
    g3 = rho*u*v
    g2 = rho*(v**2) + p
    g4 = v*(E + p)

    return (g1,g2,g3,g4)
