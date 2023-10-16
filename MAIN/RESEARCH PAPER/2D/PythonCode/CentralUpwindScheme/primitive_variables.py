"""
    Module to calculate the Primitive variables from the conserved variables
    -> Conserved variables : rho, mx, my, E
    -> Primitive variables : rho, u, v, p
"""

import numpy as np

import CentralUpwindScheme.constants as CTS

def get_primitive_variables(rho, mx, my, E) -> tuple:
    """
        Function to calculate the primitive variables from the given set of conserved variables
    """

    u = mx/rho
    v = my/rho
    p = (CTS.GAMMA-1)*(E - 0.5*rho*(u**2 + v**2))

    return (u,v,p)
