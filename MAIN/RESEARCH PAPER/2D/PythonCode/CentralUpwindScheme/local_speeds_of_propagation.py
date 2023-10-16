"""
    Module to calculate the One Sided Local Speeds of propagation
"""

import numpy as np
from math import sqrt

import CentralUpwindScheme.constants as CTS
from CentralUpwindScheme.primitive_variables import get_primitive_variables


def get_lsp_x(rho_plr, mx_plr, my_plr, E_plr) -> tuple:
    """
        Function to calculate the Local Speeds of Propagations
    """

    uN, vN, pN = get_primitive_variables(rho_plr[0], mx_plr[0], my_plr[0], E_plr[0])
    uS, vS, pS = get_primitive_variables(rho_plr[1], mx_plr[1], my_plr[1], E_plr[1])
    uE, vE, pE = get_primitive_variables(rho_plr[2], mx_plr[2], my_plr[2], E_plr[2])
    uW, vW, pW = get_primitive_variables(rho_plr[3], mx_plr[3], my_plr[3], E_plr[3])

    row, col = rho_plr[0].shape
    
    lspp = np.zeros(rho_plr[0].shape)
    lspm = np.zeros(rho_plr[0].shape)

    for i in range(row-1):
        for j in range(col):
            lspp[i][j] = max(
                0,
                uW[i+1][j] + sqrt( (CTS.GAMMA*pW[i+1][j])/rho_plr[3][i+1][j] ),
                uE[i][j] + sqrt( (CTS.GAMMA*pE[i][j])/rho_plr[2][i][j] )
            )

            lspm[i][j] = min(
                0,
                uW[i+1][j] - sqrt( (CTS.GAMMA*pW[i+1][j])/rho_plr[3][i+1][j] ),
                uE[i][j] - sqrt( (CTS.GAMMA*pE[i][j])/rho_plr[2][i][j] )
            )
    
    return (lspp, lspm)

def get_lsp_y(rho_plr, mx_plr, my_plr, E_plr) -> tuple:
    """
        Function to calculate the Local Speeds of Propagations
    """

    uN, vN, pN = get_primitive_variables(rho_plr[0], mx_plr[0], my_plr[0], E_plr[0])
    uS, vS, pS = get_primitive_variables(rho_plr[1], mx_plr[1], my_plr[1], E_plr[1])
    uE, vE, pE = get_primitive_variables(rho_plr[2], mx_plr[2], my_plr[2], E_plr[2])
    uW, vW, pW = get_primitive_variables(rho_plr[3], mx_plr[3], my_plr[3], E_plr[3])

    row, col = rho_plr[0].shape
    
    lspp = np.zeros(rho_plr[0].shape)
    lspm = np.zeros(rho_plr[0].shape)

    for i in range(row):
        for j in range(col-1):
            lspp[i][j] = max(
                0,
                vW[i][j+1] + sqrt( (CTS.GAMMA*pS[i][j+1])/rho_plr[1][i][j+1] ),
                vE[i][j] + sqrt( (CTS.GAMMA*pN[i][j])/rho_plr[0][i][j] )
            )

            lspm[i][j] = min(
                0,
                vW[i][j+1] - sqrt( (CTS.GAMMA*pS[i][j+1])/rho_plr[1][i][j+1] ),
                vE[i][j] - sqrt( (CTS.GAMMA*pN[i][j])/rho_plr[0][i][j] )
            )
    
    return (lspp, lspm)

