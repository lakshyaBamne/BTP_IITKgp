"""
    Module to calculate the Piecewise Linear Reconstruction
"""
import numpy as np

import CentralUpwindScheme.constants as CTS
import CentralUpwindScheme.utility as UTL
import CentralUpwindScheme.extend_cells as EXC

def get_slope(matrix, bc: dict, var_name: str) -> tuple:
    """
        Function to calculate the slope for a given variable
    """

    row, col = matrix.shape

    slx = np.zeros(matrix.shape)
    sly = np.zeros(matrix.shape)

    for i in range(1,row-1):
        for j in range(1,col-1):
            slx[i][j] = UTL.minmod(
                CTS.THETA*(matrix[i][j]-matrix[i-1][j]),
                0.5*(matrix[i+1][j]-matrix[i-1][j]),
                CTS.THETA*(matrix[i+1][j]-matrix[i][j])
            )

    for i in range(1,row-1):
        for j in range(1,col-1):
            sly[i][j] = UTL.minmod(
                CTS.THETA*(matrix[i][j]-matrix[i][j-1]),
                0.5*(matrix[i][j+1]-matrix[i][j-1]),
                CTS.THETA*(matrix[i][j+1]-matrix[i][j])
            )

    EXC.extend_matrix(slx, bc, var_name)
    EXC.extend_matrix(sly, bc, var_name)

    return (slx, sly)

def get_plr(matrix, bc: dict, var_name: str) -> tuple:
    """
        Function to get the Piecewise Linear Reconstruction for a variable
    """

    row, col = matrix.shape

    slx, sly = get_slope(matrix, bc, var_name)

    n = np.zeros(matrix.shape)
    s = np.zeros(matrix.shape)
    e = np.zeros(matrix.shape)
    w = np.zeros(matrix.shape)
    ne = np.zeros(matrix.shape)
    nw = np.zeros(matrix.shape)
    se = np.zeros(matrix.shape)
    sw = np.zeros(matrix.shape)

    for i in range(1, row-1):
        for j in range(1, col-1):
            n[i][j] = matrix[i][j] + 0.5*sly[i][j]
            s[i][j] = matrix[i][j] - 0.5*sly[i][j]
            e[i][j] = matrix[i][j] + 0.5*slx[i][j]
            w[i][j] = matrix[i][j] - 0.5*slx[i][j]
            ne[i][j] = matrix[i][j] + 0.5*slx[i][j] + 0.5*sly[i][j]
            nw[i][j] = matrix[i][j] - 0.5*slx[i][j] + 0.5*sly[i][j]
            se[i][j] = matrix[i][j] + 0.5*slx[i][j] - 0.5*sly[i][j]
            sw[i][j] = matrix[i][j] - 0.5*slx[i][j] - 0.5*sly[i][j]

    EXC.extend_matrix(n, bc, "PLR")
    EXC.extend_matrix(s, bc, "PLR")
    EXC.extend_matrix(e, bc, "PLR")
    EXC.extend_matrix(w, bc, "PLR")
    EXC.extend_matrix(ne, bc, "PLR")
    EXC.extend_matrix(nw, bc, "PLR")
    EXC.extend_matrix(se, bc, "PLR")
    EXC.extend_matrix(sw, bc, "PLR")

    return (n,s,e,w,ne,nw,se,sw)