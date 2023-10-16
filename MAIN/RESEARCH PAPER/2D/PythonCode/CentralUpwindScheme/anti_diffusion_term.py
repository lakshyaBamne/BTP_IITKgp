"""
    Module to calculate the Anti diffusion term for the 2D CU Scheme
"""
import numpy as np

from CentralUpwindScheme.primitive_variables import *
from CentralUpwindScheme.flux_vectors import *
from CentralUpwindScheme.utility import *

def get_star_x(rho_plr, mx_plr, my_plr, E_plr, lspp, lspm):
    """
        Function to get the intermediate term for calculating Anti Diffusion term
    """
    row, col = rho_plr[0].shape

    f1E, f2E, f3E, f4E = get_fflux(rho_plr[2], mx_plr[2], my_plr[2], E_plr[2])
    f1W, f2W, f3W, f4W = get_fflux(rho_plr[3], mx_plr[3], my_plr[3], E_plr[3])

    starx1 = np.zeros(rho_plr[0].shape)
    starx2 = np.zeros(rho_plr[0].shape)
    starx3 = np.zeros(rho_plr[0].shape)
    starx4 = np.zeros(rho_plr[0].shape)

    for i in range(row-1):
        for j in range(col):
            starx1[i][j] = (lspp[i][j]*rho_plr[3][i+1][j] - lspm[i][j]*rho_plr[2][i][j] - f1W[i+1][j] + f1E[i][j])/(lspp[i][j]-lspm[i][j])
            starx2[i][j] = (lspp[i][j]*mx_plr[3][i+1][j] - lspm[i][j]*mx_plr[2][i][j] - f2W[i+1][j] + f2E[i][j])/(lspp[i][j]-lspm[i][j])
            starx3[i][j] = (lspp[i][j]*my_plr[3][i+1][j] - lspm[i][j]*my_plr[2][i][j] - f3W[i+1][j] + f3E[i][j])/(lspp[i][j]-lspm[i][j])
            starx4[i][j] = (lspp[i][j]*E_plr[3][i+1][j] - lspm[i][j]*E_plr[2][i][j] - f4W[i+1][j] + f4E[i][j])/(lspp[i][j]-lspm[i][j])

    return (starx1, starx2, starx3, starx4)


def get_star_y(rho_plr, mx_plr, my_plr, E_plr, lspp, lspm):
    """
        Function to get the intermediate term for calculating Anti Diffusion term
    """
    row, col = rho_plr[0].shape

    g1N, g2N, g3N, g4N = get_gflux(rho_plr[0], mx_plr[0], my_plr[0], E_plr[0])
    g1S, g2S, g3S, g4S = get_gflux(rho_plr[1], mx_plr[1], my_plr[1], E_plr[1])

    stary1 = np.zeros(rho_plr[0].shape)
    stary2 = np.zeros(rho_plr[0].shape)
    stary3 = np.zeros(rho_plr[0].shape)
    stary4 = np.zeros(rho_plr[0].shape)

    for i in range(row):
        for j in range(col-1):
            stary1[i][j] = (lspp[i][j]*rho_plr[1][i][j+1] - lspm[i][j]*rho_plr[0][i][j] - g1S[i][j+1] + g1N[i][j])/(lspp[i][j]-lspm[i][j])
            stary2[i][j] = (lspp[i][j]*mx_plr[1][i][j+1] - lspm[i][j]*mx_plr[0][i][j] - g2S[i][j+1] + g2N[i][j])/(lspp[i][j]-lspm[i][j])
            stary3[i][j] = (lspp[i][j]*my_plr[1][i][j+1] - lspm[i][j]*my_plr[0][i][j] - g3S[i][j+1] + g3N[i][j])/(lspp[i][j]-lspm[i][j])
            stary4[i][j] = (lspp[i][j]*E_plr[1][i][j+1] - lspm[i][j]*E_plr[0][i][j] - g4S[i][j+1] + g4N[i][j])/(lspp[i][j]-lspm[i][j])

    return (stary1, stary2, stary3, stary4)

def get_adt_x(rho_plr, mx_plr, my_plr, E_plr, lspp, lspm):
    """
        Function to calculate the Anti diffusion term
    """
    row, col = rho_plr[0].shape

    starx1, starx2, starx3, starx4 = get_star_x(rho_plr, mx_plr, my_plr, E_plr, lspp, lspm)

    adt1 = np.zeros((row,col))
    adt2 = np.zeros((row,col))
    adt3 = np.zeros((row,col))
    adt4 = np.zeros((row,col))

    for i in range(row-1):
        for j in range(col):
            adt1[i][j] = minmod(
                rho_plr[7][i+1][j] - starx1[i][j],
                starx1[i][j] - rho_plr[6][i][j],
                rho_plr[5][i+1][j] - starx1[i][j],
                starx1[i][j] - rho_plr[4][i][j]
            )

            adt2[i][j] = minmod(
                mx_plr[7][i+1][j] - starx2[i][j],
                starx2[i][j] - mx_plr[6][i][j],
                mx_plr[5][i+1][j] - starx2[i][j],
                starx2[i][j] - mx_plr[4][i][j]
            )

            adt3[i][j] = minmod(
                my_plr[7][i+1][j] - starx3[i][j],
                starx3[i][j] - my_plr[6][i][j],
                my_plr[5][i+1][j] - starx3[i][j],
                starx3[i][j] - my_plr[4][i][j]
            )

            adt4[i][j] = minmod(
                E_plr[7][i+1][j] - starx4[i][j],
                starx4[i][j] - E_plr[6][i][j],
                E_plr[5][i+1][j] - starx4[i][j],
                starx4[i][j] - E_plr[4][i][j]
            )
    
    return (adt1, adt2, adt3, adt4)


def get_adt_y(rho_plr, mx_plr, my_plr, E_plr, lspp, lspm):
    """
        Function to calculate the Anti diffusion term
    """
    row, col = rho_plr[0].shape

    stary1, stary2, stary3, stary4 = get_star_y(rho_plr, mx_plr, my_plr, E_plr, lspp, lspm)

    adt1 = np.zeros((row,col))
    adt2 = np.zeros((row,col))
    adt3 = np.zeros((row,col))
    adt4 = np.zeros((row,col))

    for i in range(row):
        for j in range(col-1):
            adt1[i][j] = minmod(
                rho_plr[7][i][j+1] - stary1[i][j],
                stary1[i][j] - rho_plr[5][i][j],
                rho_plr[6][i][j+1] - stary1[i][j],
                stary1[i][j] - rho_plr[4][i][j]
            )

            adt2[i][j] = minmod(
                mx_plr[7][i][j+1] - stary2[i][j],
                stary2[i][j] - mx_plr[5][i][j],
                mx_plr[6][i][j+1] - stary2[i][j],
                stary2[i][j] - mx_plr[4][i][j]
            )

            adt3[i][j] = minmod(
                my_plr[7][i][j+1] - stary3[i][j],
                stary3[i][j] - my_plr[5][i][j],
                my_plr[6][i][j+1] - stary3[i][j],
                stary3[i][j] - my_plr[4][i][j]
            )

            adt4[i][j] = minmod(
                E_plr[7][i][j+1] - stary4[i][j],
                stary4[i][j] - E_plr[5][i][j],
                E_plr[6][i][j+1] - stary4[i][j],
                stary4[i][j] - E_plr[4][i][j]
            )
    
    return (adt1, adt2, adt3, adt4)

