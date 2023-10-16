"""
    Module to calculate the CU Numerical Flux vectors using the 2D CU Scheme
"""

from CentralUpwindScheme.piecewise_linear_construction import *
from CentralUpwindScheme.primitive_variables import *
from CentralUpwindScheme.local_speeds_of_propagation import *
from CentralUpwindScheme.anti_diffusion_term import *

def get_cu_flux(data: dict) -> tuple:
    """
        Function to get the CU Numerical Flux for the CU Scheme
    """

    #! STEP-1 Calculate the Piecewise Linear Reconstructions
    rho_plr = get_plr(data["rho"], data["bc"], "Density")
    mx_plr = get_plr(data["mx"], data["bc"], "Density")
    my_plr = get_plr(data["my"], data["bc"], "Density")
    E_plr = get_plr(data["E"], data["bc"], "Density")

    # Calculating flux vectors for later usage
    g1N, g2N, g3N, g4N = get_gflux(rho_plr[0], mx_plr[0], my_plr[0], E_plr[0])
    g1S, g2S, g3S, g4S = get_gflux(rho_plr[1], mx_plr[1], my_plr[1], E_plr[1])
    f1E, f2E, f3E, f4E = get_fflux(rho_plr[2], mx_plr[2], my_plr[2], E_plr[2])
    f1W, f2W, f3W, f4W = get_fflux(rho_plr[3], mx_plr[3], my_plr[3], E_plr[3])

    #! STEP-2 Calucate the Local Speeds of Propagation
    ap, am = get_lsp_x(rho_plr, mx_plr, my_plr, E_plr)
    bp, bm = get_lsp_y(rho_plr, mx_plr, my_plr, E_plr)

    #! STEP-3 Update the time step for the next iteration using CFL conditions
    update_dt(data, ap, am, bp, bm)

    #! STEP-4 Calculate the Anti Diffusion term
    adtx1, adtx2, adtx3, adtx4 = get_adt_x(rho_plr, mx_plr, my_plr, E_plr, ap, am)
    adty1, adty2, adty3, adty4 = get_adt_y(rho_plr, mx_plr, my_plr, E_plr, bp, bm)

    #! STEP-5 Calculate the CU Flux
    F1 = np.zeros(rho_plr[0].shape)
    F2 = np.zeros(rho_plr[0].shape)
    F3 = np.zeros(rho_plr[0].shape)
    F4 = np.zeros(rho_plr[0].shape)

    G1 = np.zeros(rho_plr[0].shape)
    G2 = np.zeros(rho_plr[0].shape)
    G3 = np.zeros(rho_plr[0].shape)
    G4 = np.zeros(rho_plr[0].shape)

    row, col = rho_plr[0].shape

    for i in range(row-1):
        for j in range(col):

            dist_a = ap[i][j] - am[i][j]
            prod_a = ap[i][j]*am[i][j]

            if dist_a > CTS.EPSILON:
                F1[i][j] = ( ap[i][j]*f1E[i][j] - am[i][j]*f1W[i+1][j] + prod_a*( rho_plr[3][i+1][j] - rho_plr[2][i][j] - adtx1[i][j] ) ) / dist_a
                F2[i][j] = ( ap[i][j]*f2E[i][j] - am[i][j]*f2W[i+1][j] + prod_a*( mx_plr[3][i+1][j] - mx_plr[2][i][j] - adtx2[i][j] ) ) / dist_a
                F3[i][j] = ( ap[i][j]*f3E[i][j] - am[i][j]*f3W[i+1][j] + prod_a*( my_plr[3][i+1][j] - my_plr[2][i][j] - adtx3[i][j] ) ) / dist_a
                F4[i][j] = ( ap[i][j]*f4E[i][j] - am[i][j]*f4W[i+1][j] + prod_a*( E_plr[3][i+1][j] - E_plr[2][i][j] - adtx4[i][j] ) ) / dist_a
            else:
                F1[i][j] = 0.5*( f1E[i][j] + f1W[i][j] )
                F2[i][j] = 0.5*( f2E[i][j] + f2W[i][j] )
                F3[i][j] = 0.5*( f3E[i][j] + f3W[i][j] )
                F4[i][j] = 0.5*( f4E[i][j] + f4W[i][j] )
    
    for i in range(row):
        for j in range(col-1):
            dist_b = bp[i][j] - bm[i][j]
            prod_b = bp[i][j]*bm[i][j]

            if dist_a > CTS.EPSILON:
                G1[i][j] = ( bp[i][j]*g1N[i][j] - bm[i][j]*g1S[i][j+1] + prod_b*( rho_plr[1][i][j+1] - rho_plr[0][i][j] - adty1[i][j] ) ) / dist_b
                G2[i][j] = ( bp[i][j]*g2N[i][j] - bm[i][j]*g2S[i][j+1] + prod_b*( mx_plr[1][i][j+1] - mx_plr[0][i][j] - adty2[i][j] ) ) / dist_b
                G3[i][j] = ( bp[i][j]*g3N[i][j] - bm[i][j]*g3S[i][j+1] + prod_b*( my_plr[1][i][j+1] - my_plr[0][i][j] - adty3[i][j] ) ) / dist_b
                G4[i][j] = ( bp[i][j]*g4N[i][j] - bm[i][j]*g4S[i][j+1] + prod_b*( E_plr[1][i][j+1] - E_plr[0][i][j] - adty4[i][j] ) ) / dist_b
            else:
                G1[i][j] = 0.5*( g1N[i][j] + g1S[i][j] )
                G2[i][j] = 0.5*( g2N[i][j] + g2S[i][j] )
                G3[i][j] = 0.5*( g3N[i][j] + g3S[i][j] )
                G4[i][j] = 0.5*( g4N[i][j] + g4S[i][j] )

    return (F1,F2,F3,F4,G1,G2,G3,G4)


def update_dt(data: dict, ap, am, bp, bm) -> None:
    """
        Function to update the time step for the next iteration using CFL conditions
    """
    amax = 0
    bmax = 0

    for i in range(data["Nx"]):
        for j in range(1, data["Ny"]+1):
            amax = max(amax, ap[i][j], -1*am[i][j])

    for i in range(1, data["Nx"]+1):
        for j in range(data["Ny"]):
            bmax = max(bmax, bp[i][j], -1*bm[i][j])

    data["dt"] = CTS.CFL * min( data["dx"]/amax, data["dy"]/bmax )

    if data["t"] + data["dt"] > data["time"][1]:
        data["dt"] = data["time"][1] - data["t"]

