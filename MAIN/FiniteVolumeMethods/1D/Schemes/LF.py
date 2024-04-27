"""
    Module to implement the Lax Friedrichs (LF) Scheme for
    1D Euler Equations of Gas Dynamics
"""
import math

from Initialize.InitRiemannProblem import InitProblem
from Visualize.PlotDensity import plot_density_withref, plot_density

def update_primitive_variables(RP: dict) -> None:
    """
        Function to calculate the primitive variables given the conserved variables
    """

    # we need to update the primitive variables using the conserved variables
    for i in range(1, RP["N"]+1):
        RP["u"][i] = RP["U"][1][i]/RP["U"][0][i]
        RP["p"][i] = (RP["GAMMA"]-1)*(RP["U"][2][i]-0.5*RP["U"][0][i]*RP["u"][i]**2)
        # RP["a"][i] = math.sqrt(abs(RP["GAMMA"]*RP["p"][i]/RP["U"][0][i]))
        RP["a"][i] = math.sqrt(RP["GAMMA"]*RP["p"][i]/RP["U"][0][i])

def calculate_flux_vectors(RP: dict) -> list[list]:
    """
        Function to calculate the flux vectors using the conserved variables
    """

    # Flux vectors
    F1 = [0 for i in range(RP["N"]+2)]
    F2 = [0 for i in range(RP["N"]+2)]
    F3 = [0 for i in range(RP["N"]+2)]

    for i in range(RP["N"]+2):
        F1[i] = RP["U"][1][i]
        F2[i] = RP["U"][0][i]*RP["u"][i]**2 + RP["u"][i]
        F3[i] = RP["u"][i]*(RP["U"][2][i] + RP["p"][i])

    return [F1, F2, F3]

def extend_cells(RP: dict) -> None:
    """
        Function to extend the cells
    """

    if RP["BC"] == "FREE":
        # extend conserved variables
        for i in range(3):
            RP["U"][i][0] = RP["U"][i][1]
            RP["U"][i][-1] = RP["U"][i][-2]

        # extend primitive variables
        RP["u"][0] = RP["u"][1]
        RP["u"][-1] = RP["u"][-2]

        RP["p"][0] = RP["p"][1]
        RP["p"][-1] = RP["p"][-2]

        RP["a"][0] = RP["a"][1]
        RP["a"][-1] = RP["a"][-2]
    elif RP["BC"] == "REFLECTIVE":
        for i in [0,2]:
            RP["U"][i][0] = RP["U"][i][1]
            RP["U"][i][-1] = RP["U"][i][-2]

        # extend momentum
        RP["U"][1][0] = -RP["U"][1][1]
        RP["U"][1][-1] = -RP["U"][1][-2]

        # extend primitive variables
        RP["u"][0] = -RP["u"][1]
        RP["u"][-1] = -RP["u"][-2]

        RP["p"][0] = RP["p"][1]
        RP["p"][-1] = RP["p"][-2]

        RP["a"][0] = RP["a"][1]
        RP["a"][-1] = RP["a"][-2]        
    else:
        print(f"---ERROR--- Please select a correct boundary conditions---")

def update_conserved_variables(RP: dict, FLX: list[list], dt: float, dx: float) -> None:
    """
        Function to update the conserved variables
        using Euler Forward Differences
    """

    for i in range(1, RP["N"]+1):
        RP["U"][0][i] = RP["U"][0][i] - (dt/dx)*(FLX[0][i] - FLX[0][i-1])
        RP["U"][1][i] = RP["U"][1][i] - (dt/dx)*(FLX[1][i] - FLX[1][i-1])
        RP["U"][2][i] = RP["U"][2][i] - (dt/dx)*(FLX[2][i] - FLX[2][i-1])

def fvm1d_lf(RP: dict) -> None:
    """
        Function to implement the 1D Lax Friedrichs Scheme

        -> calculates the time step (dt) using CFL conditions
        -> calculate the Flux vectors (F) for all finite volume cells
        -> calculate the Lax-Friedrichs Fluxes at the intercell interfaces
    """
    # we need to update the primitive variables using the conserved variables
    update_primitive_variables(RP)    

    amax = 0
    for i in range(1, RP["N"]+1):
        amax = max(amax, abs(RP["u"][i])+RP["a"][i])

    dx = float((RP["domx"][1]-RP["domx"][0])/RP["N"])

    #! update the time step
    dt = RP["CFL"]*dx/amax

    # Extend cells for the conserved variables
    extend_cells(RP)

    #! Calculate the Flux vectors using the conserved variables
    F = calculate_flux_vectors(RP)

    # now we can calculate the Lax Friedrichs Fluxes at the interfaces
    LF1 = [0 for i in range(RP["N"]+1)]
    LF2 = [0 for i in range(RP["N"]+1)]
    LF3 = [0 for i in range(RP["N"]+1)]

    for i in range(RP["N"]+1):
        LF1[i] = 0.5*(F[0][i]+F[0][i+1]) - 0.5*(dx/dt)*(RP["U"][0][i+1] - RP["U"][0][i])
        LF2[i] = 0.5*(F[1][i]+F[1][i+1]) - 0.5*(dx/dt)*(RP["U"][1][i+1] - RP["U"][1][i])
        LF3[i] = 0.5*(F[2][i]+F[2][i+1]) - 0.5*(dx/dt)*(RP["U"][2][i+1] - RP["U"][2][i])

    FLX = [LF1, LF2, LF3]

    update_conserved_variables(RP, FLX, dt, dx)

    return dt

def run_lf(PROBLEM: str, ORDER: str) -> None:
    """
        Function to run the scheme
    """

    # Initialize the dictionaries containing the variables
    RP = InitProblem(PROBLEM, "NORMAL", ORDER)
    RPref = InitProblem(PROBLEM, "REF", ORDER)

    # Run the Normal mode simualtion
    print(f"---LOG--- Running Riemann Problem (Normal)---")
    t = 0
    while t<RP["tfinal"]:
        dt = fvm1d_lf(RP)

        t += dt

        # Log message
        print(f"t={t}|dt={dt}")

    # Run the Reference mode simulation
    print(f"---LOG--- Running Riemann Problem (Reference)---")
    t = 0
    while t<RPref["tfinal"]:
        dt = fvm1d_lf(RPref)

        t += dt

        # Log message
        print(f"t={t}|dt={dt}")

    # Visualise the final density (normal and ref)
    plot_density_withref(RP, RPref)
