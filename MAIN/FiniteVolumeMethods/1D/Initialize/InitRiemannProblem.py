"""
    Module contains functions to initialize the riemann problem
    given the problem's name
"""

import numpy as np
import math

def InitProblem(PROBLEM: str, MODE: str, ORDER: str) -> dict:
    """
        Function takes the name of the problem as input and returns 
        a dictionary containing the required variables initialized

        PROBLEM -> [MCW, SCW, BLW, TORO-1, TORO-2, TORO-3, TORO-4, LAX, SDW, SEW]
        MODE -> [NORMAL, REF]
        ORDER -> [1, 2]
    """

    RP = {
        "PROBLEM" : PROBLEM,
        "MODE" : MODE,
        "ORDER" : ORDER,
        "GAMMA" : 1.4,
        "THETA" : 1.3,
        "EPS" : 1.0E-12
    }

    if ORDER == "1":
        RP["CFL"] = 0.9
    elif ORDER == "2":
        RP["CFL"] = 0.45
    else:
        print(f"---ERROR--- Please enter correct order of reconstruction---")

    # Initialize the variables which are different for Normal and Reference simulations
    if MODE == "NORMAL":
        RP["N"] = 1000
    elif MODE == "REF":
        RP["N"] = 4000
    else:
        print(f"---ERROR--- Please select a correct mode of operation---")

    # Initialize the variables which are different for different problems 
    if PROBLEM == "MCW":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 2.0
        RP["BC"] = "FREE"
    elif PROBLEM == "SCW":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.012
        RP["BC"] = "FREE"
    elif PROBLEM == "BLW":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.038
        RP["BC"] = "REFLECTIVE"
    elif PROBLEM == "LAX":
        RP["domx"] = (-5.0, 5.0)
        RP["tfinal"] = 1.3
        RP["BC"] = "FREE"
    elif PROBLEM == "SDW":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.012
        RP["BC"] = "FREE"
    elif PROBLEM == "SEW":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.012
        RP["BC"] = "FREE"
    elif PROBLEM == "TORO-1":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.2
        RP["BC"] = "FREE"
    elif PROBLEM == "TORO-2":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.15
        RP["BC"] = "FREE"
    elif PROBLEM == "TORO-3":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.012
        RP["BC"] = "FREE"
    elif PROBLEM == "TORO-4":
        RP["domx"] = (0.0, 1.0)
        RP["tfinal"] = 0.035
        RP["BC"] = "FREE"
    else:
        print(f"---ERROR--- Please select a correct problem for simulation---")

    # Initialize the computational grid
    RP["grid"] = InitComputationalGrid(RP)

    # Initialize the conserved variables for the given problem
    RP["U"] = InitConservedVariables(RP)

    # Initialize the useful primitive variables
    prim_vars = InitPrimitiveVariables(RP)
    
    RP["u"] = prim_vars[0]
    RP["p"] = prim_vars[1]
    RP["a"] = prim_vars[2]

    return RP
    
def InitComputationalGrid(RP: dict) -> list:
    """
        Function to initialize a computational grid
    """

    # first we need to calculate the length of an individual interval
    dx = float(( RP["domx"][1] - RP["domx"][0] ) / RP["N"])

    # grid contains the centers of the finite volume cells
    grid = [RP["domx"][0]+(i+0.5)*dx for i in range(RP["N"])]

    return grid

def InitConservedVariables(RP: dict) -> list[list]:
    """
        Function to initialize the conserved variables
    """

    U = [[0 for i in range(RP["N"]+2)] for j in range(3)]

    if RP["PROBLEM"] == "MCW":

        for i in range(1, RP["N"]+1):
            if RP["grid"][i-1] < 0.3:
                U[0][i] = 1.4
            else:
                U[0][i] = 1.0

            u = 0.1
            p = 1.0

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2

    elif RP["PROBLEM"] == "SCW":

        for i in range(1, RP["N"]+1):

            if RP["grid"][i-1] < 0.8:
                p = 1000
            else:
                p = 0.01

            U[0][i] = 1.0
            u = -19.59745

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2
    
    elif RP["PROBLEM"] == "BLW":

        for i in range(1, RP["N"]+1):

            if RP["grid"][i-1] < 0.1:
                U[0][i] = 1.0
                u = 0
                p = 1000
            elif RP["grid"][i-1]>=0.1 and RP["grid"][i-1]<=0.9:
                U[0][i] = 1.0
                u = 0
                p = 0.01
            else:
                U[0][i] = 1.0
                u = 0
                p = 100

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2
    
    elif RP["PROBLEM"] == "TORO-1":

        for i in range(1, RP["N"]+1):
            if RP["grid"][i-1] < 0.3:
                U[0][i] = 1.0
                u = 0.75
                p = 1.0
            else:
                U[0][i] = 0.125
                u = 0.0
                p = 0.1

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2
    
    elif RP["PROBLEM"] == "TORO-2":

        for i in range(1, RP["N"]+1):
            if RP["grid"][i-1] < 0.5:
                U[0][i] = 1.0
                u = -2.0
                p = 0.4
            else:
                U[0][i] = 1.0
                u = 2.0
                p = 0.4

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2

    elif RP["PROBLEM"] == "TORO-3":

        for i in range(1, RP["N"]+1):
            if RP["grid"][i-1] < 0.5:
                U[0][i] = 1.0
                u = 0.0
                p = 1000
            else:
                U[0][i] = 1.0
                u = 0.0
                p = 0.01

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2

    elif RP["PROBLEM"] == "TORO-4":

        for i in range(1, RP["N"]+1):
            if RP["grid"][i-1] < 0.4:
                U[0][i] = 5.99924
                u = 19.5975
                p = 460.894
            else:
                U[0][i] = 5.99242
                u = -6.19633
                p = 46.0950

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2

    elif RP["PROBLEM"] == "LAX":

        for i in range(1, RP["N"]+1):
            if RP["grid"][i-1] < 0:
                U[0][i] = 0.445
                u = 0.698
                p = 3.528
            else:
                U[0][i] = 0.5
                u = 0
                p = 0.571

            U[1][i] = U[0][i]*u
            U[2][i] = (p/(RP["GAMMA"]-1)) + 0.5*U[0][i]*u**2
    
    else:
        print(f"---ERROR--- Please enter correct problem for simulation---")

    return U

def InitPrimitiveVariables(RP: dict) -> list[list]:
    """
        Function to initialize the primitive variables
        u -> velocity
        p -> pressure
        a -> acoustic speed
    """

    u = [0 for i in range(RP["N"]+2)]
    p = [0 for i in range(RP["N"]+2)]
    a = [0 for i in range(RP["N"]+2)]

    for i in range(1, RP["N"]+1):
        u[i] = RP["U"][1][i]/RP["U"][0][i]
        p[i] = (RP["GAMMA"]-1)*(RP["U"][2][i]-0.5*RP["U"][0][i]*u[i]**2)
        a[i] = math.sqrt(RP["GAMMA"]*p[i]/RP["U"][0][i])

    return [u, p, a]