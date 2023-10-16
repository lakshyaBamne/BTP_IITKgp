"""
    Module to Initialize the variables for the Riemann Problems
"""

import numpy as np

import CentralUpwindScheme.constants as CTS
import CentralUpwindScheme.extend_cells as EXC
import CentralUpwindScheme.plot as PLT

def get_initial_conditions(mode: str, problem: str) -> dict:
    """
        Function returns the initial conditions for the Riemann problem
    """

    data = {}

    if problem == "MCW" and mode == "CU":
        data["domainX"] = [-0.2, 0.2]
        data["domainY"] = [0, 0.8]
        data["dx"] = float(1/250)
        data["dy"] = float(1/250)
        data["Nx"] = int( ( data["domainX"][1]-data["domainX"][0] ) / data["dx"] )
        data["Ny"] = int( ( data["domainY"][1]-data["domainY"][0] ) / data["dy"] )
        data["time"] = [0, 2]
        data["t"] = 0
        data["dt"] = 0
        data["bc"] = {
            "N" : "FREE",
            "S" : "FREE",
            "E" : "FREE",
            "W" : "FREE"
        }
    elif problem == "MCW" and mode == "REF":
        data["domainX"] = [-0.2, 0.2]
        data["domainY"] = [0, 0.8]
        data["dx"] = float(1/1000)
        data["dy"] = float(1/1000)
        data["Nx"] = int( ( data["domainX"][1]-data["domainX"][0] ) / data["dx"] )
        data["Ny"] = int( ( data["domainY"][1]-data["domainY"][0] ) / data["dy"] )
        data["time"] = [0, 2]
        data["t"] = 0
        data["dt"] = 0
        data["bc"] = {
            "N" : "FREE",
            "S" : "FREE",
            "E" : "FREE",
            "W" : "FREE"
        }
    elif problem == "2DRP-CFG3" and mode == "CU":
        data["domainX"] = [0, 1.2]
        data["domainY"] = [0, 1.2]
        # data["dx"] = float(3/2500)
        # data["dy"] = float(3/2500)
        data["dx"] = float(3/500)
        data["dy"] = float(3/500)
        data["Nx"] = int( ( data["domainX"][1]-data["domainX"][0] ) / data["dx"] )
        data["Ny"] = int( ( data["domainY"][1]-data["domainY"][0] ) / data["dy"] )
        data["time"] = [0, 1]
        data["t"] = 0
        data["dt"] = 0
        data["bc"] = {
            "N" : "FREE",
            "S" : "FREE",
            "E" : "FREE",
            "W" : "FREE"
        }
    elif problem == "2DRP-CFG3" and mode == "REF":
        data["domainX"] = [0, 1.2]
        data["domainY"] = [0, 1.2]
        data["dx"] = float(3/10000)
        data["dy"] = float(3/10000)
        data["Nx"] = int( ( data["domainX"][1]-data["domainX"][0] ) / data["dx"] )
        data["Ny"] = int( ( data["domainY"][1]-data["domainY"][0] ) / data["dy"] )
        data["time"] = [0, 1]
        data["t"] = 0
        data["dt"] = 0
        data["bc"] = {
            "N" : "FREE",
            "S" : "FREE",
            "E" : "FREE",
            "W" : "FREE"
        }
    else:
        print(f"---ERROR--- Please enter correct Problem/Mode ---")

    return data

def get_one_dimensional_grids(data: dict) -> dict:
    """
        Function to initialize the computational grid for the given problem
    """

    # first extract the required variables from the data dictionary
    x = data["domainX"]
    y = data["domainY"]

    dx = data["dx"]
    dy = data["dy"]

    Nx = int( (x[1]-x[0])/dx )
    Ny = int( (y[1]-y[0])/dy )

    # initialize the grids
    computational_grid = {
        "gridX" : [ x[0]+dx*(i+0.5) for i in range(Nx) ],
        "gridY" : [ y[0]+dy*(i+0.5) for i in range(Ny) ]
    }

    return computational_grid

def get_conserved_variables(data: dict) -> dict:
    """
        Function to initialize the conserved grids for the given problem
        -> returns a dict of numpy arrays one for each conserved variable
    """
    x = data["domainX"]
    y = data["domainY"]

    dx = data["dx"]
    dy = data["dy"]

    Nx = int( (x[1]-x[0])/dx )
    Ny = int( (y[1]-y[0])/dy )

    gridX = data["gridX"]
    gridY = data["gridY"]
    
    # conserved variable matrices
    rho = np.zeros((Nx+2, Ny+2))
    mx = np.zeros((Nx+2, Ny+2))
    my = np.zeros((Nx+2, Ny+2))
    E = np.zeros((Nx+2, Ny+2))

    for i in range(1,Nx+1):
        for j in range(1,Ny+1):

            if gridX[i-1]>1 and gridY[j-1]>1:
                rho[i][j] = 1.5
                u = 0
                v = 0
                p = 1.5
            elif gridX[i-1]<1 and gridY[j-1]>1:
                rho[i][j] = 0.5323
                u = 1.206
                v = 0
                p = 0.3
            elif gridX[i-1]<1 and gridY[j-1]<1:
                rho[i][j] = 0.138
                u = 1.206
                v = 1.206
                p = 0.029
            else:
                rho[i][j] = 0.5323
                u = 0
                v = 1.206
                p = 0.3
            
            mx[i][j] = rho[i][j]*u
            my[i][j] = rho[i][j]*v
            E[i][j] = ( p/(CTS.GAMMA-1) ) + ( rho[i][j]/2 )*( u**2 + v**2 )

    # extend the matrices according to the boundary conditions for the problem
    EXC.extend_matrix(rho, data["bc"], "Density")
    EXC.extend_matrix(mx, data["bc"], "MomentumX")
    EXC.extend_matrix(my, data["bc"], "MomentumY")
    EXC.extend_matrix(E, data["bc"], "Energy")

    conserved_variables = {
        "rho" : rho,
        "mx" : mx,
        "my" : my,
        "E" : E
    }

    return conserved_variables