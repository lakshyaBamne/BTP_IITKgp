import numpy as np
import matplotlib.pyplot as plt
import math

# Constants used in the program

GAMMA = 1.4
THETA = 1.3
EPSILON = 1.0e-12
CFL = 0.475

#! Problem specific parameters
PROBLEM = "MCW"
# PROBLEM = "SCW"
# PROBLEM = "BLW"
# PROBLEM = "LAX"
# PROBLEM = "SOS"

MODE = 1 # CU
# MODE = 2 # REF

BC = "FREE"
# BC = "REFLECTIVE"

#! Moving Contact Wave (MCW)
n = int(100)
time = (0, 2.0)
dom = (0, 1.0)

#! Blast Wave Problem (BLW)
# n = int(500)
# time = (0, 0.038)
# dom = (0, 1.0)

dx = float( (dom[1]-dom[0]) / n )

def minmod(a, b, c):
    """
        Non linear minmod limiter used for calculating slopes
    """
    if max(a,b,c) < 0:
        return max(a,b,c)
    elif min(a,b,c) > 0:
        return min(a,b,c)
    else:
        return 0
    
def make_grid():
    """
        Function to make the computational grid
    """
    grid = np.linspace(dom[0]+dx/2,dom[1]-dx/2, num=n, endpoint=True)
    return grid

grid = make_grid()

def extend_cells(U):
    """
        Function to extend the cells for the conserved variables
    """

    if BC=="FREE":
        for u in range(3):
            U[u,0] = U[u,1]
            U[u,-1] = U[u,-2]
    elif BC=="REFLECTIVE":
        pass
    else:
        print("---ERROR--- Select correct boundary conditions---")

def initialize_conserved_variables(grid):
    """
        Function to initialize the conserved variables for given problem
    """
    U = np.zeros((3, n+2))

    if PROBLEM=="MCW":

        for i in range(0,n):
            if grid[i]<=0.3:
                U[0,i+1] = 1.4 # rho
                U[1,i+1] = U[0,i+1]*0.1 # m
                U[2,i+1] = 1/(GAMMA-1) + 0.5*U[0,i+1]*(0.1**2)
            else:
                U[0,i+1] = 1.0
                U[1,i+1] = U[0,i+1]*0.1 # m
                U[2,i+1] = 1/(GAMMA-1) + 0.5*U[0,i+1]*(0.1**2)

    elif PROBLEM=="SCW":
        pass
    elif PROBLEM=="BLW":
        pass
    elif PROBLEM=="LAX":
        pass
    elif PROBLEM=="SOS":
        pass
    else:
        print("---ERROR--- Please select correct problem to solve---")

    print(U)

    extend_cells(U)

    return U

U0 = initialize_conserved_variables(grid)
U = initialize_conserved_variables(grid)

def get_old_cu_flux(U, update_time):
    """
        Function to calculte the CU Numerical flux (Old version)
    """

    # get the piecewise Linear reconstructions
    PLR = np.zeros((3,2,n+1))

    # calculate East and West approximations for the internal cells
    for u in range(0,3):
        for i in range(1, n+1):
            # calculate slopes
            slx = minmod( THETA*(U[u,i]-U[u,i-1]) , 0.5*(U[u,i+1]-U[u,i-1]) , THETA*(U[u,i+1]-U[u,i]) )
            PLR[u,0,i] = U[u,i] + 0.5*slx
            PLR[u,1,i-1] = U[u,i] - 0.5*slx
    
    # Extend the cells to get East approximation for the left boundary and West approximation for the right boundary
    if BC == "FREE":
        for u in range(0,3):
            PLR[u,0,0] = PLR[u,1,1]
            PLR[u,1,-1] = PLR[u,0,-2]
    elif BC=="REFLECTIVE":
        pass
    else:
        print("---ERROR--- Please select correct boundary conditons---")

    # store the local speeds of propagation
    ap = np.zeros(n+1)
    am = np.zeros(n+1)

    # calculate CU Numerical flux
    F = np.zeros((3, n+1))

    for i in range(n+1):
        # we need to calculate the primitive variables coming from surrounding grid points
        uE = PLR[1,0,i] / PLR[0,0,i]
        uW = PLR[1,1,i] / PLR[0,1,i]
        pE = (GAMMA-1)*( PLR[2,0,i] - 0.5*PLR[0,0,i]*(uE**2) )
        pW = (GAMMA-1)*( PLR[2,1,i] - 0.5*PLR[0,1,i]*(uW**2) )

        # Local speeds of propagation
        ap[i] = max(0, uE + math.sqrt( GAMMA*pE / PLR[0,0,i] ), uW + math.sqrt( GAMMA*pW / PLR[0,1,i] ))
        am[i] = min(0, uE - math.sqrt( GAMMA*pE / PLR[0,0,i] ), uW - math.sqrt( GAMMA*pW / PLR[0,1,i] ))

        # we need the flux vectors for calculating CU Numerical FLux
        f1E = PLR[1,0,i]
        f2E = PLR[1,0,i]*uE + pE
        f3E = uE*( PLR[2,0,i] + pE )

        f1W = PLR[1,1,i]
        f2W = PLR[1,1,i]*uW + pW
        f3W = uW*( PLR[2,1,i] + pW )

        dist = ap[i] - am[i]
        prod = ap[i]*am[i]

        if dist > EPSILON:
            F[0, i] = ( ap[i]*f1E - am[i]*f1W + prod*( PLR[0,1,i] - PLR[0,0,i] ) ) / dist
            F[1, i] = ( ap[i]*f2E - am[i]*f2W + prod*( PLR[1,1,i] - PLR[1,0,i] ) ) / dist
            F[2, i] = ( ap[i]*f3E - am[i]*f3W + prod*( PLR[1,1,i] - PLR[1,0,i] ) ) / dist
        else:
            F[0, i] = 0.5*( f1E + f1W )
            F[1, i] = 0.5*( f2E + f2W )
            F[2, i] = 0.5*( f3E + f3W )

    if update_time:
        global dt
        amax =0

        for i in range(n+1):
            amax = max( amax, ap[i], -am[i] )

        dt = CFL*dx/amax

    return F


def update_conserved_variables(U):
    """
        Function to update the conserved variables using Three stage SSPRK
    """

    # Stage-1 (if this stage alone is used we get the simple Euler forward difference integration)
    F1 = get_old_cu_flux(U, True)

    global dt, dx

    LAMBDA = dt/dx

    U1 = np.zeros((3,n+2))

    for u in range(3):
        for i in range(1, n+1):
            U1[u, i] = U[u, i] - LAMBDA*( F1[u, i] - F1[u, i-1] )

    extend_cells(U1)

    # update conserved variables
    for u in range(3):
        for i in range(n+2):
            U[u, i] = U1[u, i]

if __name__ == "__main__":
    T = 0
    dt = 0

    while T < time[1]:
        update_conserved_variables(U)
        T += dt
        print(f"T={T} | dt={dt}")

    print(f"Final time (T) = {T}")


    plt.plot(grid, U0[0,1:-1], color="black",linestyle="--", label="Final Density")
    plt.plot(grid, U[0,1:-1], color="blue", label="Final Density")
    plt.legend()
    plt.show()

