"""
    Module contains code to implement
    -> First Order Lax Friedrichs Scheme
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def plot_density(grid: list, U1: list, **kwargs) -> None:
    """
        Function to plot the density
    """

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    ax.plot(grid, U1, label="Density", color="red")
    ax.legend()

    if "TITLE" in kwargs.keys():
        ax.set_title(kwargs["TITLE"])

    plt.show()

#! Initial variables for the problem to solve
# computational domain
domx = (0.0, 1.0)

# number of intervals in the domain
n = 5000

# length of an individual interval
dx = float((domx[1] - domx[0]) / n)

# time variables
t = 0
dt = 0
# tfinal = 2.0 # MCW
tfinal = 0.012 # SCW
# tfinal = 0.038 # BLW
# tfinal = 0.2 # TORO-1
# tfinal = 0.15 # TORO-2
# tfinal = 0.012 # TORO-3
# tfinal = 0.035 # TORO-4

# constants used in the program
GAMMA = 1.4
CFL = 0.9 # for 1st order scheme
# CFL = 0.45 # for 2nd order scheme
THETA = 1.3
EPSILON = 1.0E-12

# PROBLEM = "MCW"
PROBLEM = "SCW"
# PROBLEM = "BLW"
# PROBLEM = "TORO-1"
# PROBLEM = "TORO-2"
# PROBLEM = "TORO-3"
# PROBLEM = "TORO-4"

BC = "FREE"
# BC = "REFLECTIVE"

#! Initialize the computational grid for the domain
grid = [domx[0]+(i+0.5)*dx for i in range(n)]

#! Initialize the conserved variables and the primitive variables
# conserved variables
u = [] # velocity
p = [] # pressure

U1 = [] # Density
U2 = [] # Momentum
U3 = [] # Energy

if PROBLEM == "MCW":
    for i in grid:
        if i<0.3:
            U1.append(1.4)
        else:
            U1.append(1.0)

        u.append(0.1)
        p.append(1.0)

        U2.append(U1[-1]*u[-1])
        U3.append( p[-1]/(GAMMA-1) + 0.5*U1[-1]*u[-1]**2 )
elif PROBLEM == "SCW":
    for i in grid:
        U1.append(1.0)
        u.append(-19.59745)

        if i<0.8:
            p.append(1000)
        else:
            p.append(0.01)

        U2.append(U1[-1]*u[-1])
        U3.append( p[-1]/(GAMMA-1) + 0.5*U1[-1]*u[-1]**2 )
elif PROBLEM=="BLW":
    for i in grid:
        U1.append(1.0)
        u.append(0.0)

        if i < 0.1:
            p.append(1000)
        elif i>=0.1 and i<=0.9:
            p.append(0.01)
        else:
            p.append(100)

        U2.append(U1[-1]*u[-1])
        U3.append( p[-1]/(GAMMA-1) + 0.5*U1[-1]*u[-1]**2 )
elif PROBLEM == "TORO-1":
    for i in grid:

        if i < 0.3:
            U1.append(1.0)
            u.append(0.75)
            p.append(1.0)
        else:
            U1.append(0.125)
            u.append(0.0)
            p.append(0.1)

        U2.append(U1[-1]*u[-1])
        U3.append( p[-1]/(GAMMA-1) + 0.5*U1[-1]*u[-1]**2 )
elif PROBLEM == "TORO-2":
    for i in grid:

        if i < 0.5:
            U1.append(1.0)
            u.append(-2.0)
            p.append(0.4)
        else:
            U1.append(1.0)
            u.append(2.0)
            p.append(0.4)

        U2.append(U1[-1]*u[-1])
        U3.append( p[-1]/(GAMMA-1) + 0.5*U1[-1]*u[-1]**2 )
elif PROBLEM == "TORO-3":
    for i in grid:

        if i < 0.5:
            U1.append(1.0)
            u.append(0.0)
            p.append(1000)
        else:
            U1.append(1.0)
            u.append(0.0)
            p.append(0.01)

        U2.append(U1[-1]*u[-1])
        U3.append( p[-1]/(GAMMA-1) + 0.5*U1[-1]*u[-1]**2 )
elif PROBLEM == "TORO-4":
    for i in grid:

        if i < 0.4:
            U1.append(5.99924)
            u.append(19.5975)
            p.append(460.894)
        else:
            U1.append(5.99242)
            u.append(-6.19633)
            p.append(46.0950)

        U2.append(U1[-1]*u[-1])
        U3.append( p[-1]/(GAMMA-1) + 0.5*U1[-1]*u[-1]**2 )

#! plot initial density
# plot_density(grid, U1, TITLE=f"Initial Density - {PROBLEM}")

def update_conserved_variables() -> None:
    """
        Function to update the conserved variables
    """

    # declare usage of global variables inside the function scope
    global U1, U2, U3, dt, t

    #! update the time step using CFL conditions
    u = [U2[i]/U1[i] for i in range(n)]
    p = [(GAMMA-1)*(U3[i]-0.5*U1[i]*u[i]**2) for i in range(n)]
    a = [math.sqrt(GAMMA*p[i]/U1[i]) for i in range(n)]

    # find amax bounding the wave speed throughout the domain
    amax = 0
    for i in range(n):
        amax = max(amax, abs(u[i])+a[i])

    dt = CFL*dx/amax

    # extend the cells for the conserved variables to get the intercell fluxes around each cell centroid
    if BC == "FREE":
        U1.insert(0, U1[0])
        U2.insert(0, U2[0])
        U3.insert(0, U3[0])

        U1.insert(-1, U1[-1])
        U2.insert(-1, U2[-1])
        U3.insert(-1, U3[-1])

        u.insert(0, u[0])
        u.insert(-1, u[-1])

        p.insert(0, u[0])
        p.insert(-1, u[-1])
    elif BC == "REFLECTIVE":
        U1.insert(0, U1[0])
        U2.insert(0, -U2[0])
        U3.insert(0, U3[0])

        U1.insert(-1, U1[-1])
        U2.insert(-1, -U2[-1])
        U3.insert(-1, U3[-1])

        u.insert(0, -u[0])
        u.insert(-1, -u[-1])

        p.insert(0, u[0])
        p.insert(-1, u[-1])

    # calculate the Flux variables for the centers of the intervals
    F1 = [U2[i] for i in range(n+2)]
    F2 = [U1[i]*u[i]**2 + p[i] for i in range(n+2)]
    F3 = [u[i]*(U3[i] + p[i]) for i in range(n+2)]

    # now we can calculate the Lax Friedrichs Flux at the intercell interfaces
    LLF1 = [ 0.5*((F1[i]+F1[i+1]) - (dx/dt)*(U1[i+1]-U1[i])) for i in range(n+1) ]
    LLF2 = [ 0.5*((F2[i]+F2[i+1]) - (dx/dt)*(U2[i+1]-U2[i])) for i in range(n+1) ]
    LLF3 = [ 0.5*((F3[i]+F3[i+1]) - (dx/dt)*(U3[i+1]-U3[i])) for i in range(n+1) ]

    U1n = list(np.zeros(n))
    U2n = list(np.zeros(n))
    U3n = list(np.zeros(n))

    # remove the first and last elements from conserved variables
    U1.pop(0)
    U1.pop()
    U2.pop(0)
    U2.pop()
    U3.pop(0)
    U3.pop()

    for i in range(n):
        U1n[i] = U1[i] - (dt/dx)*( LLF1[i+1] - LLF1[i] )
        U2n[i] = U2[i] - (dt/dx)*( LLF2[i+1] - LLF2[i] )
        U3n[i] = U3[i] - (dt/dx)*( LLF3[i+1] - LLF3[i] )

    return [U1n, U2n, U3n]

while t<tfinal:
    print(f"t={t}|dt={dt}")

    # get the updated conserved variables
    Un = update_conserved_variables()

    for i in range(n):
        U1[i] = Un[0][i]
        U2[i] = Un[1][i]
        U3[i] = Un[2][i]

    # update time
    t+=dt

plot_density(grid, U1, TITLE=f"Final Density - {PROBLEM}")
