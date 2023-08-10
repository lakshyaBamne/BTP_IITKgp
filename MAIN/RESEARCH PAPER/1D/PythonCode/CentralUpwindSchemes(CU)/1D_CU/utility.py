# Initial script implementing the CU numerical scheme for 1-D

import numpy as np
import matplotlib.pyplot as plt

#! GLOBAL VARIABLES
GAMMA = 1.4
THETA = 1.3

class Utility:
    """
        Utility class to encapsulate all the useful functions
    """

    def minmod(self, a, b):
        """
            minmod function for two arguments
        """
        if a>0 and b>0:
            return min(a,b)
        elif a<0 and b<0:
            return max(a, b)
        else:
            return 0

    def minmod(self, a, b, c):
        """
            minmod function for three arguments
        """
        if a>0 and b>0 and c>0:
            return min(a, b, c)
        elif a<0 and b<0 and c<0:
            return max(a, b, c)
        else:
            return 0
        
class FluxJacobian:
    """
        Class encapsulates the methods related to handling flux jacobian matrix
    """

    def find_jacobian(self, u1: float, u2: float, u3: float, gamma: float):
        """
            class constructor is used to initialize a jacobian matrix
            for the given variable values

            u1 -> Density

            u2 -> Momentum
            
            u3 -> Total energy
        """

        flux_jacobian = np.zeros((3,3))

        # row-1 of the Flux Jacobian Matrix
        flux_jacobian[0][0] = 0
        flux_jacobian[0][1] = 1
        flux_jacobian[0][2] = 0

        # row-2 of the Flux Jacobian Matrix
        flux_jacobian[1][0] = ( (gamma-3)/2 )*( u2**2/u1**2 )
        flux_jacobian[1][1] = ( (3-gamma) )*( u2/u1 )
        flux_jacobian[1][2] = gamma-1

        # row-3 of the Flux Jacobian Matrix
        flux_jacobian[2][0] = ( (gamma-1)*( u2**3/u1**3 ) ) - ( gamma*( (u2*u3)/u1**2 ) )
        flux_jacobian[2][1] = ( gamma*(u3/u1) ) - ( (3/2)*(gamma-1)*( u2**2/u1**2 ) )
        flux_jacobian[2][2] = gamma * ( u2/u1 )

        return flux_jacobian
    
    def max_eigenvalue(self, flux_jacobian):
        """
            Method returns the maximum eigenvalue of the Flux Jacobian
        """

        return max(np.linalg.eig(flux_jacobian)[0])

    def min_eigenvalue(self, flux_jacobian):
        """
            Method returns the minimum eigenvalue of the Flux Jacobian
        """

        return min(np.linalg.eig(flux_jacobian)[0])
    
    #! we need two functions to find the upper/lower bounds for the one sided
    #! local speeds of propagation

class Grid:
    """
        Class to contain the methods related to the Finite Volume Grid
    """

    # endpoints of the Computational Domain
    DOMAIN_START = 0
    DOMAIN_END = 1

    # grid information
    NUM_GRID_POINTS = None 
    DELTA_X = None
    
    #! numpy array to store the x-coordinates of the Finite Volume Grid
    GRID = None 

    def make_grid(self):
        """
            Methods makes a Finite Volume Grid following the steps :-

            1) The Grid points are initialized with ghost indices
            
            2) Initial Conditions are 
        """

        self.NUM_GRID_POINTS = int(input("Enter Number of Points on the Finite Volume Grid : "))
        self.DELTA_X = (self.DOMAIN_END-self.DOMAIN_START)/self.NUM_GRID_POINTS

        # initializing the grid of Uniform Spacing 
        # 2 extra grid points are initialized, these are for GHOST VALUES
        self.GRID = np.linspace(
            self.DOMAIN_START - self.DELTA_X, 
            self.DOMAIN_END + self.DELTA_X, 
            num=self.NUM_GRID_POINTS+2, 
            endpoint=True
        )

        # logging message
        print("\n____LOG____ GRID CREATION SUCCESSFULL\n")

    def init_variables(self, density, velocity, pressure, momemtum, energy):
        """
            Method to assign values to primitive variables on the grid values, 
            first the main grid points are assigned values
        """

        #! Initial values for the variables => density, velocity and pressure
        #! are determined by the problem of interest
        #! currently we have defined the function for the one in the Research paper

        # first we are initializing the primitive variables
        for i in range(1, len(self.GRID)-1):
            if self.GRID[i]<=0.3:
                density[i] = 1.4
                velocity[i] = 0.1
                pressure[i] = 1
            elif self.GRID[i]>0.3:
                density[i] = 1.0
                velocity[i] = 0.1
                pressure[i] = 1

        # after initializing the primitive variables we need to initialize the conserved variables
        for i in range(1, len(self.GRID)-1):
            momemtum[i] = density[i] * velocity[i]
            energy[i] = ( pressure[i] / (GAMMA-1) ) - ( (density[i] * (velocity[i]**2)) / 2 )

        # after initializing the main grid points we need to initialize the ghost values
        # for the conserved variables
        self.init_ghost_values(density, momemtum, energy)


    def init_ghost_values(self, density, momentum, energy):
        """
            Method to assign ghost values to a grid and given variables
            using the given method

            ! we are initializing the ghost values based on the FREE BOUNDARY CONDITIONS
        """

        density[0] = density[1]
        density[-1] = density[-2]

        momentum[0] = momentum[1]
        momentum[-1] = momentum[-2]

        energy[0] = energy[1]
        energy[-1] = energy[-2]



