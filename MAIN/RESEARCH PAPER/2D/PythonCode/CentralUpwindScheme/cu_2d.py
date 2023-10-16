"""
    Module to initialize the 2D Riemann Problem
"""

import numpy as np

from CentralUpwindScheme.initlaize_riemann_variables import(
    get_initial_conditions,
    get_one_dimensional_grids,
    get_conserved_variables
)

from CentralUpwindScheme.plot import plot

import CentralUpwindScheme.constants as CTS
import CentralUpwindScheme.utility as UTL

from CentralUpwindScheme.cu_flux import *

class CentralUpwindScheme:
    """
        Class representing the Riemann problem
        and run the Central Upwind Scheme
    """

    def __str__(self):
        """
            Special function to represent the class object
        """
        return f"<{self.data['problem']}({self.data['mode']})>"

    def __init__(self, problem_name, run_mode):
        """
            Constructor to initialize the class variables
        """

        # Initializing the problem variables in a dictionary
        self.data = {
            "problem" : problem_name,
            "mode" : run_mode
        }

        # Initialize the Starting conditions for the Central Upwind Scheme to run
        # get the initial conditions and merge them to the class variable dictionary
        initial_conditions = get_initial_conditions(self.data["mode"], self.data["problem"])
        self.data = UTL.merge_dicts(self.data, initial_conditions)

        # Function to initialize the One dimensional grids for the problem
        oned_grids = get_one_dimensional_grids(self.data)
        self.data = UTL.merge_dicts(self.data, oned_grids)

        # Function to initialize the conserved variable grids for the problem
        # initialize the main grid and then extend the boundaries with the given conditions
        conserved_variables = get_conserved_variables(self.data)
        self.data = UTL.merge_dicts(self.data, conserved_variables)

        self.run_cu_partial()

    def run_cu_partial(self) -> None:
        """
            Function to run the Partial Central Upwind Scheme storing the state
            at the first and the last iteration only used for static plots only
        """
        
        F1,F2,F3,F4,G1,G2,G3,G4 = get_cu_flux(self.data)

        


        
    def run_cu_complete(self) -> None:
        """
            Function to run the Complete Central Upwind Scheme storing the state 
            at every 5th iteration for animating the evolution
        """
        






