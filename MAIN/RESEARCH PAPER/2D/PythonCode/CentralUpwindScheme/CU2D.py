"""
    Module to initialize the 2D Riemann Problem
"""

import numpy as np

from CentralUpwindScheme.InitializeRiemannVariables import get_initial_conditions

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
        get_initial_conditions(self.data["mode"], self.data["problem"])



    

        




