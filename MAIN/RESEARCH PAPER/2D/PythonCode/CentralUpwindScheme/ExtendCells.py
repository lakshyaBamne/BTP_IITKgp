"""
    Module contains functions to extend the boundaries of the grid

    Boundary Conditions
    -> Neumann (Free)
    -> Solidwall
    -> Reflective
    -> Periodic
"""

def extend_matrix(matrix: list[list], bd: dict, variable: str) -> None:
    """
        Function to extend the rows and columns of a matrix
        base on the boundary conditions
    """

    

