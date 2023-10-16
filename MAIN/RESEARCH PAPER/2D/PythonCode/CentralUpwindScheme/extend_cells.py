"""
    Module contains functions to extend the boundaries of the grid

    Boundary Conditions
    -> Neumann (Free)
    -> Solidwall
    -> Reflective
    -> Periodic
"""

def extend_matrix(matrix: list[list], bc: dict, variable: str) -> None:
    """
        Function to extend the rows and columns of a matrix
        base on the boundary conditions
    """

    # Free boundary conditions on all the four boundaries
    # -> Moving Contact Wave (MCW)
    # -> 2Dimensional Riemann Problem in all 19 configurations (2DR)
    if bc["N"]=="FREE" and bc["S"]=="FREE" and bc["E"]=="FREE" and bc["W"]=="FREE":
        for i in range(len(matrix)):
            matrix[i][0] = matrix[i][1]
            matrix[0][i] = matrix[1][i]
            matrix[i][-1] = matrix[i][-2]
            matrix[-1][i] = matrix[-2][i]
    


