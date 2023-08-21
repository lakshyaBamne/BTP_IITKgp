# Python file to take input given by the C++ program and make a plot

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

def get_grid() -> list:
    """
        Function to open the file "ComputationalGrid.txt" and initialize a list
        which can act as the independent coordinate in the plot
    """

    file_name = "ComputationalDomain.txt"

    with open(file=file_name, mode='r') as f:
        grid_str = f.readline()

    grid = list(map(float, grid_str.split()))

    return grid 

def get_cons_vars(file_name) -> list[list]:
    """
        Function to extract the Conserved variables into a 2D list
    """
    cons_var = []

    with open(file=file_name, mode='r') as f:
        data = f.readlines()

    for line in data:
        one_line = list(map(float, line.split()))
        cons_var.append(one_line)

    return cons_var

def plot_cons_vars(grid, cons_var) -> None:
    """
        Function to pring the Initial and final values of the conserved variable only
        grid -> independent computational grid
        var -> dependent grid for the Conserved variable to be plotted
    """
    var_initial = cons_var[0][1:-1]
    var_final = cons_var[1][1:-1]

    fig, ax = plt.subplots(figsize=(10,10))
    
    ax.axis([-0.1, 1.1, -0.5, 4])

    ax.plot( grid, var_initial, label="Initial")
    ax.plot( grid, var_final, label="Final", linestyle="--")

    ax.legend()

    plt.show()


#! First we need to open the files to extract useful information to plot

#! 1) Get the computational grid from the file "ComputationalDomain.txt"
grid = get_grid()

#! 2) Get the conserved variables from their corresponding files
density = get_cons_vars("Density.txt")
momentum = get_cons_vars("Momentum.txt")
energy = get_cons_vars("Energy.txt")

#! Now we need to make plots to visualize these results

plot_cons_vars(grid, density)
plot_cons_vars(grid, momentum)
plot_cons_vars(grid, energy)
