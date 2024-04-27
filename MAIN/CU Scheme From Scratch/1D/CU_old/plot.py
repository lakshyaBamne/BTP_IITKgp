"""
    Python file to plot the results
"""

import numpy as np
import matplotlib.pyplot as plt

def get_grid(path):
    """
        Function to get the computational grid
    """
    with open(path, "r") as f:
        grid = list(map(float, f.readline().split(' ')[:-1]))

    return grid

def get_vars(file_name):
    """
        Function to read one file containing values of a variable in a 
        single 2D vector for all the time steps
    """
    var = []

    with open(file=file_name, mode='r') as f:
        data = f.readlines()

    for line in data:
        one_line = list(map(float, line.split()))
        var.append(one_line)

    return var

def plot(titlename):
    """
        Function to make a static plot given the results in the result folder
    """
    grid_plot = get_grid(f"CU-RESULTS/ComputationalDomain.txt")
    vars_plot = get_vars(f"CU-RESULTS/density.txt")

    grid_reference = get_grid(f"CU-REF-RESULTS/ComputationalDomain.txt")
    vars_reference = get_vars(f"CU-REF-RESULTS/density.txt")

    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot(grid_plot, vars_plot[1], label="Central Upwind Scheme", color="red")
    ax.plot(grid_reference, vars_reference[1], label="Reference", color="blue")

    ax.legend()
    plt.suptitle(titlename)

    plt.savefig("PLOTS/BLW-density")

    plt.show()

# plot the results

# plot("Moving Contact Wave")
# plot("Stationary Contact Wave")
# plot("Lax Problem")
plot("Blast Wave Problem")
# plot("Shu Osher Problem")

