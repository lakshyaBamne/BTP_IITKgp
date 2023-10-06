"""
    * @author - lakshya Bamne (20MA20029)
    * @supervisor - Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations
    ! Python file to plot the results
"""

# setting up agg to store figures as beautiful png's
import matplotlib

# pyplot provides a simple interface for the most common functionality
import matplotlib.pyplot as plt

# numpy is a numerical processing engine usefull with all kinds of data
import numpy as np
from matplotlib.ticker import (
    AutoMinorLocator
)

def read_matrix(file, n):
    """
        Function to read a single matrix instance for a given text file
        -> read last n lines from a file containing many matrices
    """
    data = []
    
    with open(file, "r") as f:
        lines = f.readlines()[-(n+1):-1]

    for line in lines:
        data.append( list(map(float, line.split()))[1:-1] )

    return data

def contour(shock_type):
    """
        Function to make a contour plot
    """
    with open(f"result/{shock_type}_ComputationalGrid.txt") as f:
        lines = f.readlines()

    x = list(map(float, lines[0].split()))
    y = list(map(float, lines[1].split()))

    # initialize the figure object and get the default axes object
    fig = plt.figure()
    ax = fig.gca()

    # getting the data to be plotted 
    X, Y = np.meshgrid(x,y)
    Z = read_matrix(f"result/{shock_type}_DENSITY.txt", 600)
    # Z = read_matrix(f"result/{shock_type}_DENSITY.txt", 10)

    # make the contour plot
    ax.contour(
        X, Y, Z, 
        colors="red", 
    )

    # title for the axes object
    ax.set_title(f"DENSITY - {shock_type}")

    # to forcefully make the axes object a square 
    ax.set_aspect('equal', adjustable='box')

    # initialize a very light background grid for easier visualization in contour plots
    ax.grid(color='green', linestyle='-', linewidth=0.1)

    # axes ticks and tick labels for good looking graphs

    # add ticks to all 4 sides of the plot
    ax.tick_params(
        top=True,
        bottom=True,
        left=True,
        right=True,
        labeltop=True,
        labelbottom=True,
        labelleft=True,
        labelright=True
    )

    ax.set_xticks(np.arange(-1,1,0.2,dtype=float))
    ax.set_yticks(np.arange(-1,1,0.2,dtype=float))

    # adding minor ticks for better visualization
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    # save the figure so it can be visualised later
    plt.savefig(f"plots/{shock_type}_DENSITY(Contour).png")

    # command to show the figure on the screen of the user
    plt.show()

def cu(shock_type):
    """
        Function to only plot the CU plot
        -> reference takes too much time to generate for some problems
    """

    fig, axd = plt.subplot_mosaic(
        [
            ['CU']
        ],
        figsize=(10,10),
        layout="constrained"
    )
    
    file_cu = f"result/{shock_type}_DENSITY.txt"

    data_cu = read_matrix(file_cu, 600)
    # data_cu = read_matrix(file_cu, 10)
    
    # CU heat map
    hm1 = axd["CU"].imshow(
        data_cu,
        cmap="jet",
        origin="lower",
        extent=[ -1,1,-1,1 ],    
    )
    axd["CU"].set_title(f"DENSITY - {shock_type}")
    fig.colorbar(hm1, ax=axd["CU"])

    plt.savefig(f"plots/{shock_type}_DENSITY(ImagePlot).png")
    plt.show()

def plot(shock_type):
    """
        Function to plot the final density
        -> cu plot
        -> contour plot
    """
    fig, axd = plt.subplot_mosaic(
        [
            ['CU', 'CONTOUR'],
        ],
        figsize=(10,10),
        layout="constrained"
    )

    file_cu = f"result/{shock_type}_DENSITY.txt"

    data_cu = read_matrix(file_cu, 600)
    
    # CU heat map
    hm1 = axd["CU"].imshow(
        data_cu,
        cmap="jet",
        origin="lower",
        extent=[ -1,1,-1,1 ],    
    )
    axd["CU"].set_title(f"DENSITY PLOT - {shock_type}")

    axd["CU"].set_aspect('equal', adjustable='box')

    fig.colorbar(hm1, ax=axd["CU"], shrink=0.5)

    # Contour plot for final density
    with open(f"result/{shock_type}_ComputationalGrid.txt") as f:
        lines = f.readlines()

    x = list(map(float, lines[0].split()))
    y = list(map(float, lines[1].split()))

    X, Y = np.meshgrid(x,y)
    Z = read_matrix(f"result/{shock_type}_DENSITY.txt", 600)

    hm2 = axd["CONTOUR"].contour(
        X, Y, Z,colors="black"
    )
    axd["CONTOUR"].set_title(f"CONTOUR PLOT - {shock_type}")

    axd["CONTOUR"].grid(color='green', linestyle='-', linewidth=0.1)
    axd["CONTOUR"].set_aspect('equal', adjustable='box')

    plt.suptitle(f"{shock_type} - Final Density plot")

    # save the plot and show on screen
    plt.savefig(f"plots/{shock_type}_DENSITY.png")
    
    # show the plots to the screen of the user
    plt.show()

"""
    Main execution part
"""
if __name__ == "__main__":
    
    # density plot and contour plot
    # for i in  range(1, 13):
    #     plot(i)

    plot(1)

    # contour plot only
    # for i in range(1,13):
    #     contour(i)

    contour(1)
