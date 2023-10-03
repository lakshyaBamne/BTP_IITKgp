"""
    * @author - lakshya Bamne (20MA20029)
    * @supervisor - Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations
    ! Python file to plot the results
"""

import matplotlib.pyplot as plt
import numpy as np

def read_matrix(file, n):
    """
        Function to read a single matrix instance for a given text file
        -> read last n lines from a file containing many matrices
    """
    data = []
    
    with open(file, "r") as f:
        lines = f.readlines()[-n:]

    for line in lines:
        data.append( list(map(float, line.split())) )

    return data

def contour(shock_type):
    """
        Function to make a contour plot
    """
    with open(f"result/{shock_type}_ComputationalGrid.txt") as f:
        lines = f.readlines()

    x = list(map(float, lines[0].split()))
    y = list(map(float, lines[1].split()))

    X, Y = np.meshgrid(x,y)
    Z = read_matrix(f"result/{shock_type}_DENSITY.txt", 600)

    fig, ax = plt.subplots()
    CS = ax.contour(X, Y, Z , colors="black")
    ax.set_title(f"DENSITY - {shock_type}")

    plt.savefig(f"plots/{shock_type}_DENSITY(Contour).png")
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
    
    # CU heat map
    hm1 = axd["CU"].imshow(
        data_cu,
        # cmap="YlGnBu",
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
            ['CU', 'CONTOUR']
        ],
        figsize=(10,10),
        layout="constrained"
    )

    file_cu = f"result/{shock_type}_DENSITY.txt"

    data_cu = read_matrix(file_cu, 600)
    
    # CU heat map
    hm1 = axd["CU"].imshow(
        data_cu,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ -1,1,-1,1 ],    
    )
    axd["CU"].set_title(f"DENSITY PLOT - {shock_type}")
    fig.colorbar(hm1, ax=axd["CU"])

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


    plt.suptitle(f"{shock_type} - Final Density plot")

    # save the plot and show on screen
    plt.savefig(f"plots/{shock_type}_DENSITY.png")
    plt.show()

"""
    Main execution part
"""
if __name__ == "__main__":
    
    for i in  range(1, 13):
        plot(i)
