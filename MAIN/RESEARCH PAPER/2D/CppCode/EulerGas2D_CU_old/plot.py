
#* @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
#* 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
#* Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure
import matplotlib.cm as cm
from matplotlib.ticker import AutoMinorLocator
import numpy as np

def get_env(mode):
    """
        Function to import the environment in which the problem is running
        -> MODE
        -> PROBLEM
        -> [xi , xf]
        -> [yi , yf]
        -> [ti , tf]
        -> Nx
        -> Ny
    """
    if mode=="CU":
        file = "env1/Environment.txt"
    else:
        file = "env2/Environment.txt"

    # dictionary storing the environment variables
    env = {}

    # open the environment file
    with open(file, "r") as f:
        lines = f.readlines()
    
    env["MODE"] = lines[0][:-1]
    env["PROBLEM"] = lines[1][:-1]
    env["domX"] = list(map( float , lines[2].split() ))
    env["domY"] = list(map( float , lines[3].split() ))
    env["time"] = list(map( float , lines[4].split() ))
    env["Nx"] = int(lines[5])
    env["Ny"] = int(lines[6])

    return env

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

def read_matrix_gen(file, n1, n2):
    """
        Function to read single matrix instance for a given text file
        -> read lines between n1 and n2
    """
    data = []

    with open(file, "r") as f:
        lines = f.readlines()[n1:n2]

    for line in lines:
        data.append( list(map(float, line.split())) )

    return data

def plot(var, env_cu, env_ref):
    """
        Function to plot the heat maps for a variable for both CU and reference
    """

    fig, axd = plt.subplot_mosaic(
        [
            ['CU', 'REF'],
        ],
        figsize=(10,10),
        layout="constrained"
    )
    
    file_cu = "result1/" + var + ".txt"
    file_ref = "result2/" + var + ".txt"

    data_cu = read_matrix(file_cu, env_cu["Ny"])
    data_ref = read_matrix(file_ref, env_ref["Ny"])

    # CU heat map
    hm1 = axd["CU"].imshow(
        data_cu,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ env_cu["domX"][0], env_cu["domX"][1], env_cu["domY"][0], env_cu["domY"][1] ],    
    )
    axd["CU"].set_title("Central Upwind Scheme")
    fig.colorbar(hm1, ax=axd["CU"])

    # Reference heat map
    hm2 = axd["REF"].imshow(
        data_ref,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ env_ref["domX"][0], env_ref["domX"][1], env_ref["domY"][0], env_ref["domY"][1] ],
    )
    axd["REF"].set_title("Reference")
    fig.colorbar(hm2, ax=axd["REF"])

    fig.suptitle(f"{env_cu['PROBLEM']} {var}")

    plt.savefig(f"plots/{env_cu['PROBLEM']}_{var}.png")
    plt.show()

def plot_cu(var, env_cu):
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
    
    file_cu = "result1/" + var + ".txt"

    data_cu = read_matrix(file_cu, env_cu["Ny"])
    
    # CU heat map
    hm1 = axd["CU"].imshow(
        data_cu,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ env_cu["domX"][0], env_cu["domX"][1], env_cu["domY"][0], env_cu["domY"][1] ],    
    )
    axd["CU"].set_title("Central Upwind Scheme")
    fig.colorbar(hm1, ax=axd["CU"])

    fig.suptitle(f"{env_cu['PROBLEM']} {var}")

    plt.savefig(f"plots/{env_cu['PROBLEM']}_{var}(CU_only).png")
    plt.show()

def plot_initial_final(var, env_cu):
    """
        Function to plot initial and final configurations for CU plot
    """

    fig, axd = plt.subplot_mosaic(
        [
            ['INI','CU'],
            ['INI','CU']
        ],
        figsize=(10,10),
        layout="constrained"
    )
    
    file_cu = "result1/" + var + ".txt"

    data_ini = read_matrix_gen(file_cu, 0, env_cu["Ny"])
    data_cu = read_matrix(file_cu, env_cu["Ny"])
    
    # CU heat map
    hm1 = axd["INI"].imshow(
        data_ini,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ env_cu["domX"][0], env_cu["domX"][1], env_cu["domY"][0], env_cu["domY"][1] ],    
    )
    axd["INI"].set_title("Central Upwind Scheme")
    fig.colorbar(hm1, ax=axd["INI"])

    # CU heat map
    hm2 = axd["CU"].imshow(
        data_cu,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ env_cu["domX"][0], env_cu["domX"][1], env_cu["domY"][0], env_cu["domY"][1] ],    
    )
    axd["CU"].set_title("Central Upwind Scheme")
    fig.colorbar(hm2, ax=axd["CU"])

    fig.suptitle(f"{env_cu['PROBLEM']} {var}")

    plt.savefig(f"plots/{env_cu['PROBLEM']}_{var}(CU_only).png")
    plt.show()

def plot_new(var, env_cu):
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

    # file_cu = f"result/{shock_type}_DENSITY.txt"
    file_cu = "result1/" + var + ".txt"

    data_cu = read_matrix(file_cu, 600)
    
    # CU heat map
    hm1 = axd["CU"].imshow(
        data_cu,
        cmap="jet",
        origin="lower",
        extent=[ -1,1,-1,1 ],    
    )
    axd["CU"].set_title(f"{var} plot")

    axd["CU"].set_aspect('equal', adjustable='box')

    fig.colorbar(hm1, ax=axd["CU"], shrink=0.5)

    # Contour plot for final density
    with open(f"result1/ComputationalGrid.txt") as f:
        lines = f.readlines()

    x = list(map(float, lines[0].split()))
    y = list(map(float, lines[1].split()))

    X, Y = np.meshgrid(x,y)
    Z = read_matrix(f"result1/Density.txt", 600)

    hm2 = axd["CONTOUR"].contour(
        X, Y, Z,colors="black"
    )
    axd["CONTOUR"].set_title(f"CONTOUR PLOT - {var}")

    axd["CONTOUR"].grid(color='green', linestyle='-', linewidth=0.1)
    axd["CONTOUR"].set_aspect('equal', adjustable='box')

    plt.suptitle(f"Final Plot for {var}")

    # save the plot and show on screen
    plt.savefig(f"plots/DensityContourPlot_{var}.png")
    
    # show the plots to the screen of the user
    plt.show()

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

#! Main plotting starts
if __name__ == "__main__":

    # get environment
    env_cu = get_env("CU")
    # env_ref = get_env("REF")

    # # plot conserved variables
    # plot("Density", env_cu, env_ref)
    # plot("MomentumX", env_cu, env_ref)
    # plot("MomentumY", env_cu, env_ref)
    # plot("Energy", env_cu, env_ref)
    
    # # plot primitive variables
    # plot("Pressure", env_cu, env_ref)
    # plot("VelocityX", env_cu, env_ref)
    # plot("VelocityY", env_cu, env_ref)


    # for plotting only the CU plot (fast)
    # plot_cu("Density", env_cu)
    plot_initial_final("Density", env_cu)
    plot_new("Density", env_cu)
