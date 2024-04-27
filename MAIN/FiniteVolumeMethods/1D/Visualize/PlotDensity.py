"""
    Module contains functions to plot the density for the simulation
"""

import matplotlib.pyplot as plt

def plot_density(RP: dict) -> None:
    """
        Function to plot the density
    """

    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    ax.plot(RP["grid"], RP["U"][0][1:-1], label=f"Density-{RP['PROBLEM']}-{RP['MODE']}", color="red")
    ax.set_title(RP["PROBLEM"])
    ax.legend()

    plt.show()

def plot_density_withref(RP: dict, RPref: dict) -> None:
    """
        Function to plot the density with reference plot
    """
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()

    ax.plot(RP["grid"], RP["U"][0][1:-1], label=f"Density-{RP['PROBLEM']}-{RP['MODE']}", color="red")
    ax.plot(RPref["grid"], RPref["U"][0][1:-1], label=f"Density-{RPref['PROBLEM']}-{RPref['MODE']}", color="black")
    ax.set_title(RP["PROBLEM"])
    ax.legend()

    fig.savefig(f"{RP['PROBLEM']}-Density-LF1.png")

    plt.show()
