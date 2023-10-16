"""
    Module to plot the results for the Central Upwind Scheme
"""

import numpy as np
import matplotlib.pyplot as plt

def plot(data: dict):
    """
        Function to plot the results
    """

    fig, axd = plt.subplot_mosaic(
        [
            [ 'Density' , 'Energy' ],
            [ 'MomentumX' , 'MomentumY' ]
        ],
        figsize=(10,10),
        layout="constrained"
    )

    rho = axd["Density"].imshow(
        data["rho"][1:-1, 1:-1],
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ data["domainX"][0], data["domainX"][1], data["domainY"][0], data["domainY"][1] ],    
    )
    axd["Density"].set_title("Density")
    fig.colorbar(rho, ax=axd["Density"])

    mx = axd["MomentumX"].imshow(
        data["mx"][1:-1, 1:-1],
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ data["domainX"][0], data["domainX"][1], data["domainY"][0], data["domainY"][1] ],    
    )
    axd["MomentumX"].set_title("Momentum X")
    fig.colorbar(mx, ax=axd["MomentumX"])

    my = axd["MomentumY"].imshow(
        data["my"][1:-1, 1:-1],
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ data["domainX"][0], data["domainX"][1], data["domainY"][0], data["domainY"][1] ],    
    )
    axd["MomentumY"].set_title("Momentum Y")
    fig.colorbar(my, ax=axd["MomentumY"])

    E = axd["Energy"].imshow(
        data["E"][1:-1, 1:-1],
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        extent=[ data["domainX"][0], data["domainX"][1], data["domainY"][0], data["domainY"][1] ],    
    )
    axd["Energy"].set_title("Energy")
    fig.colorbar(E, ax=axd["Energy"])

    fig.suptitle(f"{data['problem']} ({data['mode']})")

    plt.show()

def contour(data:dict):
    """
        Function to make a contour plot for the results for the 2D CU Scheme
    """



    