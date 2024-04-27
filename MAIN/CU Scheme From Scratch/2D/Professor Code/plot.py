import numpy as np
import matplotlib.pyplot as plt

def read_initial(problem: str) -> tuple:
    """
        Function returns the X,Y,RHO values for the initial conditions
    """
    initial_file = f"{problem}-rho-CURH000"

    with open(initial_file, 'r') as f1:
        lines = f1.readlines()

    X = []
    Y = []
    RHO = []

    for line in lines:
        line_items = [l for l in line.split(' ') if l!='' and l!='\n']

        X.append(float(line_items[0])) 
        Y.append(float(line_items[1])) 
        RHO.append(float(line_items[2])) 

    return (X,Y,RHO)


def read_output(problem: str) -> tuple:
    """
        Function to read the output for the given problem and plot the results for 2D CURH schemes
    """

    final_file = f"{problem}-rho-CURH001"

    with open(final_file, 'r') as f1:
        lines = f1.readlines()

    X = []
    Y = []
    RHO = []

    for line in lines:
        line_items = [l for l in line.split(' ') if l!='' and l!='\n']

        X.append(float(line_items[0])) 
        Y.append(float(line_items[1])) 
        RHO.append(float(line_items[2])) 

    return (X,Y,RHO)

def make_matrix(x: list, y: list, rho: list) -> list[list]:
    """
        Function to convert the x,y,z coordinates to a matrix to be plotted later as an image
    """

    x_one = sorted(list(set(x)))
    y_one = sorted(list(set(y)))

    xsize = len(x_one)
    ysize = len(y_one)

    rho_np = np.array(rho)

    return rho_np.reshape((ysize, xsize))


def plot(problem: str) -> None:
    """
        Function to plot the matrix as an image
    """
    
    # get the matrices to plot the data
    initial = read_initial(problem)
    final = read_output(problem)    

    matrix1 = make_matrix(initial[0],initial[1],initial[2])
    matrix2 = make_matrix(final[0],final[1],final[2])

    # start plotting
    fig, axd = plt.subplot_mosaic(
        [
            ['INITIAL', 'FINAL'],
        ],
        figsize=(10,10),
        layout="constrained"
    )

    hm1 = axd["INITIAL"].imshow(
        matrix1,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        # extent=[ env_cu["domX"][0], env_cu["domX"][1], env_cu["domY"][0], env_cu["domY"][1] ],    
    )
    axd["INITIAL"].set_title("Initial Conditions")
    fig.colorbar(hm1, ax=axd["INITIAL"])

    # Reference heat map
    hm2 = axd["FINAL"].imshow(
        matrix2,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        # extent=[ env_ref["domX"][0], env_ref["domX"][1], env_ref["domY"][0], env_ref["domY"][1] ],
    )
    axd["FINAL"].set_title("Final Conditions")
    fig.colorbar(hm2, ax=axd["FINAL"])

    fig.suptitle(f"{problem}")

    plt.savefig(f"Fortran-CURH-{problem}.png")
    plt.show()

def plot_final_only(problem: str) -> None:
    """
        Function to only plot the final output for the CURH Scheme for the given problem
    """

    final = read_output(problem)
    matrix = make_matrix(final[0],final[1],final[2])

    # start plotting
    fig, axd = plt.subplot_mosaic(
        [
            ['FINAL'],
        ],
        figsize=(10,10),
        layout="constrained"
    )

    hm2 = axd["FINAL"].imshow(
        matrix,
        # cmap="YlGnBu",
        cmap="jet",
        origin="lower",
        # extent=[ env_ref["domX"][0], env_ref["domX"][1], env_ref["domY"][0], env_ref["domY"][1] ],
    )
    axd["FINAL"].set_title(problem)
    fig.colorbar(hm2, ax=axd["FINAL"])

    plt.savefig(f"Fortran-CURH-{problem}-FinalOnly.png")
    plt.show()

#! Main program starts from here
problems = [
    # "MCW", 
    "CFG3", 
    # "EXP", 
    # "IMP1", # with 600x600 cells
    # "IMP2", # with 800x800 cells 
    # "KHI1", # at t=1
    # "KHI2", # at t=2.5
    # "KHI3", # at t=4
    # "RTI"
]

for problem in problems:
    # plot initial and final conditions
    plot(problem)

    # plot only final conditions
    # plot_final_only(problem)


# Plot Initial and Final Conditions
# plot("MCW")
# plot("CFG3")
# plot("EXP")
# plot("IMP1")
# plot("IMP2")
# plot("KHI1")
# plot("KHI2")
# plot("KHI3")
# plot("RTI")

# Plot only Final Conditions
# plot_final_only("MCW")
# plot_final_only("CFG3")
# plot_final_only("EXP")
# plot_final_only("IMP1")
# plot_final_only("IMP2")
# plot_final_only("KHI1")
# plot_final_only("KHI2")
# plot_final_only("KHI3")
# plot_final_only("RTI")



