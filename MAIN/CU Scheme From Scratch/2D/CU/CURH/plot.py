import matplotlib.pyplot as plt

def read_env(mode: str) -> dict:
    """
        Function to read the run environment for the problem
    """
    if mode == "CU":
        file_name = "env1/Environment.txt"
    else:
        file_name = "env2/Environment.txt"

    # dictionary to store the environment variables
    env = {}

    with open(file_name, "r") as f:
        lines = f.readlines()

    env["MODE"] = int(lines[0][:-1])
    env["PROBLEM"] = lines[1][:-1]
    env["domX"] = list(map( float , lines[2][:-1].split() ))
    env["domY"] = list(map( float , lines[3][:-1].split() ))
    env["time"] = list(map( float , lines[4][:-1].split() ))
    env["nx"] = int(lines[5][:-1])
    env["ny"] = int(lines[5][:-1])

    return env
    
def read_grids(mode: str) -> tuple[list]:
    """
        Function to read the computational grids for the problem
    """

    if mode == "CU":
        file_name = "result1/ComputationalGrid.txt"
    else:
        file_name = "result2/ComputationalGrid.txt"

    with open(file_name, "r") as f:
        lines = f.readlines()

    gridx = list(map(float, lines[0][1:-1].split(" ")[:-1]))
    gridy = list(map(float, lines[1][1:-1].split(" ")[:-1]))

    return (gridx, gridy)

def read_result_matrix(mode: str, env: dict) -> tuple[list[list]]:
    """
        Function to read the initial and final state of the density 
    """

    # set the correct file name
    if mode == "CU":
        file_name = "result1/Density.txt"
    else:
        file_name = "result2/Density.txt"

    # read the output file
    with open(file_name, "r") as f:
        lines = f.readlines()

    # half of the lines are initial and rest half are final
    density_initial = [list(map(float, line.split(" ")[:-1])) for line in lines[:env["ny"]]]
    density_final = [list(map(float, line.split(" ")[:-1])) for line in lines[env["ny"]:]]

    return (density_initial, density_final)

def plot(env: dict, grids: tuple[list], density: tuple[list[list]]) -> None:
    """
        Function to plot the Density initially and finally
    """

    fig, axd = plt.subplot_mosaic(
        [
            ['INITIAL', 'CU'],
        ],
        figsize=(10,10),
        layout="constrained"
    )

    # Initial heat map
    hm1 = axd["INITIAL"].imshow(
        density[0],
        cmap="jet",
        # origin="lower",
        extent=[ env["domX"][0], env["domX"][1], env["domY"][0], env["domY"][1] ],    
    )
    axd["INITIAL"].set_title("Initial conditions")
    fig.colorbar(hm1, ax=axd["INITIAL"])

    # Initial heat map
    hm2 = axd["CU"].imshow(
        density[1],
        cmap="jet",
        # origin="lower",
        extent=[ env["domX"][0], env["domX"][1], env["domY"][0], env["domY"][1] ],    
    )
    axd["CU"].set_title("Final conditions")
    fig.colorbar(hm1, ax=axd["CU"])

    fig.suptitle(f"{env['PROBLEM']} - Old Central Upwind Scheme")

    plt.savefig(f"plots/{env['PROBLEM']}_density.png")
    plt.show()


if __name__ == "__main__":
    env = read_env("CU")
    grids = read_grids("CU")
    density = read_result_matrix("CU", env)

    # plot initial and final density
    plot( env, grids, density )



