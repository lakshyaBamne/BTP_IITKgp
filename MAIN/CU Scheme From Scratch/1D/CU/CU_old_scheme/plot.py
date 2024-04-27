import matplotlib.pyplot as plt

def get_env(mode: str) -> dict:
    """
        Function to read the run environment for the simulations
    """
    if mode == "CU":
        file_name = "env1/Environment.txt"
    else:
        file_name = "env2/Environment.txt"

    env = {}

    with open(file_name, "r") as f:
        lines = f.readlines()

    env["MODE"] = int(lines[0][:-1])
    env["PROBLEM"] = lines[1][:-1]
    env["dom"] = list(map( float , lines[2][:-1].split() ))
    env["time"] = list(map( float , lines[3][:-1].split() ))
    env["n"] = int(lines[4][:-1])

    return env

def read_grid(mode: str) -> list:
    """
        Function to read the computational grid for the problem
    """
    if mode == "CU":
        file_name = "result1/ComputationalGrid.txt"
    else:
        file_name = "result2/ComputationalGrid.txt"

    with open(file_name, "r") as f:
        lines = f.readlines()

    grid = list(map(float, lines[0][1:-1].split(" ")[:-1]))

    return grid

def read_conserved_variables(mode: str) -> tuple[list]:
    """
        Function to read the conserved variables initial and final for plotting
    """
    
    if mode == "CU":
        file_name = "result1/density.txt"
    else:
        file_name = "result2/density.txt"

    # read the output file
    with open(file_name, "r") as f:
        lines = f.readlines()

    # half of the lines are initial and rest half are final
    density_initial = list(map(float, lines[0].split(" ")[:-1]))
    density_final = list(map(float, lines[1].split(" ")[:-1]))

    return (density_initial, density_final)

def plot(env: dict, grid: list, density: tuple[list]) -> None:
    """
        Function to plot the results for the CU Scheme
    """

    plt.plot( grid, density[0], color="black", label="Initial", linestyle="--" )
    plt.plot( grid, density[1], color="blue", label="CU" )
    plt.legend()
    
    plt.savefig(f"plots/{env['PROBLEM']}.png")
    plt.show()

if __name__ == "__main__":
    # read the output from various files for plotting
    env = get_env("CU")
    grid = read_grid("CU")
    density = read_conserved_variables("CU")

    # plot results
    plot(env, grid, density)