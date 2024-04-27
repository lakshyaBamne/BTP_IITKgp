import matplotlib.pyplot as plt

def read_env(path: str) -> dict:
    """
        Function to read the env for the problems
    """

    env = {}

    with open(path, "r") as f:
        lines = f.readlines()

    env["PROBLEM"] = lines[0][:-1]
    env["BC"] = lines[1][:-1]
    env["PPR"] = lines[2][:-1]
    env["FLUX"] = lines[3][:-1]
    env["N"] = int( lines[4][:-1] )
    env["domx"] = list(map(float, lines[5].split(" ")))

    return(env)

def read_grid(env: dict) -> list:
    """
        Function to read the computational grid for the problem
    """

    with open(f"GRIDS/{env['PROBLEM']}.txt", "r") as f:
        lines = f.readlines()

    grid = list(map(float, lines[0].split(" ")[1:-2]))

    return(grid)
    

def read_density(env: dict) -> tuple[list]:
    """
        Function to read the initial and final density for the problem
    """

    path = f"RESULTS/{env['PROBLEM']}-{env['FLUX']}.txt"

    with open(path, "r") as f:
        lines = f.readlines()

    initial_density = list(map(float, lines[0].split(" ")[1:-2]))
    final_density = list(map(float, lines[1].split(" ")[1:-2]))

    density = (initial_density, final_density)

    return(density)

def plot_one(env: dict, grid: list, density: tuple[list]) -> None:
    """
        Function to plot the initial and final densities
    """
    plt.plot(grid, density[0], label="Initial Density", linestyle="--", color="black")
    plt.plot(grid, density[1], label=f"{env['FLUX']}", color="red")

    plt.legend()

    plt.savefig(f"PLOTS/{env['PROBLEM']}-{env['FLUX']}.png")
    plt.show()


# read output
env = read_env("ENV/env.txt")
grid = read_grid(env)
density = read_density(env)

# plot output
plot_one(env, grid, density)

