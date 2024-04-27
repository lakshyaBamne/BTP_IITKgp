import matplotlib.pyplot as plt

def read_grid(path: str) -> list:
    """
        Function to read the computational grid
    """

    with open(path, "r") as f:
        lines = f.readlines()

    x = list(map(float, lines[0].split(" ")[1:-2]))

    return x

def read_density(path: str) -> tuple[list]:
    """
        Function to read the initial and final density
    """

    with open(path, "r") as f:
        lines = f.readlines()

    initial_density = list(map(float, lines[0].split(" ")[1:-2]))
    final_density = list(map(float, lines[1].split(" ")[1:-2]))

    return (initial_density, final_density)

def plot(grid: list, density: tuple[list]) -> None:
    """
        Function to plot the initial and final density
    """

    plt.plot(grid, density[0], label="t=0", linestyle="--", color="black")
    plt.plot(grid, density[1], label="t=tn", color="red")
    plt.legend()

    plt.show()

grid = read_grid("grid.txt")
density = read_density("density.txt")

plot(grid, density)
