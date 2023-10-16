# Python file to take input given by the C++ program and make a plot

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
from itertools import count

def get_grid() -> list:
    """
        Function to open the file "ComputationalGrid.txt" and initialize a list
        which can act as the independent coordinate in the plot
    """

    file_name = "ComputationalDomain.txt"

    with open(file=file_name, mode='r') as f:
        grid_str = f.readline()

    grid = list(map(float, grid_str.split()))

    return grid 

def get_vars(file_name) -> list[list]:
    """
        Function to extract the Conserved variables into a 2D list
    """
    cons_var = []

    with open(file=file_name, mode='r') as f:
        data = f.readlines()

    for line in data:
        one_line = list(map(float, line.split()))
        cons_var.append(one_line)

    return cons_var

def plot_vars(grid, cons_var) -> None:
    """
        Function to pring the Initial and final values of the conserved variable only
        grid -> independent computational grid
        var -> dependent grid for the Conserved variable to be plotted
    """
    var_initial = cons_var[0][1:-1]
    var_final = cons_var[-1][1:-1]

    fig, ax = plt.subplots(figsize=(10,10))
    
    # ax.axis([-0.5, 1.5, -0.5, 4])

    ax.plot( grid, var_initial, label="Initial", color='black')
    ax.plot( grid, var_final, label="Final", linestyle="--", color='red')

    ax.legend()

    plt.show()

def animate_vars(grid, cons_vars, name) -> None:
    """
        Function to animate the Movement of the discontinuity
    """
    fig, ax = plt.subplots(figsize=(10,10))

    ax.axis([0, 1, -0.5, 3])

    y = []

    # we need to print every line in the cons_vars 2D list one by one
    counter = count(0,1)
    def updater(i):
        ax.clear()
    
        # first plot the initial values statically
        ax.plot( grid, cons_vars[0][1:-1] , label="Initial", color='black')
    
        idx = next(counter)
        
        y.clear()
        for e in cons_vars[idx][1:-1]:
            y.append(e)

        ax.plot(grid,y, color='red', label="Updated")
        ax.legend()

    ani = FuncAnimation(fig=fig, func=updater, interval=1, frames=len(cons_vars)-2)
    ani.save(f'{name}.mp4', fps=360) # save the animation

    plt.show()        

def plot_all(grid , all_vars) -> None:
    """
        Function to plot all the variables in the same plot
    """
    plt.subplot(3,1,1)
    plt.plot( grid , all_vars[0][0][1:-1] , color='black' , label="Initial Density" )
    plt.plot( grid , all_vars[0][-1][1:-1] , color='red' , label="Updated Density")

    plt.legend()

    plt.subplot(3,1,2)
    plt.plot( grid , all_vars[1][0][1:-1] , color='black' , label="Initial Momemtum" )
    plt.plot( grid , all_vars[1][-1][1:-1] , color='blue' , label="Updated Momentum")

    plt.legend()

    plt.subplot(3,1,3)
    plt.plot( grid , all_vars[2][0][1:-1] , color='black' , label="Initial Energy" )
    plt.plot( grid , all_vars[2][-1][1:-1] , color='orange' , label="Updated Energy")

    plt.legend()

    plt.show()

def plot_both( vars_1 , vars_2 ):
    """
        Function to plot two together
    """
    file_name = "ComputationalDomain.txt"

    grid = []

    with open(file=file_name, mode='r') as f:
        data = f.readlines()

    for line in data:
        one_line = list(map(float, line.split()))
        grid.append(one_line)

    plt.plot( grid[0] , vars_1[0] , color='red', label="Numerical Scheme")
    plt.plot( grid[1] , vars_2[0] , color='blue', label="Reference")

    plt.legend()

    plt.show()    

#! First we need to open the files to extract useful information to plot

#! 1) Get the computational grid from the file "ComputationalDomain.txt"
# grid = get_grid()

#! 2) Get the conserved variables from their corresponding files
density1 = get_vars("Density1.txt")
momentum1 = get_vars("Momentum1.txt")
energy1 = get_vars("Energy1.txt")

density2 = get_vars("Density2.txt")
momentum2 = get_vars("Momentum2.txt")
energy2 = get_vars("Energy2.txt")

#! 3) Get the Primitive variavles to plot
# velocity = get_vars("Velocity2.txt")
# pressure = get_vars("Pressure2.txt")

#! Now we need to make plots to visualize these results

# all plots in the same figure
# plot_all(grid, [density, momentum, energy])

# static plots
# plot_vars(grid, density1)
# plot_vars(grid, momentum)
# plot_vars(grid, energy)
# plot_vars(grid, velocity)
# plot_vars(grid, pressure)

# animations
# animate_vars(grid, density1, "density")
# animate_vars(grid, momentum, "momentum")
# animate_vars(grid, energy, "energy")
# animate_vars(grid, velocity, "velocity")
# animate_vars(grid, pressure, "pressure")

plot_both([density1[-1][1:-1]], [density2[-1][1:-1]])


