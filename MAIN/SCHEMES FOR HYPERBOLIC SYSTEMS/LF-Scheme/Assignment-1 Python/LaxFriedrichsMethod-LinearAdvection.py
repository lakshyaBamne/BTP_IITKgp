# Using the Lax-Friedrich scheme for modeling the evolution of Linear Advection Hyperbolic PDE
###############################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
from itertools import count

###############################################################################################
# Taking appropriate input from the user about the initial conditions and properties
###############################################################################################

# domain of observation
domainEndPoints = tuple(map(float, input("Enter the staring and ending points for the Domain (separated by space) : ").split()))

# fractional ratio
timeSpaceRatio = float(input("Enter the Time Space ratio (\u0394t/\u0394x) : "))

# number of equal intervals in the domain
numIntervals = int(input("Enter the number of equal intervals to be made in the domain : "))

# Final time at which we have to stop the evolution
finalTime = float(input("Enter the final time at which to stop evolution : "))

###############################################################################################
# Calculating important values to be used in the modelling
###############################################################################################

# Length of one space interval
spaceIntervalLen = ( domainEndPoints[1] - domainEndPoints[0] ) / numIntervals

# Length of one time interval
timeIntervalLen = spaceIntervalLen * timeSpaceRatio

# number of time steps for which we want to run the evolution (iteration)
numTimeSteps = int(finalTime / timeIntervalLen)

###############################################################################################
# Initializing the domain and corresponding values at time t=0
###############################################################################################

# number of discrete points in the domain x comes from the number of intervals
numPointsDomain = numIntervals + 1

x = np.linspace(domainEndPoints[0], domainEndPoints[1], num=numPointsDomain, endpoint=True)
# y = np.copy(x)

y = np.sin(x)

# for i in range(len(x)):
#     if x[i]<(1/3) and x[i]>(-1/3):
#         # value should be 1
#         y[i] = 1.0
#     else:
#         # value should be 0
#         y[i] = -1.0

###############################################################################################
# Now we can start the iterations to observe the evolution of the Hyperbolic PDE overtime
###############################################################################################

# creating two arrays to store the values as the equation evolves
y_final = np.copy(y) 
y_temp = np.copy(y)

y_all = []

for i in range(numTimeSteps):
    # outer loop would run equal to number of time steps
    for j in range(1, len(x)-1):
        # two boundary points are initialized separately using the periodic boundary condition
        # for rest of the points the formula from LF scheme is used directly
        # y_final[j] = ((y_temp[j-1] + y_temp[j+1])/2) + (timeSpaceRatio/2)*(y_temp[j+1]-y_temp[j-1])
        y_final[j] = ((y_temp[j-1] + y_temp[j+1])/2) - (timeSpaceRatio/2)*(y_temp[j+1]**2/2-y_temp[j-1]**2/2)
    y_final[0] = y_final[-2]
    y_final[-1] = y_final[1]

    y_temp = np.copy(y_final)

    if i%5==0:
        y_all.append(y_final)

# now we have the initial and the final state of the equation in arrays y and y_final respectively
# we want to plot both on the same graph to observe the evolution of the equation

myPlot, ax = plt.subplots()

ax.plot(x,y,label='Initial State')
ax.plot(x,y_final, label='Final State')

ax.legend()

plt.show()

# code to animate the plot, taking the values of points at each iteration as one frame
counter = count(0,1)
def animate_vars(x, y) -> None:
    """
        Animate the Lax Friedrichs Scheme
    """
    fig, ax = plt.subplots(figsize=(10,10))
    ax.axis([-3,3,-1.5,1.5])

    idx = next(counter)

    def updater(i):
        """
            updater function to update the animation
        """
        ax.clear()
        ax.plot(x, y[idx], color="red")

    ani = FuncAnimation(fig=fig, func=updater, interval=1,frames=len(y))
    ani.save(f"burger.gif", fps=360)

    plt.show()

# animate_vars(x, y_all)

