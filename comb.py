import tensorflow as tf
import numpy as np
from matplotlib import pyplot as plt


def waterDiffusivity(t):
    """
    waterDiffusivity returns the diffusivity of water at a given temperature
    in mm^2/s it is accurate between 0 < t < 550 degC
    """
    #  Diffusivity source: DOI: 10.1063/1.555718
    TD_ICE = 0.1328     # mm^2/s Diffusivity T=0C  P=0.5MPa
    TD_VAP25 = 0.1456   # mm^2/s Diffusivity T=25  P=0.5MPa
    TD_VAP150 = 0.1725  # mm^2/s Diffusivity T=150 P=0.5MPa
    TD_VAP200 = 6.942   # mm^2/s Diffusivity T=200 P=0.5MPa
    TD_VAP550 = 25.58   # mm^2/s Diffusivity T=550 P=0.5MPa

    if t >= 0 and t < 25:
        return (TD_ICE - TD_VAP25) / (0 - 25) * t + TD_ICE
    if t >= 25 and t < 150:
        return (TD_VAP25 - TD_VAP150) / (25 - 150) * (t - 25) + TD_VAP25
    if t >= 150 and t < 200:
        return (TD_VAP150 - TD_VAP200) / (150 - 200) * (t - 150) + TD_VAP150
    if t >= 200 and t <= 550:
        return (TD_VAP200 - TD_VAP550) / (200 - 550) * (t - 200) + TD_VAP200
    if t < 0:
        return TD_ICE
        print("WARN: waterDiffusivity is wrong at temperature {0}".format(t))
    if t > 550:
        return (TD_VAP200 - TD_VAP550) / (200 - 550) * (t - 200) + TD_VAP200
        print("WARN: waterDiffusivity is wrong at temperature {0}".format(t))


def saveArrayt(a, t, meltTemp=0, dunits="mm", tunits="ms", step=1):
    """
    saveArrayt prints the array to a file with some formatting.
        a......... the array to be printed
        t......... the time in tunits
        meltTemp.. the ice melt temperature (for the ring)
        dunits.... the distance units
    """
    a = np.asarray(a)  # make sure there is no funny business going on...
    length = len(a)
    r = getMeltRadius(array=a, meltTemp=meltTemp)
    x = np.arange(0, length*step, step)
    y = np.arange(0, length*step, step)
    X, Y = np.meshgrid(x, y)
    circle2 = plt.Circle(
        (length/2*step, length/2*step), r*step,
        color='r', linestyle='dashed', fill=False, linewidth=2.0)
    ax = plt.gca()
    ax.add_artist(circle2)
    cp1 = plt.contourf(X, Y, a.copy(), 20)

    plt.title(
        'Contour Plot of Tempeture in ice' +
        'after t = {0:.2f}{1} and r = {2:.2f}{3}'.format(
                    t, tunits,        r*step, dunits))
    plt.xlabel('x ({0})'.format(dunits))
    plt.ylabel('y ({0})'.format(dunits))
    plt.colorbar(cp1)
    plt.savefig("out-{0:.2f}.png".format(t))
    plt.close()


def getMeltRadius(array, meltTemp):
    """Finds the melt radius given the array"""
    length = len(array)
    xc = int(length / 2)  # Centerpoint
    for y in range(0, length - 1):
        if array[xc, y] > meltTemp:
            return xc - y


if __name__ == "__main__":
    # Assumptions and Constants
    w = h = 100    # Size of plot in mm
    dx = dy = 0.1  # Position step size in mm NOTE: gridsize = w * dx if square
    steps = 1000   # Number of time steps
    t = 0          # Innitial time t_0 in ms
    D = 0.1328     # Minimum diffusivity (in this case, ice) in mm^2/s
    Umelt = 10      # Ice melt temp in C    Source:TODO
    Uhot = 500     # Hot core temp in C    Source:TODO
    Uice = -20     # Ice temp in C         Source:TODO
    r = 16         # Root radius in mm

    # Force types for TensorFlow GPU Calculations on CUDA (Nvidia Cards)
    D = np.float32(D)
    Umelt = np.float32(Umelt)
    Uhot = np.float32(Uhot)
    Uice = np.float32(Uice)

    # Helper constants
    dx2 = dx**2
    dy2 = dy**2
    dt = dx2 * dy2 / (2 * D * (dx2 + dy2))  # This is the largest stable dt
    cx, cy = w/2., h/2.  # The center coordinates of the grid
    nx, ny = int(w/dx), int(h/dy)  # The number of boxes in each dimension
    r2 = r**2

    # Create Innitial condition arrays
    U0 = Uice * np.ones((nx, ny), dtype=np.float32)  # Ice
    for i in range(nx):
        for j in range(ny):
            p2 = (i*dx-cx)**2 + (j*dy-cy)**2  # Radius^2 from center
            if p2 < r2:
                # In the root radius...
                U0[i, j] = Uhot  # Root

    saveArrayt(a=U0, t=t, meltTemp=Umelt, step=dx)  # Print initial conditions
