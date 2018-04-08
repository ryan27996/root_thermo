import numpy as np
import matplotlib.pyplot as plt
# import pdb


def waterDiffusivity(t):
    """
    waterDiffusivity returns the diffusivity of water at a given temperature
    `t` in mm^2/s it is accurate between 0 < t < 550 degC in units [mm^2/s]

    This assumes a pressure of 0.5Mpa
    """
    #  Diffusivity 10.1063/1.555718
    TD_ICE = 0.1328     # mm^2/s Diffusivity T=0C P=0.5MPa
    TD_VAP25 = 0.1456   # mm^2/s Diffusivity T=25 P=0.5MPa
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


def getMeltRadius(array, melt_temp):
    """
    Finds the melt radius in cell units given the `array` and melt temperature
    `melt_temp`

    Remember to multiply by dx to get the value in the units of dx.
    """
    length = len(array)
    xc = int(length / 2)  # Centerpoint
    for y in range(0, length - 1):
        if array[xc, y]['Temp'] > melt_temp:
            return xc - y


def do_timestep(u0, u):
    # TODO: This function requires many global constats, these should be passed
    # or determined during the function call.
    # Propagate with forward-difference in time, central-difference in space
    # TODO: Consider using non constant diffusivity
    # NOTE: This does nothing due to the u = u0 line.
    # u[1:-1, 1:-1]['Temp'] = u0[1:-1, 1:-1]['Temp'] + D * dt * (
    #       (
    #         u0[2:, 1:-1]['Temp']
    #         - 2*u0[1:-1, 1:-1]['Temp']
    #         + u0[:-2, 1:-1]['Temp']
    #       )/dx2
    #       + (
    #         u0[1:-1, 2:]['Temp']
    #         - 2*u0[1:-1, 1:-1]['Temp']
    #         + u0[1:-1, :-2]['Temp']
    #       )/dy2
    #      )

    # melt_radius = int(getMeltRadius(u0, Tmelt)*dx)  # Melt Radius
    # rr = r*dx                                  # Root Radius
    # dr = melt_radius - rr
    # rr2 = rr**2
    # u = u0  # This is probably bad practice, but wutevz
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            myT = u0[i, j]['Temp']         # This should minimize array lookups
            myS = u0[i, j]['StateChange']  #
            myE = u0[i, j]['Energy']       #
            # print("{}, {}: T={} S={} E={}".format(i, j, myT, myS, myE))
            if myS:
                # print("{}, {}: T={} S={} E={}".format(i, j, myT, myS, myE))
                # This means the cell is undergoing a state change.
                # Add energy to cell
                # TODO: Find an actual value for this
                Energy_new = EnergyMeltCell/2
                u[i, j]['Energy'] = myE + Energy_new
                # print("State Change in progress at {}, {}".format(i, j))
                if u[i, j]['Energy'] >= EnergyMeltCell:
                    # The cell is melted
                    u[i, j]['StateChange'] = False
                    # print("State Change finished at {}, {}".format(i, j))

            elif myT >= Tmelt:
                # This cell is vapor, thus we can increase it's temperature
                uxx = (u0[i+1, j]['Temp'] - 2*myT + u0[i-1, j]['Temp'])/dx2
                uyy = (u0[i, j+1]['Temp'] - 2*myT + u0[i, j-1]['Temp'])/dy2
                u[i, j]['Temp'] = myT + D * (uxx + uyy)

            else:
                # This cell is ice, we increase the temperature and verify it
                # now isn't above the sublimation temperature.
                uxx = (u0[i+1, j]['Temp'] - 2*myT + u0[i-1, j]['Temp'])/dx2
                uyy = (u0[i, j+1]['Temp'] - 2*myT + u0[i, j-1]['Temp'])/dy2
                u[i, j]['Temp'] = myT + D * (uxx + uyy)
                if u[i, j]['Temp'] >= Tmelt:
                    # We are undergoing a state change now
                    # print("State Change started at {}, {}".format(i, j))
                    u[i, j]['StateChange'] = True

    for i in range(nx):
        for j in range(ny):
            p2 = (i*dx-cx)**2 + (j*dy-cy)**2
            if p2 < r2:
                # u[i+hrr, j]['Temp'] = Thot
                # u[i-hrr, j]['Temp'] = Thot
                # u[i, j+hrr]['Temp'] = Thot
                # u[i, j-hrr]['Temp'] = Thot
                u[i, j]['Temp'] = Thot
    u0 = u.copy()
    return u0, u


# TODO: Separate lib with main
# plate size, mm
w = h = 300.
# intervals in x-, y- directions, mm
dx = dy = 1
# Thermal diffusivity of steel, mm2.s-1
D = 0.1328
pice = 9.168E-7                      # Density of ice [kg/mm^3]
vCell = w**3                         # Volume of a cell [mm^3]
mCell = pice * vCell                 # Mass of Cell of ice [kg]
EIce2Vapor = 2543.0                  # [kJ/kg] at 0C
EnergyMeltCell = EIce2Vapor * mCell  # Energy to melt a cell of ice [kJ]
Tcool, Thot = -20, 500
Tmelt = 10.0

nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = dx2 * dy2 / (2 * D * (dx2 + dy2))
gridArray = np.dtype([('Temp', np.float64),
                      ('StateChange', np.bool),
                      ('Energy', np.float64)])

u0 = np.empty((nx, ny), gridArray)
u = np.empty((nx, ny), gridArray)
u0[:, :]['Temp'] = Tcool
u0[:, :]['StateChange'] = False
u0[:, :]['Energy'] = 0.0


# Initial conditions - ring of inner radius r, width dr centred at (cx,cy) (mm)
r, cx, cy = 16, w/2., h/2.
r2 = r**2
# for i in range(nx):
#     for j in range(ny):
#         p2 = (i*dx-cx)**2 + (j*dy-cy)**2
#         if p2 < r2:
#             u0[i, j] = Thot

hrr = int(r*2*dx)
for i in range(nx):
    for j in range(ny):
        p2 = (i*dx-cx)**2 + (j*dy-cy)**2
        if p2 < r2:
            # u0[i+hrr, j]['Temp'] = Thot
            # u0[i-hrr, j]['Temp'] = Thot
            # u0[i, j+hrr]['Temp'] = Thot
            # u0[i, j-hrr]['Temp'] = Thot
            u0[i, j]['Temp'] = Thot

# TODO: Make this printing a function call
r = getMeltRadius(u0, Tmelt)  # Array, Melt_temp
x = np.arange(0, w, dx)
y = np.arange(0, h, dy)
X, Y = np.meshgrid(x, y)
fig2 = plt.figure()
circle2 = plt.Circle(
    (w/2, h/2), r*dx,
    color='r', linestyle='dashed', fill=False, linewidth=2.0)
ax = plt.gca()
# ax.add_artist(circle2)
cp1 = plt.contourf(X, Y, u0['Temp'].copy(), 20)

plt.title(
    'Contour Plot of Temperature in ice' +
    'after t = {0}{1} and r = {2}{3}'.format(
                int(0), "s",        int(r*dx), "mm"))
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.colorbar(cp1)
plt.savefig("out-{0}.png".format(str(int(0)).zfill(6)))
plt.close()

# Number of timesteps
nsteps = 101
# Output 4 figures at these timesteps
# mfig = [0, int(nsteps/3), int(2*nsteps/3), nsteps - 1]
# fignum = 0

t = 0
while t <= nsteps:

    u0, u = do_timestep(u0, u)
    t += 1
    # print(t)

    # if (t % 101) == 0:
    if True:
        r = getMeltRadius(u, Tmelt)  # Array, Melt_temp
        x = np.arange(0, w, dx)
        y = np.arange(0, h, dy)
        X, Y = np.meshgrid(x, y)
        fig2 = plt.figure()
        circle2 = plt.Circle(
            (w/2, h/2), r*dx,
            color='r', linestyle='dashed', fill=False, linewidth=2.0)
        ax = plt.gca()
        ax.add_artist(circle2)
        cp1 = plt.contourf(X, Y, u['Temp'].copy(), 20)

        plt.title(
            'Contour Plot of Temperature in ice' +
            'after t = {0}{1} and r = {2}{3}'.format(
                        int(t*dt), "s",        int(r*dx), "mm"))
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        plt.colorbar(cp1)
        plt.savefig("out-{0}.png".format(str(int(t*dt)).zfill(6)))
        plt.close()
