import numpy as np
import matplotlib.pyplot as plt
import pdb


# waterDiffusivity returns the diffusivity of water at a given temperature
# in mm^2/s it is accurate between 0 < t < 550 degC
def waterDiffusivity(t):

    #  Diffusivity 10.1063/1.555718
    TD_ICE = 0.1328  # mm^2/s Diffusivity T=0C P=0.5MPa
    TD_VAP25 = 0.1456  # mm^2/s Diffusivity T=25 P=0.5MPa
    TD_VAP150 = 0.1725  # mm^2/s Diffusivity T=150 P=0.5MPa
    TD_VAP200 = 6.942  # mm^2/s Diffusivity T=200 P=0.5MPa
    TD_VAP550 = 25.58  # mm^2/s Diffusivity T=550 P=0.5MPa

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
    length = len(array)
    xc = int(length / 2)  # Centerpoint
    for y in range(0, length - 1):
        if array[xc, y]['Temp'] > melt_temp:
            return xc - y


def do_timestep(u0, u):

    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1]['Temp'] = u0[1:-1, 1:-1]['Temp'] + D * dt * (
          (u0[2:, 1:-1]['Temp'] - 2*u0[1:-1, 1:-1]['Temp'] + u0[:-2, 1:-1]['Temp'])/dx2
          + (u0[1:-1, 2:]['Temp'] - 2*u0[1:-1, 1:-1]['Temp'] + u0[1:-1, :-2]['Temp'])/dy2
          )

    # melt_radius = int(getMeltRadius(u0, Tmelt)*dx)  # Melt Radius
    # rr = r*dx                                  # Root Radius
    # dr = melt_radius - rr
    # rr2 = rr**2
    u = u0
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            # Apply State Change
            if u0[i, j]['Temp'] >= Tmelt:
                if u0[i, j]['StateChange'] is True:
                    # add Energy
                    Energy_new = 10
                    # Energy_state = 100
                    Energy_state_change = 30
                    u[i, j]['Energy'] = u0[i, j]['Energy'] + Energy_new    # TODO add energy_new

                    if u[i, j]['Energy'] >= Energy_state_change:
                        u0[i, j]['StateChange'] = False
                        # print("State Change at {}, {}".format(i, j))

                else:
                    uxx = (u0[i+1, j]['Temp'] - 2*u0[i, j]['Temp'] + u0[i-1, j]['Temp']) / dx2
                    uyy = (u0[i, j+1]['Temp'] - 2*u0[i, j]['Temp'] + u0[i, j-1]['Temp']) / dy2
                    u[i, j]['Temp'] = u0[i, j]['Temp'] + D * (uxx + uyy)

            else:
                uxx = (u0[i+1, j]['Temp'] - 2*u0[i, j]['Temp'] + u0[i-1, j]['Temp']) / dx2
                uyy = (u0[i, j+1]['Temp'] - 2*u0[i, j]['Temp'] + u0[i, j-1]['Temp']) / dy2
                u[i, j]['Temp'] = u0[i, j]['Temp'] + dt * D * (uxx + uyy)
                if u[i, j]['Temp'] >= Tmelt:
                    u0[i, j]['StateChange'] = True
                    # u[i, j]['Temp'] = Tmelt

    for i in range(nx):
        for j in range(ny):
            p2 = (i*dx-cx)**2 + (j*dy-cy)**2
            if p2 < r2:
                print("yup")
                # u0[i+hrr, j]['Temp'] = Thot
                # u0[i-hrr, j]['Temp'] = Thot
                # u0[i, j+hrr]['Temp'] = Thot
                # u0[i, j-hrr]['Temp'] = Thot
                u[i, j]['Temp'] = Thot
    u0 = u.copy()
    return u0, u


# plate size, mm
w = h = 500.
# intervals in x-, y- directions, mm
dx = dy = 1
# Thermal diffusivity of steel, mm2.s-1
D = 0.1328

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
u0[:,:]['StateChange'] = False
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


r = getMeltRadius(u0, Tmelt)  # Array, Melt_temp
x = np.arange(0, w, dx)
y = np.arange(0, h, dy)
X, Y = np.meshgrid(x, y)
fig2 = plt.figure()
circle2 = plt.Circle(
    (w/2, h/2), r*dx,
    color='r', linestyle='dashed', fill=False, linewidth=2.0)
ax = plt.gca()
ax.add_artist(circle2)
cp1 = plt.contourf(X, Y, u0['Temp'].copy(), 20)

plt.title(
    'Contour Plot of Tempeture in ice' +
    'after t = {0}{1} and r = {2}{3}'.format(
                int(0), "s",        int(r*dx), "mm"))
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.colorbar(cp1)
plt.savefig("out-{0}.png".format(str(int(0)).zfill(6)))
plt.close()

# Number of timesteps
nsteps = 100001
# Output 4 figures at these timesteps
# mfig = [0, int(nsteps/3), int(2*nsteps/3), nsteps - 1]
# fignum = 0

t = 0
while t <= nsteps:

    u0, u = do_timestep(u0, u)
    t += 1
    # print(t)

    # if (t % 100) == 0:
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
            'Contour Plot of Tempeture in ice' +
            'after t = {0}{1} and r = {2}{3}'.format(
                        int(t*dt), "s",        int(r*dx), "mm"))
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        plt.colorbar(cp1)
        plt.savefig("out-{0}.png".format(str(int(t*dt)).zfill(6)))
        plt.close()
