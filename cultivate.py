import pdb
import argparse
import numpy as np
# import matplotlib  # For headless servers w/o x-org
# matplotlib.user('Agg')  # For headless servers w/o x-org
import matplotlib.pyplot as plt
import random
from matplotlib.patches import Circle
def avg_Temp(V):

    return Vrlx

# T = np.arange(-20,22,2)
T = [-20,20,20,20,20,20,20]
a = 0.1328E-6

def T_avg(T,alpha):
    for i in range(len(T)):
        if i = 1:
            contiune
        else:
            T[i] = alpha*((T[i] - T[i-1])/2)
    return T

n = 0
while T[-1] > -19:
    print n
    T = T_avg(T,a)
    n+=1


print T








'''
T = np.zeros((50, 50))
T_ICE = -20.0
T[0, :] = T_ICE
T[-1, :] = T_ICE
T[:, 0] = T_ICE
T[:, -1] = T_ICE


while T[50,0] < -10:
    for i in range(50):
        for j in range(50):
            if T[i, j] = -20:
                print("done.")





L = 50.0

t = 100.0


x = np.linspace(0, L, Nx+1)    # mesh points in space
dx = x[1] - x[0]
t = np.linspace(0, T, Nt+1)    # mesh points in time
dt = t[1] - t[0]
F = a*dt/dx**2
u   = zeros(Nx+1)           # unknown u at new time level
u_1 = zeros(Nx+1)           # u at the previous time level

# Set initial condition u(x,0) = I(x)
for i in range(0, Nx+1):
    u_1[i] = I(x[i])

for n in range(0, Nt):
    # Compute u at inner mesh points
    for i in range(1, Nx):
        u[i] = u_1[i] + F*(u_1[i-1] - 2*u_1[i] + u_1[i+1])

    # Insert boundary conditions
    u[0] = 0;  u[Nx] = 0

    # Update u_1 before next step
    u_1[:]= u
    '''
