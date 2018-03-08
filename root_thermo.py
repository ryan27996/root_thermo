#!/usr/bin/python
#hello

#second commet
import numpy as np
from matplotlib import pyplot as plt


# Kelvin to Celsius
def ktoc(temp):
    return temp + 273.15


# Celsius to Kelvin
def ctok(temp):
    return temp - 273.15


# NaNo3 NaNo2 KNo3 - Salt 1
# Units: Kg/m^3
def salt1_density(temp):
    return 2293.6 - 0.7497 * temp


# Salt 1 Viscocity
# Units:
# Uncertainty: 16%
# Range: 420-710 K (147-440 C)
def salt1_viscosity(temp):
    return 0.4737 - 2.297*10**(-3) * temp + \
        3.731*10**(-6) * temp**2 - \
        2.019*10**(-6) * temp**3


# Units: J/(kg K)
def salt1_specificHeat(temp):
    return 5806.0 - 10.833 * temp + 7.2413*10**(-3) * temp**2


# Really bad assumption formula
def heat_fluid_flow(density, c, flow_rate, delta_T):
    return density*flow_rate*c*delta_T


def length_of_root(k, Q, r, Delta_T):


# Variable Units
# Temperature: K
# Length:

# Le Assumptions
# Salt Temp Liq Range: 415 to 550 C
SALT_TEMP_MIN = ctok(425.0)
SALT_TEMP_MAX = ctok(550.0)
SALT_AVE_TEMP = (SALT_TEMP_MAX + SALT_TEMP_MIN)/2.0
# Ice Temp: -50 C
ICE_TEMP = ctok(-50.0)
# Fluid tube radius (Single Direction) 1.6 cm
TUBE_RADIUS = 1.6*10**(-2)
# Fluid volume flow rate 150cm^3/s
# TODO:Perhaps this value is pretty wrong, bc was random guess.
VOL_FLOW_RATE = 0.00015
# Conductor surface area based on TUBE_RADIUS
# This is the surface in which the conductive heating flows through.
#CONDUCTOR_SA = 2 * TUBE_RADIUS * TODO:Find Length

temp = np.arange(SALT_TEMP_MIN,SALT_TEMP_MAX + 1, 1)

#heat capacity array
carr = []
density = []
fluid_heat_flow = []
for t in temp:
    carr.append(salt1_specificHeat(t))
    density.append(salt1_density(t))

for i in np.arange(len(temp)):
    fluid_heat_flow.append(heat_fluid_flow(density[i], carr[i], VOL_FLOW_RATE, SALT_AVE_TEMP))


#temp vs Heat capacity
'''
plt.plot(temp,carr, label = 'Heat Capacity')
plt.legend()
plt.show()
'''


#temp vs density
'''
plt.xlabel('Temperature')
plt.ylabel('Density')
plt.plot(temp,density)
plt.show()
'''

#temp vs Fluid heat flow
# plt.xlabel('Temperature')
# plt.ylabel('Fluid Heat Flow')
plt.plot(temp, fluid_heat_flow, label="Heat Flow Rate")
plt.legend()
plt.show()
