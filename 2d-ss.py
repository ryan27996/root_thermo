#!/usr/bin/python

import pdb
import argparse
import numpy as np
# import matplotlib  # For headless servers w/o x-org
# matplotlib.user('Agg')  # For headless servers w/o x-org
import matplotlib.pyplot as plt
import random
from matplotlib.patches import Circle


# Kelvin to Celsius
def ktoc(temp):
    return temp - 273.15


# Celsius to Kelvin
def ctok(temp):
    return temp + 273.15


def circler(array, radius, center_x, center_y, temperature):
    a = array
    r = radius
    x_c = center_x
    y_c = center_y
    t = temperature

    f = 1 - r
    ddf_x = 1
    ddf_y = -2 * r
    x = 0
    y = r
    a[x_c, y_c - r] = t
    a[x_c, y_c + r] = t
    a[x_c + r, y_c] = t
    a[x_c - r, y_c] = t

    while x < y:
        if f >= 0:
            y -= 1
            ddf_y += 2
            f += ddf_y
        x += 1
        ddf_x += 2
        f += ddf_x
        a[x_c + x, y_c + y] = t
        a[x_c - x, y_c + y] = t
        a[x_c + x, y_c - y] = t
        a[x_c - x, y_c - y] = t
        a[x_c + y, y_c + x] = t
        a[x_c - y, y_c + x] = t
        a[x_c + y, y_c - x] = t
        a[x_c - y, y_c - x] = t

    return a


def circle_fill(a, radius, temp):
    for i in np.arange(0, int(radius+1), 1):
        a = circler(a, i, int(height/2), int(height/2), temp)
    return a


def test_circler():
    SIZE = 100
    a = np.zeros((SIZE, SIZE))
    a = circler(a, 5, int(SIZE/2), int(SIZE/2), 10)
    for i in np.arange(0, 6, 1):
        a = circler(a, i, int(SIZE/2), int(SIZE/2), 10)

    # plot circle
    delta = 1.0
    x = np.arange(0, 100.0, delta)
    y = np.arange(0, 100.0, delta)
    X, Y = np.meshgrid(x, y)
    plt.figure()
    plt.contourf(X, Y, a, 2)
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.show()
    # print(a)

# test_circler()


def lim_cir(r, center_x, center_y):
    circle2 = plt.Circle((center_x, center_y), r, color='r', fill=False, linewidth = 2.0, linestyle='dashed')
    ax = plt.gca()
    ax.add_artist(circle2)


def avg_Temp(V):
    Vrlx = np.zeros((height, length))
    for i in range(height):
        for j in range(length):
            if i == 0 or i == height-1 or j == 0 or j == length-1:
                Vrlx[i, j] = V[i, j]
            else:
                Vnew = .25*(Vrlx[i-1, j]+V[i+1, j]+Vrlx[i, j-1]+V[i, j+1])
                Vrlx[i, j] = Vnew
    Vrlx = circle_fill(Vrlx, R, T_FLUID)
    return Vrlx


def getMeltRadius(array, length, melt_temp):
    x = int(length / 2)
    # r1 = 0
    r2 = 0
    # for y in range(x, 0, -1):
    #     if array[x, y] <= melt_temp:
    #         r1 = y
    #         break
    for y in range(0, length - 1):
        if array[x, y] > melt_temp:
            r2 = x - y
            break
    # print(r1, r2)
    return r2


#  Example
#   ./2d-ss.py  -l 500 -m 10 -i -40 -c 500 -r 1.6 -n 1000

parser = argparse.ArgumentParser()
parser.add_argument("--length", "-l", type=int, required=True)
parser.add_argument("--melttemp", "-m", type=float, required=True)
parser.add_argument("--icetemp", "-i", type=float, required=True)
parser.add_argument("--coretemp", "-c", type=float, required=True)
parser.add_argument("--radius", "-r", type=float, required=True)
parser.add_argument("--iterations", "-n", type=int, required=True)
args = parser.parse_args()


# waterDiffusivity returns the diffusivity of water at a given temperature
# it is accurate between 0 < t < 550 degC
def waterDiffusivity(t):

    #  Diffusivity 10.1063/1.555718
    TD_ICE = 0.1328E-6  # m^2/s Diffusivity T=0C P=0.5MPa
    TD_VAP25 = 0.1456E-6  # m^2/s Diffusivity T=25 P=0.5MPa
    TD_VAP150 = 0.1725E-6  # m^2/s Diffusivity T=150 P=0.5MPa
    TD_VAP200 = 6.942E-6  # m^2/s Diffusivity T=200 P=0.5MPa
    TD_VAP550 = 25.58E-6  # m^2/s Diffusivity T=550 P=0.5MPa

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


ANWSER = 42
#  length = 100
length = args.length
height = length
#  R = 1.6  # radius (cm)
R = args.radius
# T_ICE = -40.0  # c
T_ICE = args.icetemp
# T_FLUID = 500.0  # c
T_FLUID = args.coretemp
T_MELT = args.melttemp
ITERATIONS = args.iterations


T = np.zeros((height, length))
for i in range(height):
    for j in range(length):
        T[i, j] = T_ICE

T[0, :] = T_ICE
T[-1, :] = T_ICE
T[:, 0] = T_ICE
T[:, -1] = T_ICE
T = circle_fill(T, R, T_FLUID)

t = 0
while t <= ITERATIONS:
    T = avg_Temp(T)
    t += 1
    # print t

    if (t % 100) == 0:
        # print(t, T, sep=",")
        r = getMeltRadius(T, length, T_MELT)  # Array, Length, Melt_temp
        f = "out-"
        f += str(t).zfill(6)
        f += "-"
        f += str(r)
        f += ".png"
        x = length/2.0
        y = height/2.0

        lim_cir(r, x, y)
        delta = 1.0
        x = np.arange(0, length, delta)
        y = np.arange(0, height, delta)
        X, Y = np.meshgrid(x, y)
        fig2 = plt.figure()
        circle2 = plt.Circle((length/2.0, height/2.0), r, color='r', fill=False, linewidth = 2.0, linestyle='dashed')
        ax = plt.gca()
        ax.add_artist(circle2)
        cp1 = plt.contourf(X, Y, T, 20)
        plt.title('Contour Plot of Temperature in ice')
        plt.xlabel('x (cm)')
        plt.ylabel('y (cm)')
        plt.colorbar(cp1)
        # plt.legend()
        plt.savefig(f)
        plt.close()


# x = length/2.0
# y = height/2.0
#
# '''
#
# circle1 = plt.Circle((x, y), 16.0, color='r')
#
# '''
#
#
#
#
#
#
# # T = np.flipud(T)
#
# delta = 1.0
# x = np.arange(0, length, delta)
# y = np.arange(0, height, delta)
# X, Y = np.meshgrid(x, y)
#
# fig2 = plt.figure()
# # ax.add_artist(circle1)
# cp1 = plt.contourf(X, Y, T, 25)
# plt.title('Contour Plot of Temperature in ice')
# plt.xlabel('x (mm)')
# plt.ylabel('y (mm)')
# plt.colorbar(cp1)
# # plt.legend()
# plt.savefig('thermo.png')
