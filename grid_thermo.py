#!/usr/bin/python

import pdb
import numpy as np
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
    for i in np.arange(0,int(radius+1),1):
        a = circler(a, i, int(height/2), int(height/2), temp)
    return a

def test_circler():
    SIZE = 100
    a = np.zeros((SIZE, SIZE))
    a = circler(a, 5, int(SIZE/2), int(SIZE/2), 10)
    for i in np.arange(0,6,1):
        a = circler(a, i, int(SIZE/2), int(SIZE/2), 10)

    #plot circle
    delta = 1.0
    x = np.arange(0, 100.0, delta)
    y = np.arange(0, 100.0, delta)
    X, Y = np.meshgrid(x, y)
    fig2 = plt.figure()
    cp1 = plt.contourf(X, Y, a, 2)
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.show()
    # print(a)

# test_circler()




def avg_Temp(V):
    Vrlx=np.zeros((height, length))
    for i in range(height):
        for j in range(length):
            if i == 0 or i == height-1 or j == 0 or j ==length-1:
                Vrlx[i,j]=V[i,j]
            else:
                Vnew = .25*(V[i-1,j]+V[i+1,j]+V[i,j-1]+V[i,j+1])
                Vrlx[i,j]=Vnew
    Vrlx = circle_fill(Vrlx,R,T_FLUID)
    return Vrlx

ANWSER = 42
length = 1000
height = 1000
R = 16.0 # radius (mm)
T_ICE = -40.0 #c
T_FLUID = 500.0 #c


T = np.zeros((height, length))
for i in range(height):
    for j in range(length):
        T[i, j] = T_ICE

T[0, :] = T_ICE
T[-1, :] = T_ICE
T[:, 0] = T_ICE
T[:, -1] = T_ICE
T = circle_fill(T,R,T_FLUID)

t=0
while t<=5000:
    T = avg_Temp(T)
    t+=1
    # print t



    if (t % 10) == 0:
        # print(t, T, sep=",")
        f = "out"
        f += str(t)
        f += ".png"
        x = length/2.0
        y = height/2.0

        delta = 1.0
        x = np.arange(0, length, delta)
        y = np.arange(0, height, delta)
        X, Y = np.meshgrid(x, y)

        fig2 = plt.figure()
        # ax.add_artist(circle1)
        cp1 = plt.contourf(X, Y, T, 25)
        plt.title('Contour Plot of Tempeture in ice')
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        plt.colorbar(cp1)
        # plt.legend()
        plt.savefig(f)






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
# plt.title('Contour Plot of Tempeture in ice')
# plt.xlabel('x (mm)')
# plt.ylabel('y (mm)')
# plt.colorbar(cp1)
# # plt.legend()
# plt.savefig('thermo.png')
