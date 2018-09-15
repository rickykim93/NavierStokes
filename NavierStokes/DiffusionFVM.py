# -*- coding: utf-8 -*-
"""
FILE DiffusionFVM.py
Solves Diffusion Equation using Finite Volume Method

Name: Kyu Mok Kim
Student Number: 998745381

MIE1210 Project 2 #1

@author: kimkyu4

This code is designed to solve the diffusion equation,
using FVM method with iterative TDMA.
Refer to my report for detailed explanation.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

np.set_printoptions(threshold=np.inf)

#TODO: Finish this part later

if __name__ == '__main__':

    # Set the grid dimensions
    lx = 1.
    ly = 1.
    # Set the grid resolution (number of nodes including boundaries)
    nx = 5
    ny = 5
    # Compute the grid spacing
    dx = lx / (nx)
    dy = ly / (ny)
    idx = np.arange(0, nx * ny).reshape((nx, ny))
    B = np.zeros((nx, ny)).flatten()
    # Create lists to hold the non-zero entries in coo format
    data = np.ones(2 * nx + 2 * (ny - 2)).tolist()
    rows = np.concatenate((idx[:, [0, -1]].flatten(), idx[[0, -1], 1:-1].flatten())).tolist()
    cols = list(rows)
    A = sp.coo_matrix((data, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()
    # The central part only have diffusion term.
    a_c = 2 * 20 / dx + 2 * 20 / dy
    a_e = -20 / dx
    a_w = -20 / dx
    a_n = -20 / dy
    a_s = -20 / dy

    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            data.extend([a_c, a_e, a_w, a_n, a_s])
            rows.extend(np.ones(5) * [idx[i, j]])
            cols.extend([idx[i, j], idx[i, j + 1], idx[i, j - 1], idx[i - 1, j], idx[i + 1, j]])
    # The west boundary part, Tw=10 degree
    Tw = 10
    a_c = 3 * 20 / dx + 2 * 20 / dy
    a_e = -20 / dx
    a_n = -20 / dy
    a_s = -20 / dy

    for i in range(1, nx - 1):
        data.extend([a_c, a_e, a_n, a_s])
        rows.extend(np.ones(4) * [idx[i, 0]])
        cols.extend([idx[i, 0], idx[i, 0 + 1], idx[i - 1, 0], idx[i + 1, 0]])
        B[idx[i, 0]] = 2 * 20 * Tw / dx  # aw boundary Tw=10
    # The east boundary part, Te=100 degree
    Te = 100
    a_c = 3 * 20 / dx + 2 * 20 / dy
    a_w = -20 / dx
    a_n = -20 / dy
    a_s = -20 / dy

    for i in range(1, nx - 1):
        data.extend([a_c, a_w, a_n, a_s])
        rows.extend(np.ones(4) * [idx[i, -1]])
        cols.extend([idx[i, -1], idx[i, -2], idx[i - 1, -1], idx[i + 1, -1]])

        B[idx[i, -1]] = 2 * 20 * Te / dx  # ae boundary Te=100
    # The south boundary part, dT/dy=0 as=0
    a_c = 2 * 20 / dx + 20 / dy
    a_w = -20 / dx
    a_e = -20 / dx
    a_n = -20 / dy

    for j in range(1, ny - 1):
        i = nx - 1
        data.extend([a_c, a_w, a_e, a_n])
        rows.extend(np.ones(4) * [idx[i, j]])
        cols.extend([idx[i, j], idx[i, j - 1], idx[i, j + 1], idx[i - 1, j]])
    # The southwestern corner, dT/dy=0 as=0,Tw=10 degree
    Tw = 10
    a_c = 3 * 20 / dx + 20 / dy
    a_e = -20 / dx
    a_n = -20 / dy

    data.extend([a_c, a_e, a_n])
    rows.extend(np.ones(3) * [idx[nx - 1, 0]])
    cols.extend([idx[nx - 1, 0], idx[nx - 1, 0 + 1], idx[nx - 2, 0]])

    B[idx[nx - 1, 0]] = 2 * 20 * Tw / dx  # aw
    A = sp.coo_matrix((data, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()
    # The southeastern corner, dT/dy=0 as=0,Te=100 degree
    Te = 100
    a_c = 3 * 20 / dx + 20 / dy
    a_w = -20 / dx
    a_n = -20 / dy

    data.extend([a_c, a_w, a_n])
    rows.extend(np.ones(3) * [idx[nx - 1, -1]])
    cols.extend([idx[nx - 1, -1], idx[nx - 1, -2], idx[nx - 2, -1]])

    B[idx[nx - 1, -1]] = 2 * 20 * Te / dx  # ae
    A = sp.coo_matrix((data, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()
    # k*dT/dy=h(Tinf-T) >>> k*(Tn-Tp)/(dy/2)=h(T_ext-Tn)>>> Tn=[h*dy/(2k+h*dy)]*Tinf+[2k/(2k+h*dy)]*Tp
    T_ext = 300
    h = 10
    k = 20
    Sp = (h) / (1 + h * dy / 2 / k)

    Su = (h) / (1 + h * dy / 2 / k) * T_ext
    # The north boundary part
    a_c = 2 * 20 / dx + 20 / dy - Sp
    a_w = -20 / dx
    a_e = -20 / dx
    a_s = -20 / dy

    for j in range(1, ny - 1):
        i = 0
        data.extend([a_c, a_w, a_e, a_s])
        rows.extend(np.ones(4) * [idx[i, j]])
        cols.extend([idx[i, j], idx[i, j - 1], idx[i, j + 1], idx[i + 1, j]])

        B[idx[i, j]] = Su
    # The northwestern corner, Tw=10 degree
    Tw = 10
    a_c = 3 * 20 / dx + 20 / dy - Sp
    a_e = -20 / dx
    a_s = -20 / dy

    data.extend([a_c, a_e, a_s])
    rows.extend(np.ones(3) * [idx[0, 0]])
    cols.extend([idx[0, 0], idx[0, 1], idx[1, 0]])

    B[idx[0, 0]] = 2 * 20 * Tw / dx + Su  # aw & an
    # The northeastern corner, Te=100 degree
    Te = 100
    a_c = 3 * 20 / dx + 20 / dy - Sp
    a_w = -20 / dx
    a_s = -20 / dy

    data.extend([a_c, a_w, a_s])
    rows.extend(np.ones(3) * [idx[0, -1]])
    cols.extend([idx[0, -1], idx[0, -2], idx[1, -1]])

    B[idx[0, ny - 1]] = 2 * 20 * Te / dx + Su  # ae & an
    A = sp.coo_matrix((data, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()
    T = spla.spsolve(A, B)
    T.reshape(nx, ny)
    residual_80 = A * T - B
    np.savetxt("T_Coo_80.csv", T, delimiter=",")

    import math

    r2 = residual_80 ** 2
    r2 = math.sqrt(np.sum(r2))
    print (r2)
    T = T.reshape(nx, ny)
    print (T)

    # In[105]:

    T = T.reshape(nx, ny)
    xlist = np.linspace(0.0, 1.0, nx)
    ylist = np.linspace(0.0, 1.0, ny)
    X, Y = np.meshgrid(xlist, ylist)
    plt.figure()
    cp = plt.contourf(X, Y, T)
    plt.colorbar(cp)
    plt.title('80*80 mesh Contours Plot')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()
