# -*- coding: utf-8 -*-
"""
FILE Advection-DiffusionFVM.py
Solves Advection-Diffusion Equation using Finite Volume Method

Name: Kyu Mok Kim
Student Number: 998745381

MIE1210 Project 3 #4

@author: kimkyu4

This code is designed to solve the advection-diffusion equation,
using FVM method with Scipy.
Refer to my report for detailed explanation.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

np.set_printoptions(threshold=np.inf)


# Method: 1 for central difference, 2 for upwind
def AdvecSolver(method, nx, ny):
    # Set the grid dimensions
    lx = 1.
    ly = 1.

    # Compute the grid spacing
    dx = lx / nx
    dy = ly / ny

    # Initialize the temperature array to zeros
    d = np.zeros((nx, ny)).flatten()

    # Create the index map, i.e. assign grid nodes to matrix rows
    idx = np.arange(0, nx * ny).reshape((nx, ny))

    # Create lists to hold the non-zero entries in coo format
    data = np.ones(2 * nx + 2 * (ny - 2)).tolist()
    rows = np.concatenate((idx[:, [0, -1]].flatten(), idx[[0, -1], 1:-1].flatten())).tolist()
    cols = list(rows)

    # Parameters
    g = 5  # gamma
    ux = np.zeros(shape=(nx, ny))
    uy = np.zeros(shape=(nx, ny))
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            y = -((i - 1) - (nx - 1) / 2)
            x = (j - 1) - (ny - 1) / 2
            r = np.sqrt(np.power(y, 2) + np.power(x, 2))
            if ((j - 1) - (ny - 1) / 2) > 0:
                theta = np.arctan(y / x)
            else:
                theta = np.pi + np.arctan(y / x)
            ux[i - 1, j - 1] = -r * np.sin(theta)
            uy[i - 1, j - 1] = r * np.cos(theta)
    rho = 1.2  # density of air (kg/m^3)
    De = Dw = g / dx
    Dn = Ds = g / dy
    Fe = Fw = rho * ux
    Fn = Fs = rho * uy
    # Pe=Fe/De

    # Boundary Conditions
    Tw = Tn = 100
    Ts = Te = 0

    if method == 1:  # Central Differencing
        a_e = (De - Fe / 2)
        a_w = (Dw + Fw / 2)
        a_n = (Dn - Fn / 2)
        a_s = (Ds + Fs / 2)

    if method == 2:  # Upwind Differencing
        a_e = np.zeros(shape=(nx, ny))
        a_w = np.zeros(shape=(nx, ny))
        a_n = np.zeros(shape=(nx, ny))
        a_s = np.zeros(shape=(nx, ny))
        for j in range(1, ny):
            for i in range(1, nx):
                a_e[i - 1, j - 1] = (De + max(0, -Fe[i - 1, j - 1]))
                a_w[i - 1, j - 1] = (Dw + max(Fw[i - 1, j - 1], 0))
                a_n[i - 1, j - 1] = (Dn + max(0, -Fn[i - 1, j - 1]))
                a_s[i - 1, j - 1] = (Ds + max(Fs[i - 1, j - 1], 0))

    # Middle Points (Within Boundaries)
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            a_c = a_w[i, j] + a_e[i, j] + (Fe[i, j] - Fw[i, j]) + a_n[i, j] + a_s[i, j] + (Fn[i, j] - Fs[i, j])
            data.extend([a_c, -a_e[i, j], -a_w[i, j], -a_n[i, j], -a_s[i, j]])
            rows.extend([idx[i, j]] * 5)
            cols.extend([idx[i, j], idx[i, j + 1], idx[i, j - 1], idx[i - 1, j], idx[i + 1, j]])

    # West Boundary condition
    for i in range(1, nx - 1):
        a_c = a_e[i - 1, 0] + (Fe[i - 1, 0] - Fw[i - 1, 0]) + a_n[i - 1, 0] + a_s[i - 1, 0] + (
                    Fn[i - 1, 0] - Fs[i - 1, 0]) + (2 * Dw + Fw[i - 1, 0])
        data.extend([a_c, -a_e[i - 1, 0], -a_n[i - 1, 0], -a_s[i - 1, 0]])
        rows.extend([idx[i, 0]] * 4)
        cols.extend([idx[i, 0], idx[i, 1], idx[i - 1, 0], idx[i + 1, 0]])
        d[idx[i, 0]] = (2 * Dw + Fw[i - 1, 0]) * Tw

    # North Boundary condition
    for j in range(1, ny - 1):
        if method == 1:
            a_c = a_w[0, j - 1] + a_e[0, j - 1] + (Fe[0, j - 1] - Fw[0, j - 1]) + a_s[0, j - 1] + (
                        Fn[0, j - 1] - Fs[0, j - 1]) + (2 * Dn - Fn[0, j - 1])
        if method == 2:
            a_c = a_w[0, j - 1] + a_e[0, j - 1] + (Fe[0, j - 1] - Fw[0, j - 1]) + a_s[0, j - 1] + (
                        Fn[0, j - 1] - Fs[0, j - 1]) + 2 * Dn
        data.extend([a_c, -a_w[0, j - 1], -a_e[0, j - 1], -a_s[0, j - 1]])
        rows.extend([idx[0, j]] * 4)
        cols.extend([idx[0, j], idx[0, j - 1], idx[0, j + 1], idx[1, j]])
        if method == 1:
            d[idx[0, j]] = (2 * Dn - Fn[0, j - 1]) * Tn
        if method == 2:
            d[idx[0, j]] = 2 * Dn * Tn
    # East Boundary condition
    for i in range(1, nx - 1):
        if method == 1:
            a_c = a_w[i - 1, -1] + (Fe[i - 1, -1] - Fw[i - 1, -1]) + a_n[i - 1, -1] + a_s[i - 1, -1] + (
                        Fn[i - 1, -1] - Fs[i - 1, -1]) + (2 * De - Fe[i - 1, -1])
        if method == 2:
            a_c = a_w[i - 1, -1] + (Fe[i - 1, -1] - Fw[i - 1, -1]) + a_n[i - 1, -1] + a_s[i - 1, -1] + (
                        Fn[i - 1, -1] - Fs[i - 1, -1]) + 2 * De
        data.extend([a_c, -a_w[i - 1, -1], -a_n[i - 1, -1], -a_s[i - 1, -1]])
        rows.extend([idx[i, -1]] * 4)
        cols.extend([idx[i, -1], idx[i, -2], idx[i - 1, -1], idx[i + 1, -1]])
        if method == 1:
            d[idx[i, -1]] = (2 * De - Fe[i - 1, -1]) * Te
        if method == 2:
            d[idx[i, -1]] = 2 * De * Te

    # South Boundary condition
    for j in range(1, ny - 1):
        a_c = a_w[nx - 1, j - 1] + a_e[nx - 1, j - 1] + (Fe[nx - 1, j - 1] - Fw[nx - 1, j - 1]) + a_n[nx - 1, j - 1] + (
                    Fn[nx - 1, j - 1] - Fs[nx - 1, j - 1]) + (2 * Ds + Fs[nx - 1, j - 1])
        data.extend([a_c, -a_w[nx - 1, j - 1], -a_e[nx - 1, j - 1], -a_n[nx - 1, j - 1]])
        rows.extend([idx[nx - 1, j]] * 4)
        cols.extend([idx[nx - 1, j], idx[nx - 1, j - 1], idx[nx - 1, j + 1], idx[nx - 2, j]])
        d[idx[-1, j]] = (2 * Ds + Fs[nx - 1, j - 1]) * Ts

    # South West Corner
    a_c = a_e[nx - 1, 0] + (Fe[nx - 1, 0] - Fw[nx - 1, 0]) + a_n[nx - 1, 0] + (Fn[nx - 1, 0] - Fs[nx - 1, 0]) + (
                2 * Ds + Fs[nx - 1, 0]) + (2 * Dw + Fw[nx - 1, 0])
    data.extend([a_c, -a_e[nx - 1, 0], -a_n[nx - 1, 0]])
    rows.extend([idx[nx - 1, 0]] * 3)
    cols.extend([idx[nx - 1, 0], idx[nx - 1, 1], idx[nx - 2, 0]])
    d[idx[nx - 1, 0]] = (2 * Dw + Fw[nx - 1, 0]) * Tw

    # South East Corner
    a_c = a_w[nx - 1, -1] + (Fe[nx - 1, -1] - Fw[nx - 1, -1]) + a_n[nx - 1, -1] + (Fn[nx - 1, -1] - Fs[nx - 1, -1]) + (
                2 * Ds + Fs[nx - 1, -1]) + (2 * De - Fe[nx - 1, -1])
    data.extend([a_c, -a_w[nx - 1, -1], -a_n[nx - 1, -1]])
    rows.extend([idx[nx - 1, -1]] * 3)
    cols.extend([idx[nx - 1, -1], idx[nx - 1, -2], idx[nx - 2, -1]])

    # North East Corner
    if method == 1:
        a_c = a_w[0, -1] + (Fe[0, -1] - Fw[0, -1]) + a_s[0, -1] + (Fn[0, -1] - Fs[0, -1]) + (2 * Dn - Fn[0, -1]) + (
                    2 * De - Fe[0, -1])
    if method == 2:
        a_c = a_w[0, -1] + (Fe[0, -1] - Fw[0, -1]) + a_s[0, -1] + (Fn[0, -1] - Fs[0, -1]) + 2 * Dn + 2 * De
    data.extend([a_c, -a_w[0, -1], -a_s[0, -1]])
    rows.extend([idx[0, -1]] * 3)
    cols.extend([idx[0, -1], idx[0, -2], idx[1, -1]])
    if method == 1:
        d[idx[0, ny - 1]] = (2 * Dn - Fn[0, -1]) * Tn
    if method == 2:
        d[idx[0, ny - 1]] = 2 * Dn * Tn

    # North West Corner
    if method == 1:
        a_c = a_e[0, 0] + (Fe[0, 0] - Fw[0, 0]) + a_s[0, 0] + (Fn[0, 0] - Fs[0, 0]) + (2 * Dn - Fn[0, 0]) + (
                    2 * Dw + Fw[0, 0])
    if method == 2:
        a_c = a_e[0, 0] + (Fe[0, 0] - Fw[0, 0]) + a_s[0, 0] + (Fn[0, 0] - Fs[0, 0]) + 2 * Dn + (2 * Dw + Fw[0, 0])
    data.extend([a_c, -a_e[0, 0], -a_s[0, 0]])
    rows.extend([idx[0, 0]] * 3)
    cols.extend([idx[0, 0], idx[0, 1], idx[1, 0]])
    if method == 1:
        d[idx[0, 0]] = (2 * Dw + Fw[0, 0]) * Tw + (2 * Dn - Fn[0, 0]) * Tn
    if method == 2:
        d[idx[0, 0]] = (2 * Dw + Fw[0, 0]) * Tw + 2 * Dn * Tn

    # Construct coo sparse matrix, and convert it to csr format for better efficiency
    A = sp.coo_matrix((data, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()

    # Solve the system. Pass in the T values as the rhs vector (zero everywhere except boundaries)
    T = spla.spsolve(A, d)
    T = T.reshape(nx, ny)
    return T


# Helper function to reshape the large matrix to smaller by taking average for easier error calculation
def matavg(a, shape):
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def FindConv(method):
    Tff = AdvecSolver(method, 320, 320)
    Tf = AdvecSolver(method, 160, 160)
    Tc = AdvecSolver(method, 80, 80)
    Tfff = matavg(Tff, (160, 160))
    Tffc = matavg(Tfff, (80, 80))
    ect = 0
    eft = 0
    for c in range(1, 81):
        for d in range(1, 81):
            ect = ect + (Tffc[c - 1, d - 1] - Tc[c - 1, d - 1]) ** 2
    ec = (1 / 80 / 80 * ect) ** 0.5

    for e in range(1, 161):
        for f in range(1, 161):
            eft = eft + (Tfff[e - 1, f - 1] - Tf[e - 1, f - 1]) ** 2
    ef = (1 / 160 / 160 * eft) ** 0.5

    conv = np.log(ec / ef) / np.log(160 / 80)
    return conv


def main():
    nx = ny = 320
    method = 2  # 1 for CD, 2 for Upwind
    if method == 1:
        m = "CD"
    if method == 2:
        m = "Upwind"
    T = (AdvecSolver(method, nx, ny))

    import matplotlib.pyplot as plt
    # Create an x, y grid and graph
    cp = plt.imshow(T)
    plt.colorbar(cp)
    plt.title('%d by %d Contour Plot for %s' % (nx, ny, m))
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

    conv1 = FindConv(1)
    conv2 = FindConv(2)
    print('Order of Convergence for CD is:', conv1)
    print('Order of Convergence for Upwind is:', conv2)


if __name__ == '__main__':
    import time

    start = time.time()
    main()
    end = time.time()
    print('Execution time =', end - start)