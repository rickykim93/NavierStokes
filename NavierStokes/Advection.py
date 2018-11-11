# -*- coding: utf-8 -*-
"""
FILE Advection-DiffusionFVM.py
Solves Advection-Diffusion Equation using Finite Volume Method

Name: Kyu Mok Kim
Student Number: 998745381

MIE1210 Project 3 #2 and 3

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
    ux = 2  # velocity component x
    uy = 2  # velocity compoenent y
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
        a_e = (De + max(0, -Fe))
        a_w = (Dw + max(Fw, 0))
        a_n = (Dn + max(0, -Fn))
        a_s = (Ds + max(Fs, 0))

    # Middle Points (Within Boundaries)
    a_c = a_w + a_e + (Fe - Fw) + a_n + a_s + (Fn - Fs)
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            data.extend([a_c, -a_e, -a_w, -a_n, -a_s])
            rows.extend([idx[i, j]] * 5)
            cols.extend([idx[i, j], idx[i, j + 1], idx[i, j - 1], idx[i - 1, j], idx[i + 1, j]])

    # West Boundary condition
    a_c = a_e + (Fe - Fw) + a_n + a_s + (Fn - Fs) + (2 * Dw + Fw)
    for i in range(1, nx - 1):
        data.extend([a_c, -a_e, -a_n, -a_s])
        rows.extend([idx[i, 0]] * 4)
        cols.extend([idx[i, 0], idx[i, 1], idx[i - 1, 0], idx[i + 1, 0]])
        d[idx[i, 0]] = (2 * Dw + Fw) * Tw

    # North Boundary condition
    if method == 1:
        a_c = a_w + a_e + (Fe - Fw) + a_s + (Fn - Fs) + (2 * Dn - Fn)
    if method == 2:
        a_c = a_w + a_e + (Fe - Fw) + a_s + (Fn - Fs) + 2 * Dn
    for j in range(1, ny - 1):
        data.extend([a_c, -a_w, -a_e, -a_s])
        rows.extend([idx[0, j]] * 4)
        cols.extend([idx[0, j], idx[0, j - 1], idx[0, j + 1], idx[1, j]])
        if method == 1:
            d[idx[0, j]] = (2 * Dn - Fn) * Tn
        if method == 2:
            d[idx[0, j]] = 2 * Dn * Tn
    # East Boundary condition
    if method == 1:
        a_c = a_w + (Fe - Fw) + a_n + a_s + (Fn - Fs) + (2 * De - Fe)
    if method == 2:
        a_c = a_w + (Fe - Fw) + a_n + a_s + (Fn - Fs) + 2 * De
    for i in range(1, nx - 1):
        data.extend([a_c, -a_w, -a_n, -a_s])
        rows.extend([idx[i, -1]] * 4)
        cols.extend([idx[i, -1], idx[i, -2], idx[i - 1, -1], idx[i + 1, -1]])
        if method == 1:
            d[idx[i, -1]] = (2 * De - Fe) * Te
        if method == 2:
            d[idx[i, -1]] = 2 * De * Te

    # South Boundary condition
    a_c = a_w + a_e + (Fe - Fw) + a_n + (Fn - Fs) + (2 * Ds + Fs)
    for j in range(1, ny - 1):
        data.extend([a_c, -a_w, -a_e, -a_n])
        rows.extend([idx[nx - 1, j]] * 4)
        cols.extend([idx[nx - 1, j], idx[nx - 1, j - 1], idx[nx - 1, j + 1], idx[nx - 2, j]])
        d[idx[-1, j]] = (2 * Ds + Fs) * Ts

    # South West Corner
    a_c = a_e + (Fe - Fw) + a_n + (Fn - Fs) + (2 * Ds + Fs) + (2 * Dw + Fw)
    data.extend([a_c, -a_e, -a_n])
    rows.extend([idx[nx - 1, 0]] * 3)
    cols.extend([idx[nx - 1, 0], idx[nx - 1, 1], idx[nx - 2, 0]])
    d[idx[nx - 1, 0]] = (2 * Dw + Fw) * Tw

    # South East Corner
    a_c = a_w + (Fe - Fw) + a_n + (Fn - Fs) + (2 * Ds + Fs) + (2 * De - Fe)
    data.extend([a_c, -a_w, -a_n])
    rows.extend([idx[nx - 1, -1]] * 3)
    cols.extend([idx[nx - 1, -1], idx[nx - 1, -2], idx[nx - 2, -1]])

    # North East Corner
    if method == 1:
        a_c = a_w + (Fe - Fw) + a_s + (Fn - Fs) + (2 * Dn - Fn) + (2 * De - Fe)
    if method == 2:
        a_c = a_w + (Fe - Fw) + a_s + (Fn - Fs) + 2 * Dn + 2 * De
    data.extend([a_c, -a_w, -a_s])
    rows.extend([idx[0, -1]] * 3)
    cols.extend([idx[0, -1], idx[0, -2], idx[1, -1]])
    if method == 1:
        d[idx[0, ny - 1]] = (2 * Dn - Fn) * Tn
    if method == 2:
        d[idx[0, ny - 1]] = 2 * Dn * Tn

    # North West Corner
    if method == 1:
        a_c = a_e + (Fe - Fw) + a_s + (Fn - Fs) + (2 * Dn - Fn) + (2 * Dw + Fw)
    if method == 2:
        a_c = a_e + (Fe - Fw) + a_s + (Fn - Fs) + 2 * Dn + (2 * Dw + Fw)
    data.extend([a_c, -a_e, -a_s])
    rows.extend([idx[0, 0]] * 3)
    cols.extend([idx[0, 0], idx[0, 1], idx[1, 0]])
    if method == 1:
        d[idx[0, 0]] = (2 * Dw + Fw) * Tw + (2 * Dn - Fn) * Tn
    if method == 2:
        d[idx[0, 0]] = (2 * Dw + Fw) * Tw + 2 * Dn * Tn

    # Construct coo sparse matrix, and convert it to csr format for better efficiency
    A = sp.coo_matrix((data, (rows, cols)), shape=(nx * ny, nx * ny)).tocsr()

    # Solve the system. Pass in the T values as the rhs vector (zero everywhere except boundaries)
    T = spla.spsolve(A, d)
    T = T.reshape(nx, ny)
    return T


def main():
    nx = ny = 320
    method = 2  # 1 for CD, 2 for Upwind
    T = (AdvecSolver(method, nx, ny))
    if method == 1:
        m = "CD"
    if method == 2:
        m = "Upwind"

    import matplotlib.pyplot as plt

    # For comparing with exact solution
    E = np.zeros(nx)
    for x in range(1, nx):
        E[x - 1] = T[x - 1, x - 1]
    y = np.linspace(0, 1, nx)
    plt.plot(y, E)
    # Exact Solution
    F = np.zeros(nx)
    g = 5  # gamma
    for x in range(1, nx):
        if g == 0:
            if x < (nx + 1) / 2:
                F[x - 1] = 100
        if g != 0:
            F[x - 1] = 100 * (np.exp(1.2 * 2 * (1 - x / nx) / g) - 1) / (np.exp(1.2 * 2 / g) - 1)
    plt.plot(y, F)
    plt.axis([0, 1, 0, 120])
    plt.title('%d by %d %s Graph' % (nx, ny, m))
    plt.xlabel('x (m)')
    plt.ylabel('Temperature (Celcius)')
    plt.show()

    # Create an x, y grid and graph
    cp = plt.imshow(T)
    plt.colorbar(cp)
    plt.title('%d by %d Contour Plot for %s' % (nx, ny, m))
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()


if __name__ == '__main__':
    import time

    start = time.time()
    main()
    end = time.time()
    print('Execution time =', end - start)
