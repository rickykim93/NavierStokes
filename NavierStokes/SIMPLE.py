# -*- coding: utf-8 -*-
"""
FILE SIMPLE lid driven cavity.py
Solves lid driven cavity problem using SIMPLE

Name: Kyu Mok Kim
Student Number: 998745381

MIE1210 Project 4

@author: kimkyu4

This code is designed to solve Navier Stokes,
using FVM method with SIMPLE algorithm.
Refer to my report for detailed explanation.
"""

#TODO: convert below to sparse matrix

import numpy as np

np.set_printoptions(threshold=np.inf)


def LidCavity(nx, ny, Re):
    # variables
    Re = Re  # Reynolds
    nx = nx  # mesh size
    ny = ny
    lx = 1  # cavity dimensions (m)
    ly = 1  # for caivity
    # ly=1/4  #for backstep flow
    ulid = 1  # lid velocity (m/s)
    rho = 1  # density (kg/m^3)
    nu = 1  # viscosity (m^2/s)
    rel = 0.5  # relaxation factor for velocities
    relp = 0.8  # relaxation factor for pressure
    dx = lx / nx  # grid spacing
    dy = ly / ny

    # initial conditions
    u = np.zeros((nx + 1, ny + 2))  # u (velocity in x direction)
    us = np.zeros((nx + 1, ny + 2))  # ustar: the guesses
    up = np.zeros((nx + 1, ny + 2))  # uprime: correction term
    v = np.zeros((nx + 2, ny + 1))  # v (velocity in y direction)
    vs = np.zeros((nx + 2, ny + 1))
    vp = np.zeros((nx + 2, ny + 1))
    p = np.zeros((nx + 2, ny + 2))  # pressure terms
    ps = np.zeros((nx + 2, ny + 2))
    pp = np.zeros((nx + 2, ny + 2))
    du = np.zeros((nx + 1, ny + 2))  # d terms: d=A/ap
    dv = np.zeros((nx + 2, ny + 1))
    uOld = np.zeros((nx + 1, ny + 2))
    vOld = np.zeros((nx + 2, ny + 1))
    # coefficients
    uae = np.zeros((nx + 1, ny + 2))
    uaw = np.zeros((nx + 1, ny + 2))
    uan = np.zeros((nx + 1, ny + 2))
    uas = np.zeros((nx + 1, ny + 2))
    uap = np.zeros((nx + 1, ny + 2))
    vae = np.zeros((nx + 2, ny + 1))
    vaw = np.zeros((nx + 2, ny + 1))
    van = np.zeros((nx + 2, ny + 1))
    vas = np.zeros((nx + 2, ny + 1))
    vap = np.zeros((nx + 2, ny + 1))
    pae = np.zeros((nx + 2, ny + 2))
    paw = np.zeros((nx + 2, ny + 2))
    pan = np.zeros((nx + 2, ny + 2))
    pas = np.zeros((nx + 2, ny + 2))
    pap = np.zeros((nx + 2, ny + 2))
    # source term
    pb = np.zeros((nx + 2, ny + 2))

    # for iterations
    usOld = np.zeros((nx + 1, ny + 2))
    vsOld = np.zeros((nx + 2, ny + 1))
    ppOld = np.zeros((nx + 2, ny + 2))
    niter = 2000  # main iteration limit-1
    miter = 100  # v,c,p iterations limit-1
    merr = 0.001  # acceptable error for v,c,p
    nerr = 0.00001  # acceptable error for SIMPLE
    serr = [1]  # error for convergence curve
    uerr = []
    verr = []
    perr = []

    # SIMPLE algorithm
    n = 0
    while (min(serr) > nerr):
        n += 1
        print('Main SIMPLE Iteration#:', n)
        # step1: solving momentum equations
        for i in range(1, nx):
            for j in range(1, ny + 1):
                # convective coefficients CD method
                uace = dy * (uOld[i, j] + uOld[i + 1, j]) / 2
                uacw = dy * (uOld[i, j] + uOld[i - 1, j]) / 2
                uacn = dx * (vOld[i, j] + vOld[i + 1, j]) / 2
                uacs = dx * (vOld[i, j - 1] + vOld[i + 1, j - 1]) / 2
                # coefficients using hybrid
                udx = dy / dx / Re
                udy = dx / dy / Re
                uae[i, j] = max(-uace, (udx - uace / 2), 0)
                uaw[i, j] = max(uacw, (udx + uacw / 2), 0)
                uan[i, j] = max(-uacn, (udy - uacn / 2), 0)
                uas[i, j] = max(uacs, (udy + uacs / 2), 0)
                uap[i, j] = (uae[i, j] + uaw[i, j] + uan[i, j] + uas[i, j]) / rel
                du[i, j] = dy / uap[i, j]
        # putting previous u as new guess
        for i in range(1, nx + 2):
            for j in range(1, ny + 3):
                us[i - 1, j - 1] = uOld[i - 1, j - 1]
        m = 0
        err = 1
        while (err > merr):
            m += 1
            print(n, '-U Iteration#:', m)
            for i in range(1, nx):
                for j in range(1, ny + 1):
                    us[i, j] = (uae[i, j] * us[i + 1, j] + uaw[i, j] * us[i - 1, j] + uan[i, j] * us[i, j + 1] + uas[
                        i, j] * us[i, j - 1] - dy * (ps[i + 1, j] - ps[i, j] + relp * uap[i, j] * uOld[i, j])) / uap[
                                   i, j]
            if m > 1:
                err = error(us, usOld, us, usOld, nx, ny)
                print(n, '-u Error:', err)
                if n > 10:
                    uerr = np.append(uerr, err)
            for i in range(1, nx):
                for j in range(1, ny + 1):
                    usOld[i, j] = us[i, j]
            if m > miter:
                break
                # same thing for v-velocity
        for i in range(1, nx + 1):
            for j in range(1, ny):
                # convective coefficients CD method
                vace = dy * (uOld[i, j] + uOld[i, j + 1]) / 2
                vacw = dy * (uOld[i - 1, j] + uOld[i - 1, j + 1]) / 2
                vacn = dx * (vOld[i, j] + vOld[i, j + 1]) / 2
                vacs = dx * (vOld[i, j] + vOld[i, j - 1]) / 2
                # coefficients using hybrid
                vdx = dy / dx / Re
                vdy = dx / dy / Re
                vae[i, j] = max(-vace, (vdx - vace / 2), 0)
                vaw[i, j] = max(vacw, (vdx + vacw / 2), 0)
                van[i, j] = max(-vacn, (vdy - vacn / 2), 0)
                vas[i, j] = max(vacs, (vdy + vacs / 2), 0)
                vap[i, j] = (vae[i, j] + vaw[i, j] + van[i, j] + vas[i, j]) / rel
                dv[i, j] = dx / vap[i, j]
        # putting previous v as new guess
        for i in range(1, nx + 3):
            for j in range(1, ny + 2):
                vs[i - 1, j - 1] = vOld[i - 1, j - 1]
        m = 0
        err = 1
        while (err > merr):
            m += 1
            print(n, '-V Iteration#:', m)
            for i in range(1, nx + 1):
                for j in range(1, ny):
                    vs[i, j] = (vae[i, j] * vs[i + 1, j] + vaw[i, j] * vs[i - 1, j] + van[i, j] * vs[i, j + 1] + vas[
                        i, j] * vs[i, j - 1] - dy * (ps[i, j + 1] - ps[i, j] + relp * vap[i, j] * vOld[i, j])) / vap[
                                   i, j]
            if m > 1:
                err = error(vs, vsOld, vs, vsOld, nx, ny)
                print(n, '-v Error:', err)
                if n > 10:
                    verr = np.append(verr, err)
            for i in range(1, nx + 1):
                for j in range(1, ny):
                    vsOld[i, j] = vs[i, j]
            if m > miter:
                break

        # step 2:solving pressure correction equation
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                pae[i, j] = du[i, j] * dy
                paw[i, j] = du[i - 1, j] * dy
                pan[i, j] = dv[i, j] * dx
                pas[i, j] = dv[i, j - 1] * dx
                pap[i, j] = pae[i, j] + paw[i, j] + pan[i, j] + pas[i, j]
                pb[i, j] = dy * (us[i - 1, j] - us[i, j]) + dx * (vs[i, j - 1] - vs[i, j])
        if n == 1:
            pap[1, 1] = 1
            pae[1, 1] = 0
            paw[1, 1] = 0
            pan[1, 1] = 0
            pas[1, 1] = 0
            pb[1, 1] = 0
        m = 0
        err = 1
        while (err > merr):
            m += 1
            print(n, '-P Iteration #', m)
            for i in range(1, nx + 1):
                for j in range(1, ny + 1):
                    pp[i, j] = (pae[i, j] * pp[i + 1, j] + paw[i, j] * pp[i - 1, j] + pan[i, j] * pp[i, j + 1] + pas[
                        i, j] * pp[i, j - 1] + pb[i, j]) / pap[i, j]
            if m > 1:
                err = error(pp, ppOld, pp, ppOld, nx, ny)
                print(n, '-p Error:', err)
                if n > 10:
                    perr = np.append(perr, err)
            for i in range(1, nx + 1):
                for j in range(1, ny + 1):
                    ppOld[i, j] = pp[i, j]
            if m > miter:
                break

        # step 3:correcting pressure and velocities
        for i in range(1, nx):
            for j in range(1, ny + 1):
                up[i, j] = du[i, j] * (pp[i, j] - pp[i + 1, j])
        for i in range(1, nx + 1):
            for j in range(1, ny):
                vp[i, j] = dv[i, j] * (pp[i, j] - pp[i, j + 1])
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                p[i, j] = ps[i, j] + pp[i, j] * relp
        for i in range(1, nx):
            for j in range(1, ny + 1):
                u[i, j] = us[i, j] + up[i, j]
        for i in range(1, nx + 1):
            for j in range(1, ny):
                v[i, j] = vs[i, j] + vp[i, j]

        # step 4: check for convergence
        if n > 10:
            err = error(u, uOld, v, vOld, nx, ny)
            print(n, '-SIMPLE Error is:', err)
            serr = np.append(serr, err)

        # step 5: loop
        for i in range(1, nx + 2):
            for j in range(1, ny + 3):
                uOld[i - 1, j - 1] = u[i - 1, j - 1]
        for i in range(1, nx + 3):
            for j in range(1, ny + 2):
                vOld[i - 1, j - 1] = v[i - 1, j - 1]
        for i in range(1, nx + 3):
            for j in range(1, ny + 3):
                ps[i - 1, j - 1] = p[i - 1, j - 1]

        # Boundary Conditions
        uOld[0, :] = 0  # -uOld[1,:] #west
        vOld[0, :] = 0
        p[0, :] = p[1, :]
        uOld[-1, :] = 0  # -uOld[-2,:] #east
        vOld[-1, :] = 0
        p[-1, :] = p[-2, :]
        uOld[:, 0] = 0  # south
        vOld[:, 0] = 0  # -vOld[:,1]
        p[:, 0] = p[:, 1]
        uOld[:, -1] = ulid  # north
        vOld[:, -1] = 0  # -vOld[:,-2];
        p[:, -1] = p[:, -2]
        """
        #for a step in cavity
        for i in range (1,nx//3+1):
            for j in range (1,ny//3+1):
                uOld[nx-i-1,j-1]=0
                vOld[nx-i-2,j-1]=0
                ps[nx-i,j-1]=0

        #for backwardstep/fully developed flow/couette flow
        for j in range (1,ny+1):
            if j>(ny+1)/2:
                uOld[0,j-1]=0
                vOld[0,j-1]=0
            uOld[0,j-1]=ulid
            vOld[0,j-1]=0
            """
        if n > niter:
            break

    return uOld, vOld, p, serr, uerr, verr, perr


def error(u, uOld, v, vOld, nx, ny):
    udiff = np.zeros((nx, ny))
    vdiff = np.zeros((nx, ny))
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            udiff[i - 1, j - 1] = abs(u[i - 1, j - 1] - uOld[i - 1, j - 1])
            vdiff[i - 1, j - 1] = abs(v[i - 1, j - 1] - vOld[i - 1, j - 1])
    umax = max(map(max, udiff))
    vmax = max(map(max, vdiff))
    maxerror = max(umax, vmax)
    return maxerror


def main():
    nx = 65
    ny = nx
    Re = 100
    u, v, p, serr, uerr, verr, perr = (LidCavity(nx, ny, Re))

    import matplotlib.pyplot as plt

    # Create an x, y grid and graph
    cp = plt.pcolor(u.transpose())
    plt.colorbar(cp)
    plt.title('%d by %d u-vel Contour Plot Re= %d' % (nx, ny, Re))
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

    cp = plt.pcolor(v.transpose())
    plt.colorbar(cp)
    plt.title('%d by %d v-vel Contour Plot Re= %d' % (nx, ny, Re))
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

    cp = plt.pcolor(p.transpose())
    plt.colorbar(cp)
    plt.title('%d by %d p Contour Plot Re= %d' % (nx, ny, Re))
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    uu = np.zeros((nx, ny))
    vv = np.zeros((nx, ny))
    for i in range(1, nx):
        for j in range(1, ny):
            uu[i - 1, j - 1] = u[i, j + 1]
            vv[i - 1, j - 1] = v[i + 1, j]
    plt.axis([-0.2, 1.2, 0, 1.2])
    plt.streamplot(x, y, uu.transpose(), vv.transpose())
    plt.title('%d by %d velocity streamline Re= %d' % (nx, ny, Re))
    plt.show()

    y = np.linspace(0, 1, nx + 2)
    plt.plot(y, v[:, 33])
    """plt.plot([0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1],[0,0.09233,0.10091,0.10890,0.12317,0.16077,0.17507,0.17527,0.05454,-0.24533,-0.22445,-0.16914,-0.10313,-0.08864,-0.07391,-0.05906,0],'ro')"""
    plt.axis([0, 1, min(v[:, 33]), max(v[:, 33])])
    plt.title('%d by %d Graph v-vel through horizontal center line Re= %d' % (nx, ny, Re))
    plt.xlabel('x (m)')
    plt.ylabel('velocity (m/s)')
    plt.show()

    serr = np.delete(serr, 0)
    sit = np.linspace(10, 10 + len(serr), len(serr))
    uit = np.linspace(10, 10 + len(serr), len(uerr))
    vit = np.linspace(10, 10 + len(serr), len(verr))
    pit = np.linspace(10, 10 + len(serr), len(perr))
    plt.semilogy(sit, serr)
    plt.semilogy(uit, uerr)
    plt.semilogy(vit, verr)
    plt.semilogy(pit, perr)
    plt.axis([10, 10 + len(serr), 0, max(max(serr), max(perr))])
    plt.title('%d by %d convergence graph Re= %d' % (nx, ny, Re))
    plt.legend(['SIMPLE', 'u-vel', 'v-vel', 'pressure'], loc='upper right')
    plt.xlabel('Iteration #')
    plt.ylabel('Error (log scale)')
    plt.show()


if __name__ == '__main__':
    import time

    start = time.time()
    main()
    end = time.time()
    print('Execution time =', end - start)