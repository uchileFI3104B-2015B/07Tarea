
# coding: utf-8

# In[3]:

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


# P1 Tarea 7

def crea_v_y_rho(xr, cr, tr, vr, rhor, A, numptos, dx, gamma):
    '''Inicializa los vectores v y rho
    con las condiciones iniciales y de borde
    ademas de inicializar cr que depende de rhor'''
    for i in range(numptos):
        xr[i] = i*dx
        if xr[i] > 0.1:
            rhor[i] = 1
            if xr[i] > 1.0/3:
                cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
            elif xr[i] <= 1.0/3:
                cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
        elif xr[i] <= 0.1:
            rhor[i] = 1 + 0.0808*(1 + np.cos(10 * np.pi * xr[i]))
            if xr[i] > 1.0/3:
                cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
            elif xr[i] <= 1.0/3:
                cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
    vr[0] = 0  # primero
    vr[-1] = 0  # ultimo
    return xr, cr, vr, rhor
    '''no devuelve tr porque no hay
    condiciones ni de borde ni iniciales para el tiempo'''


def formulas_x_t_v_rho(xr, cr, tr, vr, rhor, xr_n, tr_n, vr_n,
                       rhor_n, j, par=True, prim_iter=True):
    '''Asigna valores a cada componente de
    cada vector segun las formulas para cada uno'''
    if par is True:
        if prim_iter is True:
            xr_n[j] = ((xr[j-1] * (cr[j] - vr[j]) + (tr[j] -
                        tr[j-1]) * (vr[j-1] * cr[j] + cr[j] ** 2 -
                        vr[j] * cr[j] - vr[j-1] * vr[j]) + (cr[j] +
                        vr[j-1]) * xr[j]) / (vr[j-1] + 2 * cr[j] - vr[j]))
            tr_n[j] = ((-xr[j-1] + vr[j-1] * (tr[j-1] - tr[j]) +
                        cr[j] * (tr[j-1] + tr[j]) + xr[j]) / (vr[j-1] +
                        2 * cr[j] - vr[j]))
            rhor_n[j] = ((rhor[j-1] * rhor[j] * (cr[j] + vr[j-1] +
                          cr[j] - vr[j])) / (rhor[j-1] * cr[j] + rhor[j] *
                          cr[j-1]))
            vr_n[j] = ((rhor[j-1] * cr[j] * (cr[j-1] + vr[j-1]) +
                        cr[j-1] * rhor[j] * (vr[j] - cr[j])) / (rhor[j-1] *
                        cr[j] + rhor[j] * cr[j-1]))
        else:
            xr_n[j] = ((0.5 * (1/(vr[j-1] + cr[j-1] + 2 * cr[j] -
                        vr[i] + cr[i]))) * ((tr[j-1] + tr[i]) *
                        (vr[j-1] * vr[j] - vr[j-1] * cr[j] + vr[j] ** 2 +
                        vr[j] * cr[j-1] - cr[j-1] * cr[j] - cr[j] ** 2) +
                        2 * xr[j-1] * (- vr[j] + cr[j] - vr[i] + cr[i]) +
                        (tr[j-1] * (vr[i] - cr[i]) + 2 * xr[i] + vr[i] *
                        tr[i] - cr[i] * tr[i]) * (vr[j-1] + vr[j] +
                        cr[j-1] + cr[j])))
            tr_n[j] = ((1 / (vr[j-1] + cr[j-1] + 2 * cr[j] - vr[i] +
                        cr[i])) * (-2 * xr[j-1] + tr[j-1] * (vr[j-1] +
                        vr[j] + cr[j-1] + cr[j]) + tr[i]*(vr[j] -
                        cr[j] + vr[i] - cr[i]) + 2 * xr[i]))
            rhor_n[j] = ((-1.0 * (rhor[j-1] + rhor[j]) * (vr[i] *
                          (rhor[i] + rhor[j]) - rhor[i] * (cr[i] +
                          cr[j])) + (rhor[i] + rhor[j]) * (-rhor[j-1] *
                          (cr[j-1] + cr[j]) - vr[j-1] * (rhor[j-1] +
                          rhor[j]))) / ((cr[i] + cr[j]) * (rhor[j-1] +
                          rhor[j]) + (cr[j-1] + cr[j]) * (rhor[i] +
                          rhor[j])))
            vr_n[j] = (((cr[j-1] + cr[j]) / (rhor[j-1] + rhor[j])) *
                         rhor_n[j] + (rhor[j-1] * (cr[j-1] + cr[j]) +
                         vr[j-1] * (rhor[j-1] + rhor[j])) / (rhor[j-1] +
                         rhor[j]))
    else:
        if prim_iter is True:
            xr_n[j] = ((xr[j] * (cr[j+1] - vr[j+1]) + (tr[j+1] -
                        tr[j]) * (vr[j] * cr[j+1] + cr[j+1] ** 2 -
                        vr[j+1] * cr[j+1] - vr[j] * vr[j+1]) +
                        (cr[j+1] + vr[j]) * xr[j+1]) / (vr[j] +
                        2 * cr[j+1] - vr[j+1]))
            tr_n[j] = ((-xr[j] + vr[j] * (tr[j] - tr[j+1]) +
                        cr[j+1] * (tr[j] + tr[j+1]) + xr[j+1]) /
                        (vr[j] + 2 * cr[j+1] - vr[j+1]))
            rhor_n[j] = ((rhor[j] * rhor[j+1] * (cr[j] + vr[j] + cr[j+1] -
                          vr[j+1])) / (rhor[j] * cr[j+1] +
                          rhor[j+1] * cr[j]))
            vr_n[j] = ((rhor[j] * cr[j+1] * (cr[j] + vr[j]) + cr[j] *
                        rhor[j+1] * (vr[j+1] - cr[j+1])) / (rhor[j] *
                        cr[j+1] + rhor[j+1] * cr[j]))
        else:
            xr_n[j] = ((0.5 * (1 / (vr[j] + cr[j] + 2 * cr[j] - vr[j+1] +
                        cr[j+1]))) * ((tr[j] + tr[j+1]) * (vr[j] * vr[j] -
                        vr[j] * cr[j] + vr[j] ** 2 + vr[j] * cr[j] -
                        cr[j] * cr[j] - cr[j] ** 2) + 2 * xr[j] * (- vr[j] +
                        cr[j] - vr[j+1] + cr[j+1]) + (tr[j] * (vr[j+1] -
                        cr[j+1]) + 2 * xr[j+1] + vr[j+1] * tr[j+1] -
                        cr[j+1] * tr[j+1]) * (vr[j] + vr[j] + cr[j] + cr[j])))
            tr_n[j] = ((1 / (vr[j] + cr[j] + 2 * cr[j] - vr[j+1] +
                        cr[j+1])) * (- 2 * xr[j] + tr[j] * (vr[j] + vr[j] +
                        cr[j] + cr[j]) + tr[j+1] * (vr[j] - cr[j] +
                        vr[j+1] - cr[j+1]) + 2 * xr[j+1]))
            rhor_n[j] = ((-1.0 * (rhor[j] + rhor[j]) * (vr[j+1] * (rhor[j+1] +
                          rhor[j]) - rhor[j+1] * (cr[j+1] + cr[j])) +
                          (rhor[j+1] + rhor[j]) * (-rhor[j] * (cr[j] +
                          cr[j]) - vr[j] * (rhor[j] + rhor[j]))) / ((cr[j+1] +
                          cr[j]) * (rhor[j] + rhor[j]) + (cr[j] + cr[j]) *
                          (rhor[j+1] + rhor[j])))
            vr_n[j] = (((cr[j] + cr[j]) / (rhor[j] + rhor[j])) * rhor_n[j] +
                         (rhor[j] * (cr[j] + cr[j]) + vr[j] * (rhor[j] +
                         rhor[j])) / (rhor[j] + rhor[j]))


def tstep(xr, cr, tr, vr, rhor, A, numptos, dx, i, prim_iter,
          gamma):  # contiene condiciones iniciales y de borde tambien
    '''Avanza un paso temporal las soluciones xr, cr, vr y rhor'''
    xr_n, tr_n, vr_n, rhor_n = (np.zeros(numptos), np.zeros(numptos),
                                np.zeros(numptos), np.zeros(numptos))
    # Por si acaso
    vr[0] = 0  # primero
    vr[-1] = 0  # ultimo
    if i % 2 == 0:  # pares
        if prim_iter is True:
            for j in range(1, numptos):  # no toma el indice 0
                formulas_x_t_v_rho(xr, cr, tr, vr, rhor, xr_n, tr_n,
                                   vr_n, rhor_n, j)
                for i in range(numptos):
                    if xr[i] > 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
                    elif xr[i] <= 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
        else:
            for j in range(1, numptos):
                formulas_x_t_v_rho(xr, cr, tr, vr, rhor, xr_n, tr_n,
                                   vr_n, rhor_n, j, prim_iter=False)
                for i in range(numptos):
                    if xr[i] > 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
                    elif xr[i] <= 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
    elif i % 2 != 0:  # impares
        if prim_iter is True:
            for j in range(0, numptos-1):  # no toma el indice -1
                formulas_x_t_v_rho(xr, cr, tr, vr, rhor, xr_n, tr_n,
                                   vr_n, rhor_n, j, par=False)
                for i in range(numptos):
                    if xr[i] > 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
                    elif xr[i] <= 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
        else:
            for j in range(0, numptos-1):
                formulas_x_t_v_rho(xr, cr, tr, vr, rhor, xr_n, tr_n,
                                   vr_n, rhor_n, j, par=False,
                                   prim_iter=False)
                for i in range(numptos):
                    if xr[i] > 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
                    elif xr[i] <= 0.1:
                        if xr[i] > 1.0/3:
                            cr[i] = np.sqrt(A[1] * gamma * rhor[i] ** (gamma-1))
                        elif xr[i] <= 1.0/3:
                            cr[i] = np.sqrt(A[0] * gamma * rhor[i] ** (gamma-1))
    return xr_n, tr_n, vr_n, rhor_n


def prim_iter(numpasost, numptos, xr, cr, tr, vr, rhor, A, dx,
              prim_iter, gamma):
    '''Itera sobre la grilla en un tiempo fijo
    para la primera iteracion'''
    prim_iter = True
    X = []
    T = []
    RHO = []
    V = []
    C = []
    for i in range(1, numpasost):
        a, b, c, d = tstep(xr, cr, tr, vr, rhor, A, numptos, dx,
                           i, prim_iter, gamma)
        xr, tr, vr, rhor = a, b, c, d
        C.append(cr)
        X.append(xr)
        T.append(tr)
        RHO.append(rhor)
        V.append(vr)
    X, T, RHO, V = (np.array(X), np.array(T), np.array(RHO),
                    np.array(V))
    return X, T, RHO, V, C


def segda_iter(numpasost, numptos, xr, cr, tr, vr, rhor, A, dx,
               prim_iter, gamma):
    '''Itera sobre la grilla en un tiempo fijo
    para las iteraciones despues de la primera'''
    prim_iter = False
    X = []
    T = []
    RHO = []
    V = []
    C = []
    for i in range(1, numpasost):
        a, b, c, d = tstep(xr, cr, tr, vr, rhor, A, numptos, dx,
                           i, prim_iter, gamma)
        xr, tr, vr, rhor = a, b, c, d
        C.append(cr)
        X.append(xr)
        T.append(tr)
        RHO.append(rhor)
        V.append(vr)
    X, T, RHO, V = (np.array(X), np.array(T), np.array(RHO),
                    np.array(V))
    return X, T, RHO, V, C


# Main setup

A = [4, 1]
gamma = 5.0/3.0
numsol = 5
numptos = 10001
numpasost = 5
dx = 1.0/10000
xr = np.zeros(numptos)
cr = np.zeros(numptos)
tr = np.zeros(numptos)
vr = np.zeros(numptos)
rhor = np.zeros(numptos)
Xsol, Tsol, RHOsol, Vsol, Csol = [], [], [], [], []
a, b, c, d = crea_v_y_rho(xr, cr, tr, vr, rhor, A, numptos,
                          dx, gamma)
xr, cr, vr, rhor = a, b, c, d
X, T, RHO, V, C = prim_iter(numpasost, numptos, xr, cr, tr, vr, rhor,
                            A, dx, prim_iter, gamma)
Xsol.append(X)
Tsol.append(T)
RHOsol.append(RHO)
Vsol.append(V)
Csol.append(C)
xr2, tr2, vr2, rhor2, cr2 = (Xsol[0][-1], Tsol[0][-1], Vsol[0][-1],
                             RHOsol[0][-1], Csol[0][-1])
for i in range(numsol):
    X, T, RHO, V, C = segda_iter(numpasost, numptos, xr2, cr2, tr2, vr2,
                                 rhor2, A, dx, prim_iter, gamma)
    xr2, cr2, tr2, vr2, rhor2 = (X[0][-1], C[0][-1], T[0][-1], V[0][-1],
                                 RHO[0][-1])
    Xsol.append(X)
    Tsol.append(T)
    RHOsol.append(RHO)
    Vsol.append(V)
print Xsol

# No pasa el pep8 por la cantidad de parentesis y
# no puedo arreglarlo


# In[ ]:



