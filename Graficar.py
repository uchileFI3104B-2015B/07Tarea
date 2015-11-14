#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Este script grafica cargando datos entregados por Tarea07.py
'''
from __future__ import division
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import interpolate
Datos= np.load('Datos.npy')

X=[]
T=[]
V=[]
Rho=[]
for i in range(len(Datos)):
    for j in range(0,len(Datos[i][0]),100):
        X.append((Datos[i][0])[j])
        T.append((Datos[i][1])[j])
        V.append((Datos[i][2])[j])
        Rho.append((Datos[i][3])[j])


fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
plt.tripcolor(X, T, Rho, shading='gouraud')
plt.colorbar()
plt.title('Densidad del gas en el plano posicion-tiempo')
plt.xlabel('Posicion')
plt.ylabel('Tiempo')
plt.savefig('Densidad_plano.jpg')
plt.show()

fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
plt.tripcolor(X, T, V, shading='gouraud')
plt.colorbar()
plt.title('Velocidad del gas en el plano posicion-tiempo')
plt.xlabel('Posicion')
plt.ylabel('Tiempo')
plt.savefig('Velocidad_plano.jpg')
plt.show()


#fig3 = plt.figure(3)
#fig3.clf()
#ax3 = fig3.add_subplot(111)
#f = interpolate.interp2d(X, T, Rho, kind='cubic')
#xnew = np.arange(0, 1.001, 0.001)
#ynew = np.arange(0, 1.001, 0.001)
#znew = f(xnew, ynew)
#plt.tripcolor(xnew, ynew,znew, shading='gouraud')
#ax.pcolormesh(xnew, ynew,znew )
#plt.show()
