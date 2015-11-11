#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Este script genera plots a partir de los datos generados por fluido.py
'''

import bisect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from scipy import interpolate

# Recupera y ordena datos
Data = np.load("Datos.npy")

x = []
t = []
d = []
v = []

x00 = np.linspace(0, 1, 10001)
d00 = np.ones(10001)
d00[0:1001] += 0.0350*(1+np.cos(10*np.pi*(x00[0:1001].copy())))
v00 = np.zeros(10001)

x02, x04, x06, x08, x10 = [], [], [], [], []
d02, d04, d06, d08, d10 = [], [], [], [], []
v02, v04, v06, v08, v10 = [], [], [], [], []
D02, D04, D06, D08, D10 = {}, {}, {}, {}, {}

for i in range(147):
    for j in range(0, 10001, 100):
        x.append(Data[i][0][j])
        t.append(Data[i][1][j])
        v.append(Data[i][2][j])
        d.append(Data[i][3][j])

for i in range(len(Data)):
    for j in range(10001):
        if 0.199 <= Data[i][1][j] <= 0.201:
            D02[Data[i][0][j]] = (Data[i][0][j], Data[i][2][j], Data[i][3][j])
            bisect.insort(x02, Data[i][0][j])
    for j in range(10001):
        if 0.399 <= Data[i][1][j] <= 0.401:
            D04[Data[i][0][j]] = (Data[i][0][j], Data[i][2][j], Data[i][3][j])
            bisect.insort(x04, Data[i][0][j])
    for j in range(10001):
        if 0.599 <= Data[i][1][j] <= 0.601:
            D06[Data[i][0][j]] = (Data[i][0][j], Data[i][2][j], Data[i][3][j])
            bisect.insort(x06, Data[i][0][j])
    for j in range(10001):
        if 0.799 <= Data[i][1][j] <= 0.801:
            D08[Data[i][0][j]] = (Data[i][0][j], Data[i][2][j], Data[i][3][j])
            bisect.insort(x08, Data[i][0][j])
    for j in range(10001):
        if 0.999 <= Data[i][1][j] <= 1.001:
            D10[Data[i][0][j]] = (Data[i][0][j], Data[i][2][j], Data[i][3][j])
            bisect.insort(x10, Data[i][0][j])

for i in range(len(x02)):
    d02.append(D02[x02[i]][2])
    v02.append(D02[x02[i]][1])
for i in range(len(x04)):
    d04.append(D04[x04[i]][2])
    v04.append(D04[x04[i]][1])
for i in range(len(x06)):
    d06.append(D06[x06[i]][2])
    v06.append(D06[x06[i]][1])
for i in range(len(x08)):
    d08.append(D08[x08[i]][2])
    v08.append(D08[x08[i]][1])
for i in range(len(x10)):
    d10.append(D10[x10[i]][2])
    v10.append(D10[x10[i]][1])

x02 = np.array(x02)
x04 = np.array(x04)
x06 = np.array(x06)
x08 = np.array(x08)
x10 = np.array(x10)
d02 = np.array(d02)
d04 = np.array(d04)
d06 = np.array(d06)
d08 = np.array(d08)
d10 = np.array(d10)
v02 = np.array(v02)
v04 = np.array(v04)
v06 = np.array(v06)
v08 = np.array(v08)
v10 = np.array(v10)

# Plots
plt.clf()
fig1 = plt.figure(1)
plt.tripcolor(x, t, d, shading='gouraud')
plt.colorbar()
plt.title('Densidad del fluido sobre el espacio y el tiempo')
plt.xlabel('x')
plt.ylabel('t')
# plt.savefig('Dxy.jpg')
plt.show()

fig2 = plt.figure(2)
plt.tripcolor(x, t, v, shading='gouraud')
plt.colorbar()
plt.title('Velocidad del fluido sobre el espacio y el tiempo')
plt.xlabel('x')
plt.ylabel('t')
# plt.savefig('vxy.jpg')
plt.show()

plt.clf()
plt.figure(3)
plt.plot(x00, d00, color='r')
plt.plot(x00, v00+np.ones(len(v00)), color='b')
plt.title('Velocidad y densidad del fluido, t = 0.0')
plt.xlabel('x')
plt.ylabel('densidad (rojo) y velocidad + 1 (azul)')
plt.savefig('t00.eps')
# plt.show()

plt.clf()
plt.figure(4)
plt.plot(x02, d02, color='r')
plt.plot(x02, v02+np.ones(len(v02)), color='b')
plt.title('Velocidad y densidad del fluido, t = 0.2')
plt.xlabel('x')
plt.ylabel('densidad (rojo) y velocidad + 1 (azul)')
plt.savefig('t02.eps')
# plt.show()

plt.clf()
plt.figure(5)
plt.plot(x04, d04, color='r')
plt.plot(x04, v04+np.ones(len(v04)), color='b')
plt.title('Velocidad y densidad del fluido, t = 0.4')
plt.xlabel('x')
plt.ylabel('densidad (rojo) y velocidad + 1 (azul)')
plt.savefig('t04.eps')
# plt.show()

plt.clf()
plt.figure(6)
plt.plot(x06, d06, color='r')
plt.plot(x06, v06+np.ones(len(v06)), color='b')
plt.title('Velocidad y densidad del fluido, t = 0.6')
plt.xlabel('x')
plt.ylabel('densidad (rojo) y velocidad + 1 (azul)')
plt.savefig('t06.eps')
# plt.show()

plt.clf()
plt.figure(7)
plt.plot(x08, d08, color='r')
plt.plot(x08, v08+np.ones(len(v08)), color='b')
plt.title('Velocidad y densidad del fluido, t = 0.8')
plt.xlabel('x')
plt.ylabel('densidad (rojo) y velocidad + 1 (azul)')
plt.savefig('t08.eps')
# plt.show()

plt.clf()
plt.figure(8)
plt.plot(x10, d10, color='r')
plt.plot(x10, v10+np.ones(len(v10)), color='b')
plt.title('Velocidad y densidad del fluido, t = 1.0')
plt.xlabel('x')
plt.ylabel('densidad (rojo) y velocidad + 1 (azul)')
plt.savefig('t10.eps')
# plt.show()
