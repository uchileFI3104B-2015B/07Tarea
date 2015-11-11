#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Este script define la clase GasIdeal e implementa una ejecución del mismo
sobre sus curvas características.
'''

import bisect
import numpy as np
import matplotlib.pyplot as plt


class GasIdeal(object):
    '''
    Esta clase crea el objeto que contiene la información de un gas ideal
    unidimensional dentro de la región (0,1). En cada paso contiene una
    curva (x,t) con la información de velocidad v y densidad d en cada punto.
    '''
    def __init__(self):
        '''
        Inicializa gas con condiciones iniciales especificadas
        '''
        self.x = np.linspace(0, 1, 10001)
        self.t = np.zeros(10001)
        self.d = np.ones(10001)
        self.d[0:1001] += 0.0350*(1+np.cos(10*np.pi*(self.x[0:1001].copy())))
        self.v = np.zeros(10001)
        self.c = (5.0/3.0)*np.power((self.d.copy()), 2.0/3.0)
        self.c[0:3334] = 4*(self.c[0:3334].copy())
        self.paso_par = True
        self.N = 10001

    def cargar_condicion_inicial(self, datos):
        '''
        Carga una condición inicial diferente. Permite testear el algoritmo.
        '''
        self.x = datos[0]
        self.t = datos[1]
        self.v = datos[2]
        self.d = datos[3]
        self.c = (5.0/3.0)*np.power((self.d.copy()), 2.0/3.0)
        self.c[0:3334] = 4*(self.c[0:3334].copy())
        self.paso_par = True
        self.N = 10001

    def borde_izq_step(self, Q, ver_convergencia):
        '''
        Da un paso hacia el borde izquierdo (x=0)
        '''
        xQ, tQ, vQ, dQ, cQ = Q
        xR = 0
        vR = 0
        tR = (xR-xQ)/(vQ-cQ) + tQ
        dR = (vR-vQ)*dQ/cQ + dQ
        if dR <= 0:
            print "Alerta: densidad negativa o nula: x =", xR, "t =", tR,
            print "v = ", vR, "densidad =", dR
        cR = 4*(5.0/3.0)*(dR**(2.0/3.0))
        if abs(vR) > cR:
            print "Alerta: v mayor que c"
        if not ver_convergencia:
            return (xR, tR, vR, dR, cR)
        # Iterar hasta converger, con un maximo de 1000
        for cnt in range(1000):
            tRn = 2.0*(xR-xQ)/(vQ+vR-cQ-cR) + tQ
            dRn = (vR-vQ)*(dQ+dR)/(cQ+cR) + dQ
            cRn = 4*(5.0/3.0)*(dRn**(2.0/3.0))
            if abs(dRn-dR)/abs(dR) < 10**(-8):
                if abs(tRn-tR)/abs(tR) < 10**(-8):
                    break
            tR, dR, cR = tRn, dRn, cRn
        R = (xR, tRn, vR, dRn, cRn)
        return R

    def borde_der_step(self, P, ver_convergencia):
        '''
        Da un paso hacia el borde derecho (x=1)
        '''
        xP, tP, vP, dP, cP = P
        xR = 1
        vR = 0
        tR = (xR-xP)/(vP+cP) + tP
        dR = -(vR-vP)*dP/cP + dP
        if dR <= 0:
            print "Alerta: densidad negativa o nula: x =", xR, "t =", tR,
            print "v = ", vR, "densidad =", dR
        cR = (5.0/3.0)*(dR**(2.0/3.0))
        if abs(vR) > cR:
            print "Alerta: v mayor que c"
        if not ver_convergencia:
            return (xR, tR, vR, dR, cR)
        # Iterar hasta converger, con un maximo de 1000
        for cnt in range(1000):
            tRn = 2.0*(xR-xP)/(vP+vR+cP+cR) + tP
            dRn = -(vR-vP)*(dP+dR)/(cP+cR) + dP
            cRn = (5.0/3.0)*(dRn**(2.0/3.0))
            if abs(dRn-dR)/abs(dR) < 10**(-8):
                if abs(tRn-tR)/abs(tR) < 10**(-8):
                    break
            tR, dR, cR = tRn, dRn, cRn
        R = (xR, tRn, vR, dRn, cRn)
        return R

    def bulk_step(self, P, Q, ver_convergencia):
        '''
        Da un paso en el interior del intervalo
        '''
        xQ, tQ, vQ, dQ, cQ = Q
        xP, tP, vP, dP, cP = P
        tR = (xQ-(vQ-cQ)*tQ-xP+(vP+cP)*tP)/(vP+cP+cQ-vQ)
        xR = (vP+cP)*(tR-tP)+xP
        dR = (vP-vQ+cQ+cP)/(cQ/dQ+cP/dP)
        if dR <= 0:
            print "Alerta: densidad negativa o nula: x =", xR, "t =", tR,
            print "v = ", vR, "densidad =", dR
        vR = (cQ/dQ)*(dR-dQ)+vQ
        if xR <= 1.0/3.0:
            cR = 4*(5.0/3.0)*(dR**(2.0/3.0))
        else:
            cR = (5.0/3.0)*(dR**(2.0/3.0))
        if not ver_convergencia:
            R = (xR, tR, vR, dR, cR)
            if xR > 1:
                R = self.borde_der_step(P, ver_convergencia)
            if xR < 0:
                R = self.borde_izq_step(Q, ver_convergencia)
            return R
        # Iterar hasta converger, con un maximo de 1000
        for cnt in range(1000):
            tRn = ((xQ - ((vQ+vR)/2.0 - (cQ+cR)/2.0)*tQ - xP +
                    ((vP+vR)/2.0 + (cP+cR)/2.0)*tP)/(((vP+vR)+(cP+cR)+(cQ+cR) -
                                                      (vQ+vR))/2.0))
            xRn = (vP+vR+cP+cR)*(tRn-tP)/2.0 + xP
            dRn = ((vP-vQ+dP*(cP+cR)/(dP+dR)+dQ*(cQ+cR)/(dQ+dR)) /
                   ((cP+cR)/(dP+dR)+(cQ+cR)/(dQ+dR)))
            vRn = (dRn-dQ)*(cQ+cR)/(dQ+dR) + vQ
            if xRn <= 1.0/3.0:
                cRn = 4*(5.0/3.0)*dRn**(2.0/3.0)
            else:
                cRn = (5.0/3.0)*dRn**(2.0/3.0)
            if abs(xRn-xR)/abs(xR) < 10**(-5):
                if abs(tRn-tR)/abs(tR) < 10**(-5):
                    if abs(dRn-dR)/abs(dR) < 10**(-5):
                        if abs(vRn-vR)/(abs(vR)+0.001) < 10**(-5):
                            break
            xR, tR, vR, dR, cR = xRn, tRn, vRn, dRn, cRn
        R = (xRn, tRn, vRn, dRn, cRn)
        if xRn > 1:
            R = self.borde_der_step(P, ver_convergencia)
        if xRn < 0:
            R = self.borde_izq_step(Q, ver_convergencia)
        if abs(vRn) > cRn:
            print "Alerta: v mayor que c"
        return R

    def avanza_t(self, ver_convergencia):
        '''
        Avanza la curva con los datos del gas una iteración.
        '''
        # Si el paso que toca es par
        if self.paso_par:
            D = {}
            L = []
            # Calcular los puntos de la siguiente curva
            for i in range(self.N - 1):
                P = (self.x[i], self.t[i], self.v[i], self.d[i], self.c[i])
                Q = (self.x[i+1], self.t[i+1], self.v[i+1], self.d[i+1],
                     self.c[i+1])
                R = self.bulk_step(P, Q, ver_convergencia)
                if R[0] > 1 or R[0] < 0:
                    print "x fuera de rango (paso par)", i
                    print "P =", P
                    print "Q =", Q
                    print "R =", R
                if R[0] in D:
                    continue
                D[R[0]] = R
                bisect.insort(L, R[0])
            self.N = len(L) + 1
            # Ordenar los puntos (si el método resulta estable, esto es
            # redundante)
            for i in range(self.N - 1):
                self.x[i] = L[i]
                self.t[i] = D[L[i]][1]
                self.v[i] = D[L[i]][2]
                self.d[i] = D[L[i]][3]
                self.c[i] = D[L[i]][4]
            self.paso_par = False
        # Si el paso que toca es impar
        else:
            # Calcular los puntos de la siguiente curva:
            # Borde del lado derecho
            P = (self.x[self.N - 2], self.t[self.N - 2], self.v[self.N - 2],
                 self.d[self.N - 2], self.c[self.N - 2])
            R = self.borde_der_step(P, ver_convergencia)
            self.x[self.N - 1] = R[0]
            self.t[self.N - 1] = R[1]
            self.v[self.N - 1] = R[2]
            self.d[self.N - 1] = R[3]
            self.c[self.N - 1] = R[4]
            # Interior
            D = {}
            L = []
            for i in range(self.N - 2, 0, -1):
                P = (self.x[i-1], self.t[i-1], self.v[i-1], self.d[i-1],
                     self.c[i-1])
                Q = (self.x[i], self.t[i], self.v[i], self.d[i], self.c[i])
                R = self.bulk_step(P, Q, ver_convergencia)
                if R[0] > 1 or R[0] < 0:
                    print "x fuera de rango (paso impar)", i
                    print "P =", P
                    print "Q =", Q
                    print "R =", R
                if R[0] in D or R[0] == 1.0 or R[0] == 0.0:
                    continue
                D[R[0]] = R
                bisect.insort(L, R[0])
            self.N = len(L) + 2
            # Ordenar los puntos (si el método resulta estable, esto es
            # redundante)
            for i in range(self.N - 2, 0, -1):
                self.x[i] = L[i-1]
                self.t[i] = D[L[i-1]][1]
                self.v[i] = D[L[i-1]][2]
                self.d[i] = D[L[i-1]][3]
                self.c[i] = D[L[i-1]][4]
            # Borde del lado izquierdo
            Q = (self.x[0], self.t[0], self.v[0], self.d[0], self.c[0])
            R = self.borde_izq_step(Q, ver_convergencia)
            self.x[0] = R[0]
            self.t[0] = R[1]
            self.v[0] = R[2]
            self.d[0] = R[3]
            self.c[0] = R[4]
            self.paso_par = True

    def get_estado(self):
        '''
        Recupera estado actual del gas
        '''
        estado = [0, 0, 0, 0]
        if self.paso_par:
            estado[0] = self.x[0:self.N].copy()
            estado[1] = self.t[0:self.N].copy()
            estado[2] = self.v[0:self.N].copy()
            estado[3] = self.d[0:self.N].copy()
        else:
            estado[0] = self.x[0:(self.N - 1)].copy()
            estado[1] = self.t[0:(self.N - 1)].copy()
            estado[2] = self.v[0:(self.N - 1)].copy()
            estado[3] = self.d[0:(self.N - 1)].copy()
        return estado

if __name__ == '__main__':
    G = GasIdeal()
    E = G.get_estado()
    Data = []
    Data.append(E)
    t_act = 0
    counter = 0
    counter2 = 0
    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(111)
    while t_act <= 1.2:
        G.avanza_t(True)
        E = G.get_estado()
        t_act = E[1][5000]
        counter += 1
        counter2 += 1
        if counter == 100:
            print "Tiempo al centro del intervalo (0,1):", t_act
            print "Numero de posiciones: ", G.N
            Data.append(E)
            counter = 0
        if counter2 == 1000:
            ax3.plot(E[0], E[3])
            counter2 = 0
    plt.show()
    np.save('Datos.npy', Data)
