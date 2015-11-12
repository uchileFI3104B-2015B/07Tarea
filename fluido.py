# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Este script define la clase Fluido, y un #Main que realiza la iteración para
crear los datos requeridos
"""

import numpy as np
import bisect

class Fluido(object):
    '''
    Define un fluido compresible en 1D, que contiene una curva (x,y),
    la velocidad v en cada paso, y la densidad rho=d en cada paso.
    '''
    def __init__(self):
        '''
        Inicializa el fluido con las condiciones iniciales especificadas.
        '''
        RUT = 0.0370
        self.x = np.linspace(0, 1, 10001)
        self.t = np.zeros(10001)
        self.d = np.ones(10001)
        self.d[0:1001] += RUT*(1+np.cos(10*np.pi*(self.x[0:1001].copy())))
        self.v = np.zeros(10001)
        self.c = np.sqrt((5.0/3.0)*np.power((self.d.copy()), 2.0/3.0))
        self.c[0:3334] = 2*(self.c[0:3334].copy())
        self.paso_par = True
        self.N = 10001
    
    def paso_izquierdo(self, Q, ver_convergencia):
        '''
        Da un paso hacia el borde izquierdo (x=0)
        '''
        xQ, tQ, vQ, dQ, cQ = Q
        xR = 0
        vR = 0
        tR = (xR-xQ)/(vQ-cQ) + tQ
        dR = (vR-vQ)*dQ/cQ + dQ
        cR = np.sqrt(4*(5.0/3.0)*(dR**(2.0/3.0)))
        if not ver_convergencia:
            return (xR, tR, vR, dR, cR)
        for cnt in range(1000):
            tRn = 2.0*(xR-xQ)/(vQ+vR-cQ-cR) + tQ
            dRn = (vR-vQ)*(dQ+dR)/(cQ+cR) + dQ
            cRn = np.sqrt(4*(5.0/3.0)*(dRn**(2.0/3.0)))
            if abs(dRn-dR)/abs(dR) < 10**(-8):
                if abs(tRn-tR)/abs(tR) < 10**(-8):
                    break
            tR, dR, cR = tRn, dRn, cRn
        R = (xR, tRn, vR, dRn, cRn)
        return R
    
    def paso_derecho(self, P, ver_convergencia):
        '''
        Da un paso hacia el borde derecho (x=1)
        '''
        xP, tP, vP, dP, cP = P
        xR = 1
        vR = 0
        tR = (xR-xP)/(vP+cP) + tP
        dR = -(vR-vP)*dP/cP + dP
        cR = np.sqrt((5.0/3.0)*(dR**(2.0/3.0)))
        if not ver_convergencia:
            return (xR, tR, vR, dR, cR)
        for cnt in range(1000):
            tRn = 2.0*(xR-xP)/(vP+vR+cP+cR) + tP
            dRn = -(vR-vP)*(dP+dR)/(cP+cR) + dP
            cRn = np.sqrt((5.0/3.0)*(dRn**(2.0/3.0)))
            if abs(dRn-dR)/abs(dR) < 10**(-8):
                if abs(tRn-tR)/abs(tR) < 10**(-8):
                    break
            tR, dR, cR = tRn, dRn, cRn
        R = (xR, tRn, vR, dRn, cRn)
        return R
        
    def paso_principal(self, P, Q, ver_convergencia):
        '''
        Da un paso en el interior del intervalo.
        '''
        xQ, tQ, vQ, dQ, cQ = Q
        xP, tP, vP, dP, cP = P
        tR = (xQ-(vQ-cQ)*tQ-xP+(vP+cP)*tP)/(vP+cP+cQ-vQ)
        xR = (vP+cP)*(tR-tP)+xP
        dR = (vP-vQ+cQ+cP)/(cQ/dQ+cP/dP)
        vR = (cQ/dQ)*(dR-dQ)+vQ
        if xR <= 1.0/3.0:
            cR = np.sqrt(4*(5.0/3.0)*(dR**(2.0/3.0)))
        else:
            cR = np.sqrt((5.0/3.0)*(dR**(2.0/3.0)))
        if not ver_convergencia:
            R = (xR, tR, vR, dR, cR)
            if xR > 1:
                R = self.paso_derecho(P, ver_convergencia)
            if xR < 0:
                R = self.paso_izquierdo(Q, ver_convergencia)
            return R
        for cnt in range(1000):
            tRn = ((xQ - ((vQ+vR)/2.0 - (cQ+cR)/2.0)*tQ - xP +
                    ((vP+vR)/2.0 + (cP+cR)/2.0)*tP)/(((vP+vR)+(cP+cR)+(cQ+cR) -
                                                      (vQ+vR))/2.0))
            xRn = (vP+vR+cP+cR)*(tRn-tP)/2.0 + xP
            dRn = ((vP-vQ+dP*(cP+cR)/(dP+dR)+dQ*(cQ+cR)/(dQ+dR)) /
                   ((cP+cR)/(dP+dR)+(cQ+cR)/(dQ+dR)))
            vRn = (dRn-dQ)*(cQ+cR)/(dQ+dR) + vQ
            if xRn <= 1.0/3.0:
                cRn = np.sqrt(4*(5.0/3.0)*dRn**(2.0/3.0))
            else:
                cRn = np.sqrt((5.0/3.0)*dRn**(2.0/3.0))
            if abs(xRn-xR)/abs(xR) < 10**(-5):
                if abs(tRn-tR)/abs(tR) < 10**(-5):
                    if abs(dRn-dR)/abs(dR) < 10**(-5):
                        if abs(vRn-vR)/(abs(vR)+0.001) < 10**(-5):
                            break
            xR, tR, vR, dR, cR = xRn, tRn, vRn, dRn, cRn
        R = (xRn, tRn, vRn, dRn, cRn)
        if xRn > 1:
            R = self.paso_derecho(P, ver_convergencia)
        if xRn < 0:
            R = self.paso_izquierdo(Q, ver_convergencia)
        return R
        
    def avance_t(self, ver_convergencia):
        '''
        Avanza la curva con los datos del fluido en un paso temporal,
        siguiendo las diferencias entre pasos pares e impares.
        '''
        if self.paso_par:
            D = {}
            L = []
            for i in range(self.N - 1):
                P = (self.x[i], self.t[i], self.v[i], self.d[i], self.c[i])
                Q = (self.x[i+1], self.t[i+1], self.v[i+1], self.d[i+1],
                     self.c[i+1])
                R = self.paso_principal(P, Q, ver_convergencia)
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
            for i in range(self.N - 1): # Ordenar los puntos
                self.x[i] = L[i]
                self.t[i] = D[L[i]][1]
                self.v[i] = D[L[i]][2]
                self.d[i] = D[L[i]][3]
                self.c[i] = D[L[i]][4]
            self.paso_par = False
        else:
            # Borde derecho
            P = (self.x[self.N - 2], self.t[self.N - 2], self.v[self.N - 2],
                 self.d[self.N - 2], self.c[self.N - 2])
            R = self.paso_derecho(P, ver_convergencia)
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
                R = self.paso_principalp(P, Q, ver_convergencia)
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
            for i in range(self.N - 2, 0, -1): #Ordenar los puntos
                self.x[i] = L[i-1]
                self.t[i] = D[L[i-1]][1]
                self.v[i] = D[L[i-1]][2]
                self.d[i] = D[L[i-1]][3]
                self.c[i] = D[L[i-1]][4]
            # Borde izquierdo
            Q = (self.x[0], self.t[0], self.v[0], self.d[0], self.c[0])
            R = self.paso_izquierdo(Q, ver_convergencia)
            self.x[0] = R[0]
            self.t[0] = R[1]
            self.v[0] = R[2]
            self.d[0] = R[3]
            self.c[0] = R[4]
            self.paso_par = True
            
    def estado_actual(self):
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

