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
        # Iterarcion hasta converger, con un maximo de 1000
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
        # Iterarcion hasta converger, con un maximo de 1000
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
        Da un paso en el interior del intervalo
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
        # Iterarcion hasta converger, con un maximo de 1000
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
        
    