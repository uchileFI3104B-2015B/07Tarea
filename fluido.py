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
    la velocidad v en cada paso, y la densidad rho en cada paso.
    '''
    def __init__(self):
        '''
        Inicializa el fluido con las condiciones iniciales especificadas.
        '''
        RUT = 0.0370
        self.x = np.linspace(0, 1, 10001)
        self.t = np.zeros(10001)
        self.rho = np.ones(10001)
        self.rho[0:1001] += RUT*(1+np.cos(10*np.pi*(self.x[0:1001].copy())))
        self.v = np.zeros(10001)
        self.c = np.sqrt((5.0/3.0)*np.power((self.rho.copy()), 2.0/3.0))
        self.c[0:3334] = 2*(self.c[0:3334].copy())
        self.paso_par = True
        self.N = 10001