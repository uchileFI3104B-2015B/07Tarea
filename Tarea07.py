#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Este script es para resolver las ecuaciones de un fluido compresible en una dimension
sin viscosidad ni gravedad con el metodo de las caracteristicas
'''
from __future__ import division
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#Funcion de iniciacion
def inicializa(gamma):
    for i in range(Pasos_x):
        if X[i]<=0.1:
            Rho[i] = 1+0.0350*(1+np.cos(10*np.pi*X[i]))
        elif X[i]>0.1:
            Rho[i] = 1
    for i in range(Pasos_x):
        if X[i]<=1/3:
            C[i]=np.sqrt(gamma*4*Rho[i]**(gamma-1))
        elif X[i]>1/3:
            C[i]=np.sqrt(gamma*1*Rho[i]**(gamma-1))
    return Rho, C

#Funcion para los puntos intermedios
def intermedios(X,T,V,Rho,C,paso,Pasos_x,gamma=5/3):
    X_plus = np.array(())
    V_plus = np.array(())
    T_plus = np.array(())
    C_plus = np.array(())
    Rho_plus = np.array(())

    ##agregar bordes
    Tr, Xr, Vr, Rhor, Cr = bordes(X,T,V,Rho,C,paso,0,0,gamma=5/3)
    X_plus = np.append(X_plus, Xr)
    V_plus = np.append(V_plus, Vr)
    T_plus = np.append(T_plus, Tr)
    Rho_plus = np.append(Rho_plus, Rhor)
    C_plus = np.append(C_plus, Cr)

    if paso=="par":
        #print "entre par"
        for i in range(Pasos_x-1):
            Tr = (X[i+1]-X[i]-(V[i+1]-C[i+1])*T[i+1] + (V[i]+C[i])*T[i]) / (V[i]+C[i]-V[i+1]+C[i+1])
            Xr = (V[i]+C[i])*(Tr-T[i]) + X[i]
            Vr = (-Rho[i+1]+Rho[i]+ (V[i+1]*Rho[i+1]/C[i+1]) + (V[i]*Rho[i]/C[i]) ) / (Rho[i+1]/C[i+1] + Rho[i]/C[i])
            Rhor = (Vr - V[i+1])*Rho[i+1]/C[i+1] + Rho[i+1]
            ###print "Tr", X[i+1],X[i],V[i+1],C[i+1],T[i+1],V[i],C[i],T[i]#Tr esta dando raro
            if Xr<=1/3:
                Cr=np.sqrt(np.absolute(gamma*4*Rhor**(gamma-1)))
            elif Xr>1/3:
                Cr=np.sqrt(np.absolute(gamma*1*Rhor**(gamma-1)))
            for j in range(5):
                #print "que sueño",(V[i]+Vr+C[i]+Cr-V[i+1]-Vr+C[i+1]+Cr)
                Tr_plus = 2*( X[i+1]-X[i]+(V[i]+Vr+C[i]+Cr)*(T[i] / 2) - (V[i+1]+Vr-C[i+1]-Cr)*(T[i+1] / 2) ) / (V[i]+Vr+C[i]+Cr-V[i+1]-Vr+C[i+1]+Cr)
                Xr_plus = (1/2)*(V[i]+Vr+C[i]+Cr)*(Tr_plus-T[i])+X[i]
                Vr_plus =( Rho[i] - Rho[i+1] + V[i+1]*(Rho[i+1]+Rhor)/(C[i+1]+Cr) + V[i]*(Rho[i]+Rhor)/(C[i]+Cr) )/ ( ((Rho[i+1]+Rhor)/(C[i+1]+Cr)) + ((Rho[i]+Rhor)/(C[i]+Cr)) )
                Rhor_plus = (Vr_plus-V[i+1])*(Rho[i+1]+Rhor)/(C[i+1]+Cr) +Rho[i+1]
                if Rhor_plus <0:
                    continue
                #print "ror",Rhor_plus
                if Xr_plus<=1/3:
                    Cr_plus=np.sqrt(np.absolute(gamma*4*Rhor_plus**(gamma-1)))
                elif Xr_plus>1/3:
                    Cr_plus=np.sqrt(np.absolute(gamma*1*Rhor_plus**(gamma-1)))
                #if np.isnan(Tr_plus)==True or np.isnan(Xr_plus)==True or np.isnan(Vr_plus)==True or np.isnan(Rhor_plus)==True or np.isnan(Cr)==True :
                    #print "se hizo cero",Rhor,Rho[i+1],(Vr - V[i+1])
                if np.absolute(Xr_plus-Xr)/Xr <10**(-10):
                    Tr, Xr, Vr, Rhor, Cr = Tr_plus, Xr_plus, Vr_plus, Rhor_plus, Cr_plus
                    #print "quebraré", j
                    break

                Tr, Xr, Vr, Rhor, Cr = Tr_plus, Xr_plus, Vr_plus, Rhor_plus, Cr_plus

            #bordes
            if Xr<=0 or Xr>=1:
                continue
#                Tr, Xr, Vr, Rhor, Cr = bordes(X,T,V,Rho,C,paso,i,Xr,gamma=5/3)##no usar ese punto y cambiarlo por el borde
            if Rhor_plus <0:
                continue

            X_plus = np.append(X_plus, Xr)
            V_plus = np.append(V_plus, Vr)
            T_plus = np.append(T_plus, Tr)
            Rho_plus = np.append(Rho_plus, Rhor)
            C_plus = np.append(C_plus, Cr)

        Tr, Xr, Vr, Rhor, Cr = bordes(X,T,V,Rho,C,paso,-1,1,gamma=5/3)
        X_plus = np.append(X_plus, Xr)
        V_plus = np.append(V_plus, Vr)
        T_plus = np.append(T_plus, Tr)
        Rho_plus = np.append(Rho_plus, Rhor)
        C_plus = np.append(C_plus, Cr)

        #if np.array_equal(X_plus , np.sort(X_plus))==False:
            #print "el sort no es igual"
        tipo = [('X', float), ('V', float), ('T', float), ('Rho', float), ('C', float)]
        values =[]
        for k in range(len(X_plus)):
            values.append((X_plus[k], V_plus[k], T_plus[k], Rho_plus[k], C_plus[k] ))
        ToSort = np.array(values, dtype=tipo)
        np.sort(ToSort, order='X')

        for k in range(len(X_plus)):
            X_plus[k]=ToSort[k][0]
            V_plus[k]=ToSort[k][1]
            T_plus[k]=ToSort[k][2]
            Rho_plus[k]=ToSort[k][3]
            C_plus[k]=ToSort[k][4]
        Pasos_x = len(X_plus)
        return X_plus, V_plus, T_plus, Rho_plus, C_plus, Pasos_x
    elif paso=="impar":
        #print "entre impar"
        for i in range(1,Pasos_x):
            Tr = (X[i]-X[i-1]-(V[i]-C[i])*T[i] + (V[i-1]+C[i-1])*T[i-1])  / (V[i-1]+C[i-1]-V[i]+C[i])
            Xr = (V[i-1]+C[i-1])*(Tr-T[i-1]) + X[i-1]
            Vr = (-Rho[i]+Rho[i-1]+ (V[i]*Rho[i]/C[i]) + (V[i-1]*Rho[i-1]/C[i-1]) ) / (Rho[i]/C[i] + Rho[i-1]/C[i-1])
            Rhor = (Vr - V[i])*Rho[i]/C[i] + Rho[i]
            if X[i]<=1/3:
                Cr=np.sqrt(gamma*4*(Rhor)**(gamma-1))
            elif X[i]>1/3:
                Cr=np.sqrt(gamma*1*(Rhor)**(gamma-1))
            for j in range(5):
                #print "que sueño",(V[i-1]+Vr+C[i-1]+Cr-V[i]-Vr+C[i]+Cr)
                Tr_plus = 2*( X[i]-X[i-1]+(V[i-1]+Vr+C[i-1]+Cr)*(T[i-1] / 2) - (V[i]+Vr-C[i]-Cr)*(T[i] / 2) ) / (V[i-1]+Vr+C[i-1]+Cr-V[i]-Vr+C[i]+Cr)
                Xr_plus = (1/2)*(V[i-1]+Vr+C[i-1]+Cr)*(Tr_plus-T[i-1]) +X[i-1]
                Vr_plus =( Rho[i-1] - Rho[i] + V[i]*(Rho[i]+Rhor)/(C[i]+Cr) + V[i-1]*(Rho[i-1]+Rhor)/(C[i-1]+Cr) )/ ( ((Rho[i]+Rhor)/(C[i]+Cr)) + ((Rho[i-1]+Rhor)/(C[i-1]+Cr)) )
                Rhor_plus = (Vr_plus-V[i])*(Rho[i]+Rhor)/(C[i]+Cr) +Rho[i]
                if Rhor_plus <0:
                    continue
                if Xr_plus<=1/3:
                    Cr_plus=np.sqrt(gamma*4*np.absolute(Rhor_plus)**(gamma-1))
                elif Xr_plus>1/3:
                    Cr_plus=np.sqrt(gamma*1*np.absolute(Rhor_plus)**(gamma-1))
                #if np.isnan(Tr_plus)==True or np.isnan(Xr_plus)==True or np.isnan(Vr_plus)==True or np.isnan(Rhor_plus)==True or np.isnan(Cr)==True :
                    #print "se hizo cero", Tr,Xr,Vr,Rhor,Cr
                if np.absolute(Tr_plus-Tr)/Tr <10**(-10):
                    Tr, Xr, Vr, Rhor, Cr = Tr_plus, Xr_plus, Vr_plus, Rhor_plus, Cr_plus
                    #print "quebraré", j
                    break
                Tr, Xr, Vr, Rhor, Cr = Tr_plus, Xr_plus, Vr_plus, Rhor_plus, Cr_plus
            if Xr<=0 or Xr>=1:
                continue
                #Tr, Xr, Vr, Rhor, Cr = bordes(X,T,V,Rho,C,paso,i,Xr,gamma=5/3)
            if Rhor_plus <0:
                continue
            X_plus = np.append(X_plus, Xr)
            V_plus = np.append(V_plus, Vr)
            T_plus = np.append(T_plus, Tr)
            Rho_plus = np.append(Rho_plus, Rhor)
            C_plus = np.append(C_plus, Cr)

            #print i
        #print "sali impar"
        Tr, Xr, Vr, Rhor, Cr = bordes(X,T,V,Rho,C,paso,-1,1,gamma=5/3)
        X_plus = np.append(X_plus, Xr)
        V_plus = np.append(V_plus, Vr)
        T_plus = np.append(T_plus, Tr)
        Rho_plus = np.append(Rho_plus, Rhor)
        C_plus = np.append(C_plus, Cr)

        #if np.array_equal(X_plus , np.sort(X_plus))==False:
            #print "el sort no es igual"
        tipo = [('X', float), ('V', float), ('T', float), ('Rho', float), ('C', float)]
        values =[]
        for k in range(len(X_plus)):
            values.append((X_plus[k], V_plus[k], T_plus[k], Rho_plus[k], C_plus[k] ))
        ToSort = np.array(values, dtype=tipo)
        np.sort(ToSort, order='X')

        for k in range(len(X_plus)):
            X_plus[k]=ToSort[k][0]
            V_plus[k]=ToSort[k][1]
            T_plus[k]=ToSort[k][2]
            Rho_plus[k]=ToSort[k][3]
            C_plus[k]=ToSort[k][4]
        Pasos_x = len(X_plus)
        return X_plus, V_plus, T_plus, Rho_plus, C_plus, Pasos_x
#Funcion para bordes
def bordes(X,T,V,Rho,C,paso,i,Xr,gamma=5/3):
    if paso == "par" and Xr<=0:
        Xr=0
        Vr=0
        Tr = (Xr-X[i+1])/(V[i+1]-C[i+1]) + T[i+1]
        Rhor = (Vr-V[i+1])*Rho[i+1]/C[i+1] + Rho[i+1]
        Cr = np.sqrt(np.absolute(gamma*4*Rhor**(gamma-1)))
        for j in range(5):
            Tr_plus = 2.0*(Xr-X[i+1])/(V[i+1]+Vr-C[i+1]-Cr) + T[i+1]
            Rhor_plus = (Vr-V[i+1])*(Rho[i+1]+Rhor)/(C[i+1]+Cr) + Rho[i+1]
            Cr_plus = np.sqrt(np.absolute(gamma*4*Rhor_plus**(gamma-1)))
            if Tr!= 0 and np.absolute(Tr_plus-Tr)/Tr <10**(-10):
                Tr, Rhor, Cr = Tr_plus, Rhor_plus, Cr_plus
        return Tr, Xr, Vr, Rhor, Cr
    if paso == "impar" and Xr<=0:
        Xr=0
        Vr=0
        Tr = (Xr-X[i])/(V[i]-C[i]) + T[i]
        Rhor = (Vr-V[i])*Rho[i]/C[i] + Rho[i]
        Cr = np.sqrt(np.absolute(gamma*4*Rhor**(gamma-1)))
        for j in range(5):
            Tr_plus = 2.0*(Xr-X[i])/(V[i]+Vr-C[i]-Cr) + T[i]
            Rhor_plus = (Vr-V[i])*(Rho[i]+Rhor)/(C[i]+Cr) + Rho[i]
            Cr_plus = np.sqrt(np.absolute(gamma*4*Rhor_plus**(gamma-1)))
            if Tr!= 0 and np.absolute(Tr_plus-Tr)/Tr <10**(-10):
                Tr, Rhor, Cr = Tr_plus, Rhor_plus, Cr_plus
        return Tr, Xr, Vr, Rhor, Cr
    if paso == "par" and Xr>=1:
        Xr=1
        Vr=0
        Tr = (Xr-X[i])/(V[i]+C[i]) + T[i]
        Rhor = -1*(Vr-V[i])*Rho[i]/C[i] + Rho[i]
        Cr = np.sqrt(np.absolute(gamma*1*Rhor**(gamma-1)))
        for j in range(5):
            Tr_plus = 2.0*(Xr-X[i])/(V[i]+Vr+C[i]+Cr) + T[i]
            Rhor_plus = -1*(Vr-V[i])*(Rho[i]+Rhor)/(C[i]+Cr) + Rho[i]
            Cr_plus = np.sqrt(np.absolute(gamma*1*Rhor_plus**(gamma-1)))
            if Tr!= 0 and np.absolute(Tr_plus-Tr)/Tr <10**(-10):
                Tr, Rhor, Cr = Tr_plus, Rhor_plus, Cr_plus
        return Tr, Xr, Vr, Rhor, Cr
    if paso == "impar" and Xr>=1:
        Xr=1
        Vr=0
        Tr = (Xr-X[i-1])/(V[i-1]+C[i-1]) + T[i-1]
        Rhor = -1*(Vr-V[i-1])*Rho[i-1]/C[i-1] + Rho[i-1]
        Cr = np.sqrt(np.absolute(gamma*1*Rhor**(gamma-1)))
        for j in range(5):
            Tr_plus = 2.0*(Xr-X[i-1])/(V[i-1]+Vr+C[i-1]+Cr) + T[i-1]
            Rhor_plus = -1*(Vr-V[i-1])*(Rho[i-1]+Rhor)/(C[i-1]+Cr) + Rho[i-1]
            Cr_plus = np.sqrt(np.absolute(gamma*1*Rhor_plus**(gamma-1)))
            if Tr!= 0 and np.absolute(Tr_plus-Tr)/Tr <10**(-10):
                Tr, Rhor, Cr = Tr_plus, Rhor_plus, Cr_plus
        return Tr, Xr, Vr, Rhor, Cr

#Funcion para actualizar variables avanzando en tiempo
def avanzar(X,V,T,Rho,C,paso,Pasos_x):
    #print "avance"
    Siguientes_datos = intermedios(X,T,V,Rho,C,paso,Pasos_x) #atento por si no actualiza bien
    X = Siguientes_datos[0]
    V = Siguientes_datos[1]
    T = Siguientes_datos[2]
    Rho = Siguientes_datos[3]
    C = Siguientes_datos[4]
    Pasos_x = Siguientes_datos[5]
    Tiempo = T[5000]
    return X, V, T, Rho, C, Tiempo, Pasos_x
#Funcion que retorna los valores
def valores_actuales():
    pass

# Main
# Setup
Pasos_x = 10001
gamma = 5/3

X = np.linspace(0,1,Pasos_x)
T = np.zeros(Pasos_x)
V = np.zeros(Pasos_x)
Rho = np.zeros(Pasos_x)
C = np.zeros(Pasos_x)

Rho, C = inicializa(gamma)
print Rho
tiempo=0
paso="par"
Contador=0

Data = []

while tiempo<1.2:
    X, V, T, Rho, C, tiempo, Pasos_x = avanzar(X,V,T,Rho,C,paso,Pasos_x)
    if paso=="par":
        paso="impar"
    elif paso=="impar":
        paso="par"
    if Contador==100 or tiempo >1.2:
        print "tiempo", tiempo, "Pasos_x",Pasos_x, Contador
        X=np.ndarray.tolist(X)
        T=np.ndarray.tolist(T)
        V=np.ndarray.tolist(V)
        Rho=np.ndarray.tolist(Rho)
        Data.append( (X,T,V,Rho))
        Contador=0
    #Pasos_x +=1
    Contador +=1

np.save('Datos', Data)
