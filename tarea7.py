from __future__ import division
from math import *
from scipy import integrate as int
from scipy import optimize
from scipy.integrate import ode
import pyfits #modulo para leer archivos fits
import matplotlib.pyplot as p #modulo para graficar
import numpy as n #este modulo es para trabajar con matrices como en matlab
import matplotlib as mp
from mpl_toolkits.mplot3d import Axes3D

from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis, title, show
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def avanza(xp,xq,tp,tq,vp,vq,rhop,rhoq,m,A,gamma):

    '''metodo que permite calcular, dados dos puntos P y Q, los valores
    de densidad, tiempo, x, velocidad y velocidad del sonido del punto R,
    correspondiente a la interseccion de las caracteristicas.
    Se ingresan los valores x, t, v y rho de P y Q, junto con los parametros
    de la ec de estado A y gamma; y m es el orden de precisión de los datos
    '''

    cp=A*gamma*rhop**(gamma-1.0)  #vel del sonido
    cq=A*gamma*rhoq**(gamma-1.0)
    
    #iteraciones a primer orden
    tr=(xq-xp+tp*(vp+cp)-tq*(vq-cq))/(vp+cp-vq+cq)
    xr=xp+(vp+cp)*(tr-tp)
    vr=(vp*rhop/cp + vq*rhoq/cq +rhop-rhoq)/(rhop/cp + rhoq/cq)
    rhor=rhop*(1-(vr-vp)/cp)
    cr=A*gamma*rhor**(gamma-1)

    
    #iteraciones a orden superior
    for i in range(m):
        tn= (2*(xq-xp)+tp*(vp+vr+cp+cr)-tq*(vq+vr-cq-cr))/(vp+cp-vq+cq+2*cr)
        xn= xp+ 0.5*(vp+vr+cp+cr)*(tn-tp)
        vn=(vp*(rhop+rhor)/(cp+cr)+ vq*(rhoq+rhor)/(cq+cr) +rhop-rhoq)/((rhop+rhor)/(cp+cr)+(rhoq+rhor)/(cq+cr))
        rhon= rhop- (vn-vp)*(rhop+rhor)/(cp+cr)
        cn=A*gamma*rhon**(gamma-1)
        tr=tn
        xr=xn
        vr=vn
        rhor=rhon
        cr=cn
        
    return [tn,xn,vn,rhon,cn]

m=10   #orden de precision de los calculos en metodo 'avanza'

#def de los valores de A
A=n.ones(10000)

for i in range(3334):
    A[i]=4.0
    
gamma=5.0/3
x=n.linspace(0,0.9999,10000)   #division del intervalo [0,1]
rho=n.ones(10000)    
t=n.zeros(10000)
v=n.zeros(10000)

#densidad inicial
for i in range(1001):
    rho[i]=1.0+0.0977*(1.0+cos(10*pi*x[i]))

c=A*gamma*n.power(rho, (gamma-1))     #vel del sonido
    

x1=n.ones((2,10000))

rho1=n.ones((2,10000))

t1=n.ones((2,10000))

v1=n.ones((2,10000))

c1=n.ones((2,10000))

xf=n.ones((1000,10000))
rhof=n.ones((1000,10000))
tf=n.ones((1000,10000))
vf=n.ones((1000,10000))
cf=n.ones((1000,10000))

for l in range(1000):  #1000 es el numero de curvas que se obtendran
    


    for k in range(50):  #se itera 50 veces este for para obtener 100 pasos

        for i in range(9999):  #se avanza sobre el intervalo [0,1]
            b=avanza(x[i],x[i+1],t[i],t[i+1],v[i],v[i+1],rho[i],rho[i+1],m,A[i],gamma)
            t1[0,i]=b[0]
            x1[0,i]=b[1]
            v1[0,i]=b[2]
            rho1[0,i]=b[3]
            c1[0,i]=b[4]

        #calculo del borde derecho
        x1[0,9999]=1.0
        v1[0,9999]=0.0
        t1[0,9999]=(1.0-x[9999])/(v[9999]+c[9999])+t[9999]
        rho1[0,9999]=rho[9999]*(1.0 + v[9999]/c[9999])
        c1[0,9999]=A[9999]*gamma*n.power(rho1[0,9999],gamma-1.0)

        for j in range(9999):
            a=avanza(x1[0,j],x1[0,j+1],t1[0,j],t1[0,j+1],v1[0,j],v1[0,j+1],rho1[0,j],rho1[0,j+1],m,A[j],gamma)
            t1[1,j+1]=a[0]
            x1[1,j+1]=a[1]
            v1[1,j+1]=a[2]
            rho1[1,j+1]=a[3]
            c1[1,j+1]=a[4]    

        #calculo del borde izquierdo
        x1[1,0]=0.0
        v1[1,0]=0.0
        t1[1,0]=(-x1[0,0])/(v1[0,0]-c1[0,0])+t1[0,0]
        rho1[1,0]=rho1[0,0]*(1.0 - v1[0,0]/c1[0,0])
        c1[1,0]=A[0]*gamma*n.power(rho1[1,0],gamma-1.0)

        x=x1[1,:]
        t=t1[1,:]
        v=v1[1,:]
        rho=rho1[1,:]

    #se llenan estas matrices con el paso 100 de cada for
    xf[l,:]=x1[1,:]
    tf[l,:]=t1[1,:]
    vf[l,:]=v1[1,:]
    rhof[l,:]=rho1[1,:]


def graficaIntensidad(Z):

    'grafica curvas de nivel'

    im = imshow(Z,cmap=cm.RdBu)

    cset = contour(Z,n.arange(-2.0,2.0,0.1),linewidths=2,cmap=cm.Set2)

    clabel(cset,inline=True,fmt='%1.1f',fontsize=10)

    colorbar(im)

    show()

def grafica3D(X,Y,Z):

    'grafica funciones de 2 variables'

    fig = p.figure()

    fig.clf()

    ax = fig.add_subplot(111,projection='3d')

    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,cmap=cm.RdBu,linewidth=0.15)        

    fig.colorbar(surf, shrink=0.5, aspect=5)

    p.show()


V=vf+1.0

grafica3D(xf,tf,vf)
grafica3D(xf,tf,rhof)

graficaIntensidad(v1)
graficaIntensidad(rho1)

p.plot(x1[10,:],rho1[10,:], label='Densidad')
p.plot(x1[10,:],V[10,:], label='Velocidad +1')
p.xlabel('x')
p.title('Densidad y velocidad')
p.show()

