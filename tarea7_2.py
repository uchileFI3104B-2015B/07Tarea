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

    cp=A*gamma*rhop**(gamma-1.0)
    cq=A*gamma*rhoq**(gamma-1.0)

    tr=(xq-xp+tp*(vp+cp)-tq*(vq-cq))/(vp+cp-vq+cq)
    xr=xp+(vp+cp)*(tr-tp)
    vr=(vp*rhop/cp + vq*rhoq/cq +rhop-rhoq)/(rhop/cp + rhoq/cq)
    rhor=rhop*(1-(vr-vp)/cp)
    cr=A*gamma*rhor**(gamma-1)

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

m=10    

A=n.ones(10000)

for i in range(3334):
    A[i]=4.0
    
gamma=5.0/3
x=n.linspace(0,0.9999,10000)
rho=n.ones(10000)
t=n.zeros(10000)
v=n.zeros(10000)


for i in range(1001):
    rho[i]=1.0+0.0977*(1.0+cos(10*pi*x[i]))

c=A*gamma*n.power(rho, (gamma-1))    
    

x1=n.zeros((2001,10000))

rho1=n.zeros((2001,10000))

t1=n.zeros((2001,10000))

v1=n.zeros((2001,10000))

c1=n.zeros((2001,10000))


for k in range(1000):

    for i in range(9999):
        b=avanza(x[i],x[i+1],t[i],t[i+1],v[i],v[i+1],rho[i],rho[i+1],m,A[i],gamma)
        t1[2*k+1,i]=b[0]
        x1[2*k+1,i]=b[1]
        v1[2*k+1,i]=b[2]
        rho1[2*k+1,i]=b[3]
        c1[2*k+1,i]=b[4]

    x1[2*k+1,9999]=1.0
    v1[2*k+1,9999]=0.0
    t1[2*k+1,9999]=(1.0-x[9999])/(v[9999]+c[9999])+t[9999]
    rho1[2*k+1,9999]=rho[9999]*(1.0 + v[9999]/c[9999])
    c1[2*k+1,9999]=A[9999]*gamma*n.power(rho1[2*k+1,9999],gamma-1.0)

    for j in range(9999):
        a=avanza(x1[2*k+1,j],x1[2*k+1,j+1],t1[2*k+1,j],t1[2*k+1,j+1],v1[2*k+1,j],v1[2*k+1,j+1],rho1[2*k+1,j],rho1[2*k+1,j+1],m,A[j],gamma)
        t1[2*k+2,j+1]=a[0]
        x1[2*k+2,j+1]=a[1]
        v1[2*k+2,j+1]=a[2]
        rho1[2*k+2,j+1]=a[3]
        c1[2*k+2,j+1]=a[4]    

    x1[2*k+2,0]=0.0
    v1[2*k+2,0]=0.0
    t1[2*k+2,0]=(-x1[2*k+1,0])/(v1[2*k+1,0]-c1[2*k+1,0])+t1[2*k+1,0]
    rho1[2*k+2,0]=rho1[2*k+1,0]*(1.0 - v1[2*k+1,0]/c1[2*k+1,0])
    c1[2*k+2,0]=A[0]*gamma*n.power(rho1[2*k+2,0],gamma-1.0)

    x=x1[2*k+2,:]
    t=t1[2*k+2,:]
    v=v1[2*k+2,:]
    rho=rho1[2*k+2,:]

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


V=v1+1.0



p.plot(x1[10,:],rho1[10,:], label='Densidad')
p.plot(x1[10,:],V[10,:], label='Velocidad +1')
p.xlabel('x')
p.title('Densidad y velocidad')
p.show()

p.plot(x1[500,:],rho1[500,:], label='Densidad')
p.plot(x1[500,:],V[500,:], label='Velocidad +1')
p.xlabel('x')
p.title('Densidad y velocidad')
p.show()

p.plot(x1[1200,:],rho1[1200,:], label='Densidad')
p.plot(x1[1200,:],V[1200,:], label='Velocidad +1')
p.xlabel('x')
p.title('Densidad y velocidad')
p.show()

p.plot(x1[1998,:],rho1[1998,:], label='Densidad')
p.plot(x1[1998,:],V[1998,:], label='Velocidad +1')
p.xlabel('x')
p.title('Densidad y velocidad')
p.show()

x2=n.zeros((201,100))  #reducir datos para grafico 3D
t2=n.zeros((201,100))
v2=n.zeros((201,100))

for i in range(201):   #selecciona 1 de cada 10 en t y 1 de 100 en x
    for j in range(100):
        x2[i,j]=x1[10*i,100*j]
        t2[i,j]=t1[10*i,100*j]
	v2[i,j]=v1[10*i,100*j]
            


grafica3D(x2,t2,v2)
grafica3D(x2,t2,rho2)

graficaIntensidad(v2)
graficaIntensidad(rho2)

