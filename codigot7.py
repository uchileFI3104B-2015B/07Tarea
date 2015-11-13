from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

'''Rellena con las condiciones iniciales para rho'''


def inicializa_rho(rho, N_stephs, h):
    for i in range(N_steps):
        x = i*h
        if x > 0.1:
            rho[i] = 1
        else:
            rho[i] = 1 + 0.0375 * (1 + np.cos(10 * np.pi * x))


def inicializa_v(v, N_steps):

    '''velocidad nula pata t=0 en todo el espacio x
    incluyendo condiciones de borde con x=0 y x=1'''

    for i in range(N_steps):
        v[i] = 0


gamma = 5/3.
'''Dependencia con A para la presion. Presion asociada a cualquier densidad'''


def presion(x, rho):
    for i in range(0, N_steps):
        if x[i] < 0.3:
            p[i] = (rho[i]) ** gamma
        else:
            p[i] = 4 * (rho[i]) ** gamma


N_steps = 100
N_pasos_temporales = 100
t = np.zeros(N_steps)

''' x originales de la primera recta'''

x = np.linspace(0, 1, N_steps)
x[0] = 0
x[-1] = 1
x_inicial = 0
x_final = 1
rho = np.zeros(N_steps)
# h = (x_final-x_inicial) / (N_steps - 1)
h = 1/9.


c = np.zeros(N_steps)
v = np.zeros(N_steps)
p = np.zeros(N_steps)

''' Se definen arreglos dentro de funcion calculos_r en el for '''

x_r = np.zeros(N_steps)
t_r = np.zeros(N_steps)
v_r = np.zeros(N_steps)
rho_r = np.zeros(N_steps)
c_r = np.zeros(N_steps)
p_r = np.zeros(N_steps)


'''inicializo rho y v para t=0 '''

inicializa_rho(rho, N_steps, h)
inicializa_v(v, N_steps)
presion(x, rho)
# funciona bien, entrega las presiones para x entre 0 y 1 para t=0
# con las densidades rho de cada x respectivo


def calculos_r(p, x, t, v, c, rho):

        for i in range(0, N_steps):

            c2 = (gamma) * p / rho
            c = np.sqrt(c2)

            m = x[i] + (v[i] - c[i]) * (t[i-1] - t[i])
            - (v[i] - c[i]) * x[i-1] / (v[i-1] + c[i-1])
            b = 1 + (c[i] - v[i]) / (v[i-1] + c[i-1])
            x_r[i] = m / b

            d = t[i] * (v[i] - c[i]) + x[i-1] - x[i]
            - t[i-1] * (v[i-1] + c[i-1])
            a = v[i] - c[i] - v[i-1] - c[i-1]
            t_r[i] = d / a

            e = (c[i] * rho[i-1] / rho[i]) + c[i-1]
            f = rho[i-1] * (v[i-1] + c[i-1] + c[i]-v[i])
            rho_r[i] = f / e

            g = c[i] * (rho_r[i]) + rho[i] * (v[i] - c[i])
            h = rho[i]
            v_r[i] = g / h


v_solucion = np.zeros((N_pasos_temporales, N_steps))
rho_solucion = np.zeros((N_pasos_temporales, N_steps))
v_solucion[0, :] = v.copy()
# caso base : obtengo rho para t=0 para el arreglo de x entre 0 y 1
rho_solucion[0, :] = rho.copy()
for i in range(1, N_pasos_temporales):
    presion(x, rho)
    calculos_r(p, x, t, v, c, rho)
    rho = rho_r.copy()
    v = v_r.copy()
    rho_solucion[i, :] = rho.copy()
    v_solucion[i, :] = v.copy()


fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

ax.plot(x, rho_solucion[0, :], 'g')
ax.plot(x, 1 + v_solucion[0, :], 'r')
ax.set_ylim(0.98, 1.09)
ax.set_title('Grafico rho(x) y 1+v(x) v/s x en t=0')
ax.set_xlabel('x [unidades arbitrarias]')
ax.set_ylabel('rho(x)[verde] y 1+v(x)[rojo] [unidades arbitrarias]')
plt.savefig('figura0rhovsvx1.png')


fig = plt.figure(2)
fig.clf()
ax2 = fig.add_subplot(111)

ax2.plot(x, rho_solucion[0.2, :], 'g')
ax2.plot(x, 1 + v_solucion[0.2, :], 'r')
ax2.set_ylim(0.97, 1.09)
ax2.set_title('Grafico rho(x) y 1+v(x) v/s x t=0.2')
ax2.set_xlabel('x [unidades arbitrarias]')
ax2.set_ylabel('rho(x)[verde] y 1+v(x)[rojo] [unidades arbitrarias]')
plt.savefig('figura02rhovsvx1.png')

fig = plt.figure(3)
fig.clf()
ax3 = fig.add_subplot(111)
ax3.plot(x, rho_solucion[0.4, :], 'g')
ax3.plot(x, 1 + v_solucion[0.4, :], 'r')
ax3.set_ylim(-10, 10)
ax3.set_title('Grafico rho(x) y 1+v(x) v/s x en t=0.4')
ax3.set_xlabel('x [unidades arbitrarias]')
ax3.set_ylabel('rho(x)[verde] y 1+v(x)[rojo] [unidades arbitrarias]')
plt.savefig('figura04rhovs1vx.png')


fig = plt.figure(4)
fig.clf()
ax4 = fig.add_subplot(111)
ax4.plot(x, rho_solucion[0.6, :], 'g')
ax4.plot(x, 1 + v_solucion[0.6, :], 'r')
ax4.set_ylim(-10, 10)
ax4.set_title('Grafico rho(x) y 1+v(x) v/s x en t=0.6')
ax4.set_xlabel('x [unidades arbitrarias]')
ax4.set_ylabel('rho(x)[verde] y 1+v(x)[rojo] [unidades arbitrarias]')
plt.savefig('figura06rhovs1vx.png')

fig = plt.figure(5)
fig.clf()
ax5 = fig.add_subplot(111)
ax5.plot(x, rho_solucion[0.8, :], 'g')
ax5.plot(x, 1 + v_solucion[0.8, :], 'r')
ax5.set_ylim(-10, 10)
ax5.set_title('Grafico rho(x) y 1+v(x) v/s x en t=0.8')
ax5.set_xlabel('x [unidades arbitrarias]')
ax5.set_ylabel('rho(x)[verde] y 1+v(x)[rojo] [unidades arbitrarias]')
plt.savefig('figura08rhovs1vx.png')


fig = plt.figure(6)
fig.clf()
ax6 = fig.add_subplot(111)
ax6.plot(x, rho_solucion[1, :], 'g')
ax6.plot(x, 1 + v_solucion[1, :], 'r')
ax6.set_ylim(-10, 10)
ax6.set_title('Grafico rho(x) y 1+v(x) v/s x en t=1.0')
ax6.set_xlabel('x [unidades arbitrarias]')
ax6.set_ylabel('rho(x)[verde] y 1+v(x)[rojo] [unidades arbitrarias]')
plt.savefig('figura1rhovs1vx.png')

dt = 1/N_pasos_temporales

fig7 = plt.figure(7)
fig7.clf()
ax7 = fig7.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax7.pcolormesh(X, Y, v_solucion)
fig7.colorbar(ax4.pcolormesh(X, Y, v_solucion))
ax7.set_title('Grafico v(x) en plano(x,t) ')
ax7.set_xlabel('x[unidades arbitrarias]')
ax7.set_ylabel('t[seg]')
plt.savefig('figuravplano.png')


fig8 = plt.figure(8)
fig8.clf()
ax8 = fig8.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax8.pcolormesh(X, Y, rho_solucion)
fig8.colorbar(ax8.pcolormesh(X, Y, rho_solucion))
ax8.set_title('Grafico rho(x) en plano(x,t) ')
ax8.set_xlabel('x[unidades arbitrarias]')
ax8.set_ylabel('t[seg]')
plt.savefig('figurarhoplano.png')

plt.show()
plt.draw()
