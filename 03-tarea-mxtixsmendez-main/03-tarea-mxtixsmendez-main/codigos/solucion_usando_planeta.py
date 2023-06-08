from planeta import Planeta
import matplotlib.pyplot as plt
import numpy as np
vy = 0.2
condicion_inicial = [10, 0, 0, vy]

mercurio = Planeta(condicion_inicial)
# Método RK4
n = 5*10000
dt = 0.01
x = np.zeros(n)
y = np.zeros(n)
vx = np.zeros(n)
vy = np.zeros(n)
vel_ang = np.zeros(n)  # velocidad angular
t = np.zeros(n)
energia_total_rk4 = np.zeros(n)
for i in range(n):
    x[i], y[i], vx[i], vy[i] = mercurio.y_actual
    t[i] = mercurio.t_actual
    energia_total_rk4[i] = mercurio.energia_total()
    radio = np.sqrt(x[i]**2+y[i]**2)
    vel_ang[i] = radio*(vy[i]*x[i]-vx[i]*y[i])/(x[i]**2+y[i]**2)
    mercurio.avanza_rk4(dt)
afelios = 0
for i in range(1, n-1):
    pendiente_prev = (vel_ang[i-1]-vel_ang[i])/(t[i-1]-t[i])
    pendiente_sgt = (vel_ang[i+1]-vel_ang[i])/(t[i+1]-t[i])
    if pendiente_prev*pendiente_sgt < 0 and pendiente_prev < 0:
        afelios += 1
print('afelios en rk4: ', afelios)

# Velocidad angular
fig0 = plt.figure(0)
ax0 = fig0.add_subplot(111)
ax0.plot(t, vel_ang, label='Método RK4')
ax0.grid()
ax0.set_xlabel('t', fontsize=15)
ax0.set_ylabel('Velocidad angular $r\dot{\Theta}$', fontsize=15)
ax0.legend()
fig0.show()

# Órbita RK4
fig1 = plt.figure(1)
fig1.clf()
ax1 = fig1.add_subplot(111)
ax1.plot(x, y, label='Método RK4')
ax1.grid()
ax1.set_xlabel('Eje X', fontsize=15)
ax1.set_ylabel('Eje Y', fontsize=15)
ax1.legend()
fig1.savefig('orbita rk4.pdf')

# Energía total RK4
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.plot(t, energia_total_rk4, label='Método RK4')
ax2.grid()
ax2.set_xlabel('tiempo t', fontsize=15)
ax2.set_ylabel('Energía total', fontsize=15)
ax2.legend()
fig2.savefig('energia total rk4.pdf')


# Método Verlet
x = np.zeros(n)
y = np.zeros(n)
vx = np.zeros(n)
vy = np.zeros(n)
vel_ang = np.zeros(n)  # velocidad angular
t = np.zeros(n)
energia_total_verlet = np.zeros(n)
mercurio.avanza_rk4(dt)
for i in range(0, n):
    mercurio.avanza_verlet(dt)
    x[i], y[i], vx[i], vy[i] = mercurio.y_actual
    t[i] = mercurio.t_actual
    radio = np.sqrt(x[i] ** 2 + y[i] ** 2)
    vel_ang[i] = radio*(vy[i]*x[i]-vx[i]*y[i])/(x[i]**2+y[i]**2)
    energia_total_verlet[i] = mercurio.energia_total()
afelios = 0
for i in range(1, n-1):
    pendiente_prev = (vel_ang[i-1]-vel_ang[i])/(t[i-1]-t[i])
    pendiente_sgt = (vel_ang[i+1]-vel_ang[i])/(t[i+1]-t[i])
    if pendiente_prev*pendiente_sgt < 0 and pendiente_prev < 0:
        afelios += 1
print('afelios en verlet: ', afelios)
# Órbita verlet
fig3 = plt.figure(3)
fig3.clf()
ax3 = fig3.add_subplot(111)
ax3.plot(x, y, label='Método Verlet')
ax3.grid()
ax3.set_xlabel('Eje X', fontsize=15)
ax3.set_ylabel('Eje Y', fontsize=15)
ax3.legend()
fig3.savefig('orbita verlet.pdf')

# Energía total Verlet
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)
ax4.plot(t, energia_total_verlet, label='Método Verlet')
ax4.grid()
ax4.set_xlabel('tiempo t', fontsize=15)
ax4.set_ylabel('Energía total', fontsize=15)
ax4.legend()
fig4.savefig('energia total verlet.pdf')
# Método Beeman
x = np.zeros(n)
y = np.zeros(n)
vx = np.zeros(n)
vy = np.zeros(n)
vel_ang = np.zeros(n)  # velocidad angular
t = np.zeros(n)
energia_total_beeman = np.zeros(n)
mercurio.avanza_rk4(dt)
for i in range(0, n):
    mercurio.avanza_beeman(dt)
    x[i], y[i], vx[i], vy[i] = mercurio.y_actual
    t[i] = mercurio.t_actual
    radio = np.sqrt(x[i] ** 2 + y[i] ** 2)
    vel_ang[i] = radio*(vy[i]*x[i]-vx[i]*y[i])/(x[i]**2+y[i]**2)
    energia_total_beeman[i] = mercurio.energia_total()
afelios = 0
for i in range(1, n-1):
    pendiente_prev = (vel_ang[i-1]-vel_ang[i])/(t[i-1]-t[i])
    pendiente_sgt = (vel_ang[i+1]-vel_ang[i])/(t[i+1]-t[i])
    if pendiente_prev*pendiente_sgt < 0 and pendiente_prev < 0:
        afelios += 1
print('afelios en beeman: ', afelios)
# Órbita Beeman
fig5 = plt.figure(5)
fig5.clf()
ax5 = fig5.add_subplot(111)
ax5.plot(x, y, label='Método Beeman')
ax5.grid()
ax5.set_xlabel('Eje X', fontsize=15)
ax5.set_ylabel('Eje Y', fontsize=15)
ax5.legend()
fig5.savefig('orbita beeman.pdf')


# Energía total Beeman
fig6 = plt.figure(6)
ax6 = fig6.add_subplot(111)
ax6.plot(t, energia_total_beeman, label='Método Beeman')
ax6.grid()
ax6.set_xlabel('tiempo t', fontsize=15)
ax6.set_ylabel('Energía total', fontsize=15)
ax6.legend()
fig6.savefig('energia total beeman.pdf')
# Segunda parte
alpha = 10**(-2.742)
mercurio = Planeta(condicion_inicial, alpha)
n = 5*18000
dt = 0.1
# Método Beeman
x = np.zeros(n)
y = np.zeros(n)
vx = np.zeros(n)
vy = np.zeros(n)
vel_ang = np.zeros(n)  # velocidad angular
t = np.zeros(n)
energia_total_beeman = np.zeros(n)
mercurio.avanza_rk4(dt)
for i in range(0, n):
    mercurio.avanza_beeman(dt)
    x[i], y[i], vx[i], vy[i] = mercurio.y_actual
    t[i] = mercurio.t_actual
    radio = np.sqrt(x[i] ** 2 + y[i] ** 2)
    vel_ang[i] = radio*(vy[i]*x[i]-vx[i]*y[i])/(x[i]**2+y[i]**2)
    energia_total_beeman[i] = mercurio.energia_total()
afelios = 0
afelios_lista = list()
afelios_lista_tiempo = list()
for i in range(1, n-1):
    pendiente_prev = (vel_ang[i-1]-vel_ang[i])/(t[i-1]-t[i])
    pendiente_sgt = (vel_ang[i+1]-vel_ang[i])/(t[i+1]-t[i])
    if pendiente_prev*pendiente_sgt < 0 and pendiente_prev < 0:
        afelios += 1
        angulo = np.arctan(y[i]/x[i])
        afelios_lista.append(angulo)
        afelios_lista_tiempo.append(t[i])
print('afelios en beeman: ', afelios)
# Velocidad angular de precesión
fig7 = plt.figure(7)
fig7.clf()
ax7 = fig7.add_subplot(111)
ax7.plot(afelios_lista_tiempo, afelios_lista, '-o', label='Método Beeman')
ax7.grid()
ax7.set_xlabel('tiempo t', fontsize=15)
ax7.set_ylabel('Ángulo $\Theta$')
ax7.legend()
fig7.savefig('velocidad angular de precesion.pdf')
m1 = (afelios_lista[2]-afelios_lista[1])
m2 = (afelios_lista_tiempo[2]-afelios_lista_tiempo[1])
m = m1/m2
print('velocidad de precesion (pendiente)', m)
# Órbita Beeman
fig8 = plt.figure(8)
fig8.clf()
ax8 = fig8.add_subplot(111)
ax8.plot(x, y, label='Método Beeman')
ax8.grid()
ax8.set_xlabel('Eje X', fontsize=15)
ax8.set_ylabel('Eje Y', fontsize=15)
ax8.legend()
fig8.savefig('orbitales beeman con alpha no cero.pdf')
# Energía total Beeman
fig9 = plt.figure(9)
ax9 = fig9.add_subplot(111)
ax9.plot(t, energia_total_beeman, label='Método Beeman')
ax9.grid()
ax9.set_xlabel('tiempo t', fontsize=15)
ax9.set_ylabel('Energía total', fontsize=15)
ax9.legend()
fig9.savefig('energia total beeman con alpha no cero.pdf')