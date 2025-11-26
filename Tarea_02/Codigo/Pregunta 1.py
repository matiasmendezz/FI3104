import numpy as np
import matplotlib.pyplot as plt

# constantes
g = 9.8  # en m/s^2, graveda
l = 5.427  # en m, largo del péndulo
gamma = 2.427  # en s^-1, coeficiente de fricción
phi_0 = 0.0001  # en rad, ángulo en condición inicial
phi_punto_0 = 0  # en rad/s^2, velocidad angular en condición inicial


# Pregunta 1
# Parte 1
# constante solucion analitica
w = np.sqrt((g / l) - (gamma / 2) ** 2)
A = np.sqrt(phi_0 ** 2 + ((phi_punto_0 + gamma/2*phi_0)/w) ** 2)
phi = np.arctan((phi_0 * w) / (phi_punto_0 + gamma / 2 * phi_0))


def sol_analitica(t):
    """
    Función analítica solución a la ecuación diferencial
    de un péndulo simple en un medio viscoso
    d^2 phi/t^2 + d phi/dt+ (g/l) phi = 0
    """
    global w, A, phi
    output = A * np.e ** (-gamma / 2 * t) * np.sin(w * t + phi)
    return output


# Solución numérica
def funcion_a_integrar(t, y):
    """
    Funcion a integrar númericamente
    Inputs:
    t : [float], tiempo
    y: [np.array], [phi, dphi/dt]
    """
    global gamma, l, g
    output = np.array([y[1], -gamma * y[1] - (g / l) * y[0]])
    return output


def K1(fun, dt, t_n, y_n):
    """
    Constante K1 del algoritmo de RK4
    """
    output = fun(t_n, y_n)
    return output


def K2(fun, dt, t_n, y_n):
    """
    Constante K2 del algoritmo de RK4
    """
    k1_n = K1(fun, dt, t_n, y_n)
    output = fun(t_n+1/2*dt, y_n+1/2*k1_n*dt)
    return output


def K3(fun, dt, t_n, y_n):
    """
    Constante K3 del algoritmo de RK4
    """
    k2_n = K2(fun, dt, t_n, y_n)
    output = fun(t_n+1/2*dt, y_n+1/2*k2_n*dt)
    return output


def K4(fun, dt, t_n, y_n):
    """
    Constante K4 del algoritmo de RK4
    """
    k3_n = K3(fun, dt, t_n, y_n)
    output = fun(t_n+dt, y_n+k3_n*dt)
    return output


def paso_rk4(fun, dt, t_n, y_n):
    """
    Calcula el paso y_n+1 a través del método de RK4
    y_n+1=y_n+1/6*dt*(k1+2*k2+2*k3+k4)
    """
    k1_n = K1(fun, dt, t_n, y_n)
    k2_n = K2(fun, dt, t_n, y_n)
    k3_n = K3(fun, dt, t_n, y_n)
    k4_n = K4(fun, dt, t_n, y_n)
    output = y_n + 1/6*dt*(k1_n+2*k2_n+2*k3_n+k4_n)
    return output


#  Gráfico solución analítica vs numérica
#  solución analítica
t_analitica = np.linspace(0, 10, 10000)
y_analitica = sol_analitica(t_analitica)

# Solución numérica
dt = 0.001
t_rk4 = np.arange(0, 10, dt)
y_rk4 = np.zeros((len(t_rk4), 2))  # y=[dphi/dt, phi]]
y_rk4[0] = [phi_0, phi_punto_0]  # condición inicial
for i in range(1, len(t_rk4)):
    y_rk4[i] = paso_rk4(funcion_a_integrar, dt, t_rk4[i-1], y_rk4[i-1])
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(211)
ax.plot(t_analitica, y_analitica, label='solución analítica')
ax.plot(t_rk4, y_rk4[:, 0], label='solución RK4')
ax.grid()
ax.set_xlabel('t', fontsize=15)
ax.set_ylabel('$\phi (t)$', fontsize=15)
ax.legend()
ax_1 = fig.add_subplot(212)
ax_1.grid()
error = np.fabs(y_rk4[:, 0]-y_analitica)
ax_1.plot(t_analitica, error, label='error')
ax_1.set_xlabel('t')
ax_1.set_ylabel('$|\phi(t)_{analitica}-\phi_{numerica}|$')
ax_1.legend()
fig.show()
fig.savefig('grafico phi sin perturbacion pregunta 1.pdf')
# Segunda parte
# constantes solucion analitica
phi_0 = np.pi/2.427  # Punto de oscilación
w = np.sqrt((g / l) - (gamma / 2) ** 2)
A = np.sqrt(phi_0 ** 2 + ((phi_punto_0 + gamma / 2 * phi_0) / w) ** 2)
phi = np.arctan((phi_0 * w) / (phi_punto_0 + gamma / 2 * phi_0))


def sol_analitica_2(t):
    """
    Función analítica solución a la ecuación diferencial
    de un péndulo simple en un medio viscoso
    d^2 phi/t^2 + d phi/dt+ (g/l) phi = 0
    """
    global w, A, phi
    output = A * np.e ** (-gamma / 2 * t) * np.sin(w * t + phi)
    return output


# Datos t,y para la solución analítica
y_analitica = sol_analitica_2(t_analitica)
dt = 0.001  # intervalo
# Datos para solución numérica
t_rk4 = np.arange(0, 10, dt)
y_rk4 = np.zeros((len(t_rk4), 2))  # y=[dphi/dt, phi]]
y_rk4[0] = [phi_0, phi_punto_0]  # condición inicial
for i in range(1, len(t_rk4)):
    y_rk4[i] = paso_rk4(funcion_a_integrar, dt, t_rk4[i-1], y_rk4[i-1])

# Grafico segunda parte
# Comparacion entre la solucion para pequeñas oscilaciones vs la numerica
# para un punto de oscilación pi/2.427
fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
ax2.plot(t_analitica, y_analitica, label='solución analítica')
ax2.plot(t_rk4, y_rk4[:, 0], label='solución RK4')
ax2.grid()
ax2.set_xlabel('t', fontsize=15)
ax2.set_ylabel('$\phi (t)$', fontsize=15)
ax2.legend()
# fig2.show()
fig2.savefig('grafico phi con perturbacion pregunta 1.pdf')

# Energia cinetica


def energia(t):
    """
    Calcula la energia total del pendulo en un momento t
    """
    global A, w, gamma, phi, g, l
    angulo_phi = sol_analitica_2(t)
    velocidad = A*np.e**(-0.5*gamma*t)*(w*np.cos(w*t+phi))
    energia_cinetica = 0.5*(l**2)*velocidad**2
    energia_potencial = -g*l*np.cos(angulo_phi)
    return energia_potencial+energia_cinetica


def energia_numerica(y):
    """
    Calcula la energia total para la solucion numerica de phi
    el array y contiene phi y su velocidad
    """
    phi = y[0]
    velocidad = y[1]
    energia_cinetica = 0.5*l**2*velocidad**2
    energia_potencial = -g*l*np.cos(phi)
    output = energia_cinetica+energia_potencial
    return output


energia_total = energia(t_analitica)
energia_numerica_total = np.zeros(len(t_rk4))
for i in range(len(t_rk4)):
    energia_numerica_total[i] = energia_numerica(y_rk4[i])
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.plot(t_analitica, energia_total, label='analítica')
ax3.plot(t_rk4, energia_numerica_total, label='RK4')
ax3.set_xlabel('t', fontsize=15)
ax3.set_ylabel('$E_{total}$', fontsize=15)
ax3.legend()
ax3.grid()
fig3.savefig('Grafico energia con perturbacion.pdf')
fig3.show()
