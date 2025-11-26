import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
"""
Pregunta 1
Los datos GLB.Ts+dSST.csv se pasan a un excel para leerlos. Este nuev excel
se llama datos.xlsx

Forma parametrica escogida: y = a_1 * (x-u_1)^3 + a_2
Se utiliza el método de mínimos cuadrados no lineales de scipy: curve_fit

Print entregados por el codigo:
Año promedio en el que la temperatura habra cambiado dos grados es 2054
El intervalo de confianza al 68% es: [2051.552755552305, 2057.151529799679]

"""
np.random.seed(1234)
# Lectura de Datos
data = pd.read_excel('datos.xlsx', engine='openpyxl')
años = np.arange(1880, 2020, 1)
temperaturas = data.values[2:]  # Arreglo con temperaturas por año
# Convertimos todos los numeros que estan a string a float
# Excepto a los '***'
for temperatura in temperaturas:
    for i in range(len(temperatura)):
        if not temperatura[i] == '***':
            temperatura[i] = float(temperatura[i])
j_d = np.zeros(len(años))
# Generamos el arreglo con las temperaturas J-D
for i in range(len(años)):
    año = años[i]
    j_d[i] = temperaturas[i][13]

# Funciones


def f_modelo(x, a_1, u_1, a_2):
    """
    Funcion modelo
    y = a_1 * (x-u_1)^3 + a_2
    """
    y = a_1*(x-u_1)**3 + a_2
    return y


def f_m_inversa(y, a_1, u_1, a_2):
    """
    Funcion modelo inversa:
    x = u_1+ ( (y-a_2) / a_1 )^(1/3)
    """
    x = u_1 + np.cbrt((y-a_2)/a_1)
    return x


# Grafico mejor fiteo
p0 = 1, 1950, 0.2  # parametros iniciales
param_opt, param_cov = curve_fit(f_modelo, años, j_d, p0=p0)
a_1, u_1, a_2 = param_opt
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(años, j_d, 'o', label='Temperaturas J-D')
ax.plot(años, f_modelo(años, a_1, u_1, a_2), label='curve_fit')
ax.set_xlabel('Años', fontsize=15)
ax.set_ylabel('Temperaturas [°C]', fontsize=15)
ax.set_title('Gráfico con fiteo para la temperatura J-D', fontsize=20)
ax.legend()
fig.savefig('Grafico mejor fiteo P1.pdf')
# Bootstrop
N = len(j_d)
Nboot = int(N * np.log10(N)**2)
dos_grados = np.zeros(Nboot)
for i in range(Nboot):
    s = np.random.randint(low=0, high=N, size=N)
    año_opt, año_cov = curve_fit(f_modelo, años[s], j_d[s], p0=p0)
    a_1, u_1, a_2 = año_opt
    dos_grados[i] = f_m_inversa(2, a_1, u_1, a_2)

# Grafico
fig_2 = plt.figure(2)
ax_2 = fig_2.add_subplot(111)
ax_2.hist(dos_grados, bins=30, edgecolor='black')
ax_2.axvline(np.mean(dos_grados), color='r')
ax_2.set_xlabel('Años', fontsize=15)
ax_2.set_ylabel('Frecuencia', fontsize=15)
ax_2.set_title('Histograma de año con 2 grados de cambio')
fig_2.savefig('Histograma año con 2 grados de cambio P1.pdf')


# Intervalo de confianza
# Definimos el intervalo que contiene el 68 % de los valores
# como intervalo de confianza
dos_grados = np.sort(dos_grados)
limite_bajo = dos_grados[int(Nboot * 0.16)]
limite_alto = dos_grados[int(Nboot * 0.84)]
print('Año promedio en el que la temperatura habra cambiado dos grados '
      'es {}'.format(int(np.mean(dos_grados))))
print("El intervalo de confianza al 68% es: [{}, {}]".format(limite_bajo,
                                                             limite_alto))
