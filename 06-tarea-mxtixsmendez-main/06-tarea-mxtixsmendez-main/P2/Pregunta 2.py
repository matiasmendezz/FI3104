import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from lmfit import (minimize, Parameters, report_fit)
mpl.use('TkAgg')  # maña de Pycharm para abrir graficos interactivos
"""
Forma paramatrica del espectro escogida:
y = a_1 * x+ a_2 + a_3 * exp( -(x-a_4)^2/ (a_5^2) )
Contiene 5 parametros: a_1, a_2, a_3, a_4, a_5
"""


def linea(x, a_1, a_2):
    output = a_1 + a_2 * x
    return output


def gaussiana(x, a_3, a_4, a_5):
    # a_3 * exp( -(x-a_4)^2/ (a_5^2)
    output = a_3*np.exp(-(x-a_4)**2/(a_5**2))
    return output


def f_modelo(x, a_1, a_2, a_3, a_4, a_5):
    output_1 = linea(x, a_1, a_2)
    output_2 = gaussiana(x, a_3, a_4, a_5)
    output = output_1 + output_2
    return output


# Lectura de Datos
data = np.loadtxt('espectro.dat')
n = 122
x = np.zeros(n)  # longitud de onda
y = np.zeros(n)  # flujo
for i in range(len(data)):
    x[i] = data[i][0]
    y[i] = data[i][1]
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
ax.plot(x, y, label='Datos')
ax.legend()
ax.set_xlabel(r'Wavelength [$\AA$]', fontsize=15)
ax.set_ylabel(r'$F_\nu$ [erg $s^{-1}$ $Hz^{-1}$ $cm^{-2}$]', fontsize=15)
fig.savefig('Grafico datos.pdf')

# Estimacion errores gaussianos
i = 0
while x[i] <= 6550:
    i += 1
x_recta = x[0:i+1]  # Espectro desde 6460 hasta 6550 Hz
y_recta = y[0:i+1]
print('Estimacion de la linea desde {} [\AA] hasta {} [\AA]'.format(x[0], x[i]))
params = [1, 0.5]  # ojimetro
param_opt, param_cov = curve_fit(linea, x_recta, y_recta, p0=params)
a_1, a_2 = param_opt
# Grafico linea y datos
fig_2 = plt.figure(2)
fig_2.clf()
ax_2 = fig_2.add_subplot(111)
ax_2.plot(x, y, label='Datos')
ax_2.plot(x_recta, linea(x_recta, a_1, a_2), label='curve_fit')
ax_2.set_xlabel(r'Wavelength [$\AA$]', fontsize=15)
ax_2.set_ylabel(r'$F_\nu$ [erg $s^{-1}$ $Hz^{-1}$ $cm^{-2}$]', fontsize=15)
ax_2.legend()
fig_2.savefig('Grafico recta en espectro.pdf')

# Estimacion del sigma
error_recta = y_recta-linea(x_recta, a_1, a_2)
sigma = np.std(error_recta)
print('sigma obtenido: ', sigma)
# Grafico del error
x_hist = np.arange(-1.5, 1.5, 0.1)
fig_3 = plt.figure(3)
fig_3.clf()
ax_3 = fig_3.add_subplot(111)
ax_3.hist(error_recta, bins=10, density=True, edgecolor='black', label='Datos')
ax_3.set_xlabel(r'Error de medición [$\AA$] ', fontsize=15)
ax_3.set_ylabel('Frecuencia', fontsize=15)
ax_3.legend()
fig_3.savefig('Histograma.pdf')

# Fiteo final
# Asumimos un error constante para cada dato
error = np.zeros(len(x))
for i in range(len(error)):
    error[i] = sigma


def residual(params, x, y, yerr):
    a_1 = params['a_1']
    a_2 = params['a_2']
    a_3 = params['a_3']
    a_4 = params['a_4']
    a_5 = params['a_5']
    modelo = a_1 + a_2 * x + a_3*np.exp(-(x-a_4)**2/(a_5**2))
    return (y-modelo)/yerr


# Parametros iniciales
params = Parameters()
params.add('a_1', value=1.3)
params.add('a_2', value=0.01)
params.add('a_3', value=1)
params.add('a_4', value=6565)
params.add('a_5', value=10)

fiteo = minimize(residual, params, args=(x, y, error))
report_fit(fiteo)
a_1 = fiteo.params['a_1'].value
a_2 = fiteo.params['a_2'].value
a_3 = fiteo.params['a_3'].value
a_4 = fiteo.params['a_4'].value
a_5 = fiteo.params['a_5'].value

# Grafico fiteo final
fig_4 = plt.figure(4)
fig.clf()
ax_4 = fig_4.add_subplot(111)
ax_4.plot(x, y, label='Datos')
ax_4.plot(x, f_modelo(x, a_1, a_2, a_3, a_4, a_5), label='lmfit')
ax_4.legend()
ax_4.set_xlabel(r'Wavelength [$\AA$]', fontsize=15)
ax_4.set_ylabel(r'$F_\nu$ [erg $s^{-1}$ $Hz^{-1}$ $cm^{-2}$]', fontsize=15)
fig_4.savefig('Grafico fiteo espectro.pdf')

""""
Resultados entregados:


Estimacion de la linea desde 6460.0 [\AA] hasta 6550.7 [\AA]
sigma obtenido:  6.012839265154404e-19
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 121
    # data points      = 122
    # variables        = 5
    chi-square         = 143.941784
    reduced chi-square = 1.23027166
    Akaike info crit   = 30.1773236
    Bayesian info crit = 44.1974288
[[Variables]]
    a_1:  8.8769e-17 +/- 6.8187e-18 (7.68%) (init = 1.3)
    a_2:  7.8026e-21 +/- 1.0395e-21 (13.32%) (init = 0.01)
    a_3: -1.0068e-17 +/- 4.3879e-19 (4.36%) (init = 1)
    a_4:  6563.22332 +/- 0.16321141 (0.00%) (init = 6565)
    a_5:  4.60781481 +/- 0.23426078 (5.08%) (init = 10)
[[Correlations]] (unreported correlations are < 0.100)
    C(a_1, a_2) = -1.000
    C(a_3, a_5) =  0.547
"""
