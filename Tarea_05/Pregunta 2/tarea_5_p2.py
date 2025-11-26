import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


"""
El codigo toma casi 4 horas en compilarse.
W(x) = 1 / (5*sqrt(pi)) * ( 2e^-(X-1.5)^2 + 3e^-(x+1)^2 )
Tira un error que presumo que es del ultimo elemento del arreglo
histograma_promedio.
La barra de error no se nota pero si esta incluida. Para n=10 y
N_muestra = 100_000 se grafica rapido y es mas visible la barra de error.
"""


def W(x):
    const = 1/(5*np.sqrt(np.pi))
    output = const*(2*np.e**(-(x-1.5)**2) + 3*np.e**(-(x+1)**2))
    return output


normalizacion = quad(W, -100, 100)[0]
# print(normalizacion)
x = np.arange(0, 10, 0.001)
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(x, W(x)/normalizacion)


# Metropolis
# x_p = x_0 + d * r

def paso_metropolis(x_0, d=0.1):
    aprueba = 0
    x_p = x_0 + d * np.random.uniform(low=-1, high=1)
    if W(x_p) / W(x_0) > np.random.uniform(low=0, high=1):
        aprueba = 1
        x_0 = x_p
    return x_0, aprueba


def genera_muestra_metropolis(x_0, N_muestra, d=0.1):
    muestra_metropolis = np.zeros(N_muestra)
    muestra_metropolis[0] = x_0
    apruebos = 1
    for i in range(1, N_muestra):
        x_p, aprueba = paso_metropolis(muestra_metropolis[i-1], d=d)
        apruebos += aprueba
        muestra_metropolis[i] = x_p
    # Tasa de aprobacion
    tasa = apruebos / N_muestra
    return muestra_metropolis, tasa


# Datos
np.random.seed(271)
N_muestra = 10_000_000
d = 3.75
x_0 = 1.8
n = 100  # Repeticiones de exp
histogramas = list()
aprobaciones = np.zeros(n)
bin_edges = np.arange(-8, 8, 0.1)
bin_centers = np.zeros(len(bin_edges))
# Creamos arreglo con los centros de cada caja

for i in range(len(bin_edges)-1):
    prom = (bin_edges[i] + bin_edges[i+1])*0.5
    bin_centers[i] = prom

# Generamos n histogramas que se guardan en la lista histogramas
for i in range(n):
    muestra, aprobacion = genera_muestra_metropolis(x_0, N_muestra, d=d)
    histograma = plt.hist(muestra, bins=bin_edges, density=True)[0]
    print('Se agrego histograma numero {}'.format(i))
    histogramas.append(histograma)
    print(aprobacion)
    aprobaciones[i] = aprobacion

# Creamos una lista de listas
# Por cada caja se guardara en una lista los valores de los histogramas
# en esa caja
lista_cajas = list()
for i in range(len(bin_centers)):
    lista_cajas.append(list())
print('se pasan a guardar en las listas')
for i in range(len(histogramas)):
    histograma = histogramas[i]
    for j in range(len(histograma)):
        lista_cajas[j].append(histograma[j])

histograma_promedio = np.zeros(len(bin_centers))
histograma_std = np.zeros(len(bin_centers))

for i in range(len(histograma_promedio)):
    promedio_caja_i = np.mean(lista_cajas[i])
    std_caja_i = np.std(lista_cajas[i])
    histograma_promedio[i] = promedio_caja_i
    histograma_std[i] = std_caja_i
# Se comienza a graficar
plt.clf()
plt.plot(x, W(x)/normalizacion, 'r', label='Funcion original')
plt.plot(bin_centers, histograma_promedio, drawstyle='steps-mid',
         label='Metropolis')
plt.errorbar(bin_centers, histograma_promedio, yerr=histograma_std, ls='None',
             marker='None')
plt.xlabel('x', fontsize=20)
plt.ylabel('W(x)', fontsize=20)
plt.legend()
plt.title('Histograma promedio con 100 repeticiones', fontsize=20)
plt.savefig('Grafico histograma curva p2 corregido.pdf')
plt.show()
