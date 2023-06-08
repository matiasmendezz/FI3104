import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
"""
Calculo de integrales para encontrar el centro de masa de la interseccion entre
el toro de ecuacion: z^2(sqrt(x^2+y^2)-3)^1 <= 1 y el cilindro de ecuacion
(x-2)^2+z^2<=1. Con densidad de masa: rho(x,y,z)=0.5*(x^2+y^2+z^2). Se calculan
cuatro integrales por el metodo de Monte Carlos para varias  dimensiones.
"""

np.random.seed(432)


def densidad(x, y, z):
    '''
    Calcula la densidad de sólido en el punto (x, y, z)
    Este punto debe estar adentro de la intersección del toro.
    '''
    rho = 0.5*(x**2+y**2+z**2)
    return rho


def esta_dentro(x, y, z):
    """
    Verifica si el punto (x, y, z) esta en la interseccion del toro
    con el cilindro. En caso de estar adentro retorna el valor booleano True.
    """
    toro = z**2+(np.sqrt(x**2+y**2)-3)**2
    cilindro = (x-2)**2 + z**2
    if toro <= 1 and cilindro <= 1:
        return True
    else:
        return False


def calcula_centro_de_masa(N):
    """
    Calcula el centro de masa de la intersección del toro con el cilindro
    con una muestra de N puntos. Para eso calcula 4 integrales a traves del
    metodo de Monte Carlo. Retorna las coordenadas del centro de
    masa ademas de la masa del objeto.
    """
    x = np.random.uniform(low=-4, high=4, size=N)
    y = np.random.uniform(low=-4, high=4, size=N)
    z = np.random.uniform(low=-1, high=1, size=N)
    volumen = 8.0 * 8.0 * 2.0
    '''
    Funciones a integrar:
    f_m = rho(x,y,z)
    f_x = x*rho(x,y,z)
    f_y = y*rho(x,y,z)
    f_z = z*rho(x,y,z)
    '''
    sum_f_m = sum_f_x = sum_f_y = sum_f_z = 0
    for i in range(N):
        x_i = x[i]
        y_i = y[i]
        z_i = z[i]
        if esta_dentro(x_i, y_i, z_i) is True:
            sum_f_m += densidad(x_i, y_i, z_i)
            sum_f_x += x_i * densidad(x_i, y_i, z_i)
            sum_f_y += y_i * densidad(x_i, y_i, z_i)
            sum_f_z += z_i * densidad(x_i, y_i, z_i)

    m = volumen * sum_f_m / N
    int_f_x = volumen * sum_f_x / N
    int_f_y = volumen * sum_f_y / N
    int_f_z = volumen * sum_f_z / N
    x_cm = int_f_x / m
    y_cm = int_f_y / m
    z_cm = int_f_z / m
    return m, x_cm, y_cm, z_cm


def repite_calculo_centro_de_masa(N, n):
    '''
    Repite el experimento de calcular el centro de masas  n veces. Con una
    muetra aleatoria de N puntos. Retorna el promedio de la posición del centro
    de masa y su desviaión estandar.
    '''
    m = np.zeros(n)
    x_cm = np.zeros(n)
    y_cm = np.zeros(n)
    z_cm = np.zeros(n)
    for i in range(n):
        m[i], x_cm[i], y_cm[i], z_cm[i] = calcula_centro_de_masa(N)
    prom_m = np.mean(m)
    prom_x_cm = np.mean(x_cm)
    prom_y_cm = np.mean(y_cm)
    prom_z_cm = np.mean(z_cm)
    err_m = np.std(m)
    err_x_cm = np.std(x_cm)
    err_y_cm = np.std(y_cm)
    err_z_cm = np.std(z_cm)
    prom = np.array([prom_m, prom_x_cm, prom_y_cm, prom_z_cm])
    err = np.array([err_m, err_x_cm, err_y_cm, err_z_cm])
    return prom, err


n = 1000  # repiticiones

arreglo_N = [100, 500, 1_000, 2000, 3000, 4000, 5_000, 10_000, 20_000, 30_000,
             40_000,  50_000, 100_000, 500_000]
arreglo_std_x = np.zeros(len(arreglo_N))
arreglo_std_y = np.zeros(len(arreglo_N))
arreglo_std_z = np.zeros(len(arreglo_N))

"""
"""
for i in range(len(arreglo_N)):
    print(arreglo_N[i])
    N = int(arreglo_N[i])
    prom, err = repite_calculo_centro_de_masa(N, n)
    print('Para N {:_} puntos:'.format(N))
    print('          masa m = {}'.format(prom[0]))
    print('          X_cm = {}'.format(prom[1]))
    print('          Y_cm = {}'.format(prom[2]))
    print('          Z_cm = {}'.format(prom[3]))
    arreglo_std_x[i] = err[1]
    arreglo_std_y[i] = err[2]
    arreglo_std_z[i] = err[3]


# Grafico desviacion estandar en x
fig_x = plt.figure(1)
fig_x.clf()
ax_x = fig_x.add_subplot(111)
ax_x.plot(arreglo_N, arreglo_std_x)
ax_x.set_ylabel(r'$\sigma_{Xcm}$', fontsize=20)
ax_x.set_xlabel('N', fontsize=20)
ax_x.set_title('Desviación estándar de $X_{cm}$ en función de N',
               fontsize=20)
fig_x.savefig('desv est en x.pdf')
# Grafico desviacion estandar en y
fig_y = plt.figure(2)
fig_y.clf()
ax_y = fig_y.add_subplot(111)
ax_y.plot(arreglo_N, arreglo_std_y)
ax_y.set_ylabel(r'$\sigma_{Ycm}$', fontsize=20)
ax_y.set_xlabel('N', fontsize=20)
ax_y.set_title('Desviación estándar de $Y_{cm}$ en función de N',
               fontsize=20)
fig_y.savefig('desv est en y.pdf')

# Grafico desviacion estandar en z
fig_z = plt.figure(3)
fig_z.clf()
ax_z = fig_z.add_subplot(111)
ax_z.plot(arreglo_N, arreglo_std_z)
ax_z.set_ylabel(r'$\sigma_{Zcm}$', fontsize=20)
ax_z.set_xlabel('N', fontsize=20)
ax_z.set_title('Desviación estándar de $Z_{cm}$ en función de N',
               fontsize=20)
fig_z.savefig('desv est en z.pdf')

