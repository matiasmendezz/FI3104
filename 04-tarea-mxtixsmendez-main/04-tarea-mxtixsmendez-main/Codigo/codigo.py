import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.patches as patches
#  constantes
Q = 2.427  # carga total

#  setup
Lx = 10
Ly = 15
N_pasos_x = 101
N_pasos_y = 151
h = Lx / (N_pasos_x - 1)
area_letra = 2 * (8 * 1 + 1 * 0.5)  # en cm^2
V = np.zeros((N_pasos_x, N_pasos_y))
w = 1


def esta_en_la_letra(x, y):
    """
    Verifica que el par ordenado (x, y) este dentro de la letra m. Si esta
    adentro retorna True.
    """
    if x >= -2.5 and x <= -1.5 and y >= -3.5 and y <= 3.5:
        return True
    elif x >= 1.5 and x <= 2.5 and y >= -3.5 and y <= 3.5:
        return True
    elif x >= -1.5 and x <= -0.5 and y >= 1.5 and y <= 2.5:
        return True
    elif x >= -0.5 and x <= 0.5 and y >= 0.5 and y <= 1.5:
        return True
    elif x >= 0.5 and x <= 1.5 and y >= 1.5 and y <= 2.5:
        return True
    else:
        return False


def rho(i, j, h):
    """
    Calcula la carga de un punto (x_i, x_j), si esta dentro de la letra
    entonces la carga es Q / Area_letra. Si esta fuera la carga es cero.
    """
    x_i = -5 + i * h
    y_j = -7.5 + j * h
    if esta_en_la_letra(x_i, y_j) is True:
        output = Q / area_letra
    else:
        output = 0
    return output


def sobre_relajacion(V, i, j, h, w):
    '''
    Implementa el metodo de sobre relajacion sobre V_ij, para esto se necesita
    la matriz V, las coordenadas i y j, h y el coeficiente
    de sobrerelajacion w.
    '''
    relajacion = (V[i + 1, j] + V[i - 1, j] + V[i, j - 1] + V[i, j + 1] -
                  (h ** 2) * rho(i, j, h))
    sobre_relajacion = (1 - w) * V[i, j] + w / 4 * relajacion
    return sobre_relajacion


def una_iteracion(V, N_pasos_x, N_pasos_y, h, w):
    '''
    Función que hace una iteracion del algoritmo de sobre-relajacion para
    integrar la funcion V. La iteracion se divide en trs partes, la primera
    aplica sobre-relajacion en los puntos que no son adyacente y tampoco
    parte del segmento de CB derivativa, mientras
    que la
    segunda parte aplica el algoritmo de iteracion

    '''
    for i in range(1, N_pasos_x-1):
        for j in range(1, N_pasos_y-1):
            if -20 <= i <= 80 and j == 20:  # segmento de cb
                V[i, j] = V[i, j-1] - h   # derivada discreta
            elif j == 19:  # segmento abajo de la cb
                relajacion = (V[i+1, j] + V[i-1, j] + V[i, j-1] + V[i, j+1] -
                              h ** 2 * rho(i, j, h) + h)
                sobre_relaja = (1-w)*V[i, j]+w/4*relajacion
                V[i, j] = sobre_relaja
            elif j == 21:  # segmento de arriba de la cb
                relajacion = (V[i + 1, j] + V[i - 1, j]+V[i, j - 1] +
                              V[i, j + 1]-h ** 2 * rho(i, j, h) - h)
                sobre_relaja = (1 - w) * V[i, j] + w / 4 * relajacion
                V[i, j] = sobre_relaja
            else:  # caso general
                V[i, j] = sobre_relajacion(V, i, j, h, w)
    return V


def convergio(V, V_anterior, rtol):
    '''
    Funcion que calcula si V c/r a V anterior connverge en un radio de
    tolerancia rtol.
    '''
    not_zero = V_anterior != 0
    dif_relativa = (V[not_zero]-V_anterior[not_zero])/V_anterior[not_zero]
    output = np.fabs(dif_relativa).max() < rtol
    return output


def converger(N_pasos_x, N_pasos_y, h, w, rtol):
    '''
    Funcion que itera el metodo de sobrerelajacion sobre la funcion V
    hasta hacerla converger a un rtol determinado.
    Grafica V y ademas muestra el numero de iteraciones.
    '''
    V = np.zeros((N_pasos_x, N_pasos_y))
    V = una_iteracion(V, N_pasos_x, N_pasos_y, h, w)
    V_anterior = V.copy()
    V = una_iteracion(V, N_pasos_x, N_pasos_y, h, w)
    contador = 2
    while not convergio(V, V_anterior, rtol) and contador < 10000:
        contador += 1
        V_anterior = V.copy()
        V = una_iteracion(V, N_pasos_x, N_pasos_y, h, w)
    V = V.transpose()
    # Grafico de V en 2d
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pos = ax.imshow(V, origin='lower', extent=[-5, 5, -7.5, 7.5], cmap='hot')
    ax.axhline(0, linestyle=':')
    ax.axvline(0, linestyle=':')
    ax.set_xlabel('Eje X [cm]', fontsize=15)
    ax.set_ylabel('Eje Y [cm]', fontsize=15)
    ax.set_title('Potencial Electrostático [V]', fontsize=15)
    fig.colorbar(pos, ax=ax, cmap='hot')
    rect = patches.Rectangle((-2.5, -3.5), 5, 7, linewidth=1, edgecolor='g',
                             facecolor='none')
    ax.add_patch(rect)
    ax.contour(V, origin='lower', extent=[-5, 5, -7.5, 7.5], colors='k')
    fig.savefig('iteracion con cb con rtol={} y w={}.pdf'.format(rtol, w))
    print('Para w={} y rtol={} el n° de iteraciones es: {}'.format(w, rtol,
          contador))
    """
    # Grafico 3D
    fig_1 = plt.figure()
    ax_1 = fig_1.add_subplot(111, projection='3d')
    y = np.linspace(-5, 5, N_pasos_y)
    x = np.linspace(-7.5, 7.5, N_pasos_x)
    X, Y = np.meshgrid(x, y)
    ax_1.set_xlabel('Eje X [cm]', fontsize=15)
    ax_1.set_ylabel('Eje Y [cm]', fontsize=15)
    surf = ax_1.plot_surface(X, Y, V, cmap='hot')
    fig_1.colorbar(surf, ax=ax_1)
    fig_1.savefig('grafico 3d con cb con rtol={} y w={}.pdf'.format(rtol, w))
    fig_1.show()

    """


# Comprobamos que la letra M este bien hecha
# Asignamos las cargas en la matriz V
for i in range(N_pasos_x):
    for j in range(N_pasos_y):
        V[i, j] = rho(i, j, h)

# Grafico que verifica el dibujo de M
V = V.transpose()
fig_1 = plt.figure(1)
fig_1.clf()
ax_1 = fig_1.add_subplot(111)
ax_1.imshow(V, origin='lower', extent=[-5, 5, -7.5, 7.5], cmap='hot')
ax_1.axhline(0, linestyle=':')
ax_1.axvline(0, linestyle=':')
ax_1.set_xlabel('Eje X [cm]', fontsize=15)
ax_1.set_ylabel('Eje Y [cm]', fontsize=15)
ax_1.set_title('Letra M', fontsize=15)
rect = patches.Rectangle((-2.5, -3.5), 5, 7, linewidth=1, edgecolor='g',
                         facecolor='none')
ax_1.add_patch(rect)
fig_1.savefig('letra M.pdf')


# Iteracion sin la condicion de borde derivativa
# Comprueba metodo de relajamiento sin cb
V = np.zeros((N_pasos_x, N_pasos_y))  # Reseteamos valores
for i in range(100):
    for i in range(1, N_pasos_x - 1):
        for j in range(1, N_pasos_y - 1):
            V[i, j] = sobre_relajacion(V, i, j, h, w)

V = V.transpose()
fig_2 = plt.figure(2)
fig_2.clf()
ax_2 = fig_2.add_subplot(111)
pos = ax_2.imshow(V, origin='lower', extent=[-5, 5, -7.5, 7.5], cmap='hot')
fig_2.colorbar(pos, ax=ax_2, cmap='hot')
ax_2.axhline(0, linestyle=':')
ax_2.axvline(0, linestyle=':')
ax_2.set_xlabel('Eje X [cm]', fontsize=15)
ax_2.set_ylabel('Eje Y [cm]', fontsize=15)
ax_2.set_title('Potencial Electrostático [V]', fontsize=15)
rect = patches.Rectangle((-2.5, -3.5), 5, 7, linewidth=1, edgecolor='g',
                         facecolor='none')
ax_2.add_patch(rect)
fig_2.savefig('iteracion sin cb.pdf')


# Iteracion final
# Analisis convergencia
W = [1, 1.4, 1.6, 1.8]
for i in range(len(W)):
    w = W[i]
    converger(N_pasos_x, N_pasos_y, h, w, 0.1)
    converger(N_pasos_x, N_pasos_y, h, w, 0.01)
