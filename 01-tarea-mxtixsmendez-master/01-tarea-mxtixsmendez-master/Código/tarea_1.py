import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.special import gamma

# rut: 20.243.742-7 --> ultimos 3 digitos antes del guion: 742
k = 4.742


def gamma_1(z, N=10 ** 5):
    """
    Calcula la función gamma de z, aplicando el método del trapecio.
    Este método se aplica en la integral de gamma luego
    de haber hecho el c.v y=e^(-x)
    """
    a = 10 ** (-10)  # límite inferior de la integral, este es una aproximación de 0
    b = 1  # límite superior de la integral
    dx = (b - a) / N
    x_i = np.linspace(0, 1, N + 1)  # vector x_i=a+i*dx
    x_i[0] = a  # aproximación de 0
    f_0 = math.log(1 / a) ** (z - 1)
    f_N = math.log(1 / b) ** (z - 1)
    f_i = np.log(1 / x_i) ** (z - 1)  # f_i=f(x_i)
    f_i[N] = 0
    f_i[0] = 0  # hacemos cero esto porque queremos sumar entre f_1 y f_{N-1}
    integral = (f_0 * 0.5 + f_N * 0.5 + sum(f_i)) * dx  # formula compuesta
    return integral


def precision_gamma(f, z, N=10 ** 5, epsilon=1.0e-6):
    """
    Mejora la precisión de la función gamma_1 a través
    del uso de la convergencia relativa. Retorna gamma(z) con
    una precisión de epsilon.
    """
    I_N = f(z, N)  # Integral para N segmentos usando trapecio
    I_2N = f(z, 2 * N)  # Integral para 2N segmentos usando trapecio
    tol = math.fabs((I_2N - I_N) / I_N)  # tolerancia relativa
    while True:
        if tol > epsilon:
            N = 2 * N  # iteramos denuevo para seguir mejorando la precisión
            I_N = I_2N
            I_2N = f(z, 2 * N)
            tol = math.fabs((I_2N - I_N) / I_N)
        else:  # momento en que ya se llega a la tolerancia
            break  # por lo que termina la iteración
    S_2N = (4 / 3) * I_2N - (1 / 3) * I_N  # Simpson, se retorna este valor ya que tiene mejor error de arrastre
    return S_2N


fracc = 1 / (2 ** (k / 2) * precision_gamma(gamma_1, k / 2))  # fraccion que aparece en la formula de chi^2
# fracc se calcula fuera de la función g(y) para no calcular fracc siempre que se llame a g(y)
def f(y, N=10 ** 5):
    """
    Calcula la integral de X^2, usando el método del trapecio.
    Los límites de la integral son (0,y).  
    """
    a = 0  # límite inferior de la integral
    b = y  # límite superior de la integral
    dx = (b - a) / N
    f_0 = a ** (k / 2 - 1) * math.e ** (-a / 2)
    f_N = a ** (k / 2 - 1) * math.e ** (-a / 2)
    x_i = np.linspace(0, y, N + 1)  # vector x_i=a+i*dx
    f_i = x_i ** (k / 2 - 1) * math.e ** (-x_i / 2)  # f_i=f(x_i)
    f_i[N] = 0
    f_i[0] = 0  # hacemos cero esto porque queremos sumar entre f_1 y f_{N-1}
    integral = (f_0 * 0.5 + f_N * 0.5 + sum(f_i)) * fracc * dx
    return integral


def precision_f(f, z, N=10 ** 5, epsilon=1.0e-5):
    """
    Mejora la precisión de la función f a través
    del uso de la convergencia relativa. Retorna f(z) con
    una precisión de epsilon.
    """
    I_N = f(z, N)  # Integral para N segmentos usando trapecio
    I_2N = f(z, 2 * N)  # Integral para 2N segmentos usando trapecio
    tol = math.fabs((I_2N - I_N) / I_N)  # tolerancia relativa
    while True:
        if tol > epsilon:
            N = 2 * N  # iteramos denuevo para seguir mejorando la precisión
            I_N = I_2N
            I_2N = f(z, 2 * N)
            tol = math.fabs((I_2N - I_N) / I_N)
        else:  # momento en que ya se llega a la tolerancia
            break  # por lo que termina la iteración
    S_2N = (4 / 3) * I_2N - (1 / 3) * I_N  # Simpson, se retorna este valor ya que tiene mejor error de arrastre
    return S_2N


def g(x):
    """
    g(x)=f(x)-0.95 función que contiene como raíz la solución
    de f(x)=0.95
    """
    return precision_f(f, x) - 0.95


# grafico de la funcion g(x), esto lo hacemos para encontrar un
# intervalo que encierre a la raíz
n = 100
x = np.linspace(0, 20, n)
y = f(x) - 0.95  # g(x) usamos esta definición para graficar más rápido
plt.grid(True, zorder=6)
plt.plot(x, y, 'b')
plt.axhline(0, color='0.5')
plt.title('Gráfico f(x) vs. x', fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('f(x)', fontsize=15)
plt.show()
# Despues de ver el grafico vemos que la raiz esta encerrada entre 7.5 y 12.5
# Para encontrar la raiz utilizamos el método de la bisección
a = 7.5
b = 12.5
tol = 1.0e-5  # precisión de la solución
p = (a + b) / 2  # punto medio
while math.fabs(b - a) > tol:
    """
    Algoritmo de la bisección
    """
    if g(p) * g(a) > 0:
        a = p
    elif g(p) * g(a) < 0:
        b = p
    p = (a + b) / 2
# 1e-5; p=10.668931007385254 
# 1e-4; p=10.668907165527344
# 1e-3; p=10.66864013671875
# 1e-2; p= 10.6689453125


print(p) #valor de a
