import math
import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
from scipy.integrate import trapz
k = 4.742
"""
Solución a la tarea 1 usando los métodos de las librerias
Se usa para comparar resultados
"""
def integracion_funcion_chi2(z, N=10 ** 4, epsilon=1.0e-12):
    """
    Entrega la solución numérica de la integral con una precisión
    de epsilin
    """
    a=0
    b=z
    x=np.linspace(a,b,N)
    y=chi2.pdf(x,k)
    I_N=trapz(y,x)
    
    N=2*N
    x=np.linspace(a,b,N)
    y=chi2.pdf(x,k)
    I_2N =trapz(y,x)
    tol = math.fabs((I_2N - I_N) / I_N)  # tolerancia relativa
    while True:
        if tol > epsilon:
            N = 2 * N  # iteramos denuevo para seguir mejorando la precisión
            I_N = I_2N
            x=np.linspace(a,b,N)
            y=chi2.pdf(x,k)
            I_2N =trapz(y,x)
            tol = math.fabs((I_2N - I_N) / I_N)
        else:  # momento en que ya se llega a la tolerancia
            break  # por lo que termina la iteración
    S_2N = (4 / 3) * I_2N - (1 / 3) * I_N  # Simpson, se retorna este valor ya que tiene mejor error de arrastre
    return S_2N

from scipy import optimize
def g(x):
    return integracion_funcion_chi2(x)-0.95
#optimize.bisect(g,7.5,12.5)

# print('root=',optimize.bisect(g,7.5,12.5))

plt.figure()
n=100
x=np.linspace(0,20,n)
y=np.zeros(n)
for i in range(n):
    y[i]=g(x[i])
plt.grid(True, zorder=6)
plt.plot(x,y,'b')
plt.axhline(0, color='0.5')
plt.title('g(x) vs. x',fontsize=20)
plt.xlabel('x',fontsize=15)
plt.ylabel('g(x)',fontsize=15)
plt.axvline(10.668931007385254,label='$10^{-5}$',color='green')
plt.axvline(10.668907165527344,label='$10^{-4}$',color='blue')
plt.axvline(10.66864013671875,label='$10^{-3}$',color='red')
plt.axvline(10.6689453125,label='$10^{-2}$',color='purple')
plt.legend()
plt.show()