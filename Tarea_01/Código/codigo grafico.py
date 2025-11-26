import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.special import gamma
from tarea 1.py import *

plt.figure(1)

n=100
x=np.linspace(0,20,n)
y=f(x)-0.95 #g(x), usamos esta definición para graficar más rápido
plt.grid(True, zorder=6)
plt.plot(x,y,'b')
plt.axhline(0, color='0.5')
plt.title('g(x) vs. x',fontsize=20)
plt.xlabel('x',fontsize=15)
plt.ylabel('g(x)',fontsize=15)

plt.figure(2)
x_1=np.linspace(0,10,10*20)
y_1=x_1**(k/2-1)*math.e**(-x_1)
x_2=np.linspace(1.0e-200,1,10000)
y_2=np.log(1/x_2)**(k/2-1)
plt.grid(True, zorder=10)
plt.plot(x_1,y_1,'black')
plt.fill(x_1, y_1)
plt.xlabel('x',fontsize=15)
plt.ylabel('$x^{k/2-1}e^{-x}$',fontsize=15)

plt.figure(3)
plt.grid(True, zorder=10)
plt.semilogy(x_2,y_2,'black')
plt.fill(x_2,y_2)
plt.xlabel('y',fontsize=15)
plt.ylabel('$ln(1/y)^{k/2-1}$',fontsize=15)
plt.show()


