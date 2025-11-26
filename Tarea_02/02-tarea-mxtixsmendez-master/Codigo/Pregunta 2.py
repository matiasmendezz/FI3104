import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45
from mpl_toolkits.mplot3d import Axes3D
# Constantes
sigma = 10
beta = 8/3
rho = 28
Y_0 = [1, 1, 1]
t_0 = 0


def fun_a_integrar(t, Y):
    """
    Función a integrar numéricamente
    Inputs:
    t : [float], tiempo
    y: [np.array], [x,y,z]
    """
    global sigma, beta, rho
    output = np.array([1, 1, 1])
    output[0] = sigma * (Y[1] - Y[0])
    output[1] = Y[0] * (rho - Y[2]) - Y[1]
    output[2] = Y[0]*Y[1]-beta*Y[2]
    return output


# Integracion numerica
t_max = 50
res = RK45(fun_a_integrar, t_0, Y_0, t_max, vectorized=True, rtol=1e-13)
N = 90000
t = np.zeros(N)
Y = np.zeros((N, 3))
Y[0] = Y_0
t[0] = t_0
for i in range(N):
    res.step()
    t[i] = res.t
    Y[i] = res.y


# Integracion numerica con una pequeña perturbacion
Y_0_1 = [1, 1+0.0427, 1]
res_1 = RK45(fun_a_integrar, t_0, Y_0_1, t_max, vectorized=True, rtol=1e-13)
t_1 = np.zeros(N)
Y_1 = np.zeros((N, 3))
Y_1[0] = Y_0_1
t_1[0] = t_0
for i in range(N):
    res_1.step()
    t_1[i] = res_1.t
    Y_1[i] = res_1.y

# Graficos
# Grafico 3D
# Solucion sin perturbacion
fig1 = plt.figure(1)
ax1_1 = fig1.add_subplot(111, projection='3d')
ax1_1.plot(Y[:, 0], Y[:, 1], Y[:, 2])
ax1_1.set_xlabel('x')
ax1_1.set_ylabel('y')
ax1_1.set_zlabel('z')
plt.savefig('grafico 3D sin perturbacion.pdf')
# Solucion con perturbacion
fig2 = plt.figure(2)
ax1_2 = fig2.add_subplot(111, projection='3d')
ax1_2.plot(Y_1[:, 0], Y_1[:, 1], Y_1[:, 2], 'orange')
ax1_2.set_xlabel('x')
ax1_2.set_ylabel('y')
ax1_2.set_zlabel('z')
plt.savefig('grafico 3D con perturbacion.pdf')
# Grafico coordenada por tiempo
fig3 = plt.figure(3)
ax2 = fig3.add_subplot(311)  # x vs t
ax3 = fig3.add_subplot(312)  # y vs t
ax4 = fig3.add_subplot(313)  # z vs t
ax2.plot(t, Y[:, 0], label='sin perturbación')
ax2.plot(t_1, Y_1[:, 0], label='con perturbación')
ax3.plot(t, Y[:, 1], label='sin perturbación')
ax3.plot(t_1, Y_1[:, 1], label='con perturbación')
ax4.plot(t, Y[:, 2], label='sin perturbación')
ax4.plot(t_1, Y_1[:, 2], label='con perturbación')
ax2.set_xlabel('t', fontsize=15)
ax2.set_ylabel('x', fontsize=15)
ax3.set_xlabel('t', fontsize=15)
ax3.set_ylabel('y', fontsize=15)
ax4.set_xlabel('t', fontsize=15)
ax4.set_ylabel('z', fontsize=15)
ax2.legend(loc=2, fontsize=10,  bbox_to_anchor=(0.1, 0.4, 0.5, 0.5))
ax3.legend(loc=2, fontsize=10,  bbox_to_anchor=(0.08, 0.4, 0.5, 0.5))
ax4.legend(fontsize=10,  bbox_to_anchor=(0.05, 0., 0.5, 0.5))
plt.savefig('grafico coordenadas vs tiempo.pdf')

# Graficos por plano
# Plano XY
fig4 = plt.figure(4)
ax5 = fig4.add_subplot(211)
ax5.plot(Y[:, 0], Y[:, 1], label='sin perturbación')
ax5.set_ylabel('x', fontsize=15)
ax5.set_xlabel('y', fontsize=15)
ax5.legend()

ax6 = fig4.add_subplot(212)
ax6.plot(Y_1[:, 0], Y_1[:, 1], 'green', label='con perturbación')
ax6.set_ylabel('x', fontsize=15)
ax6.set_xlabel('y', fontsize=15)
ax6.legend()
plt.savefig('Grafico plano XY.pdf')

# Plano YZ
fig5 = plt.figure(5)
ax7 = fig5.add_subplot(211)
ax7.plot(Y[:, 1], Y[:, 2], label='sin perturbación')
ax7.set_ylabel('y', fontsize=15)
ax7.set_xlabel('z', fontsize=15)
ax7.legend()

ax8 = fig5.add_subplot(212)
ax8.plot(Y_1[:, 1], Y_1[:, 2], 'green', label='con perturbación')
ax8.set_ylabel('y', fontsize=15)
ax8.set_xlabel('z', fontsize=15)
ax8.legend()
plt.savefig('Grafico plano YZ.pdf')

# Plano ZX
fig6 = plt.figure(6)
ax8 = fig6.add_subplot(211)
ax8.plot(Y[:, 2], Y[:, 0], label='sin perturbación')
ax8.set_ylabel('z', fontsize=15)
ax8.set_xlabel('x', fontsize=15)
ax8.legend()

ax9 = fig6.add_subplot(212)
ax9.plot(Y_1[:, 2], Y_1[:, 0], 'green', label='con perturbación')
ax9.set_ylabel('z', fontsize=15)
ax9.set_xlabel('x', fontsize=15)
ax9.legend()
plt.savefig('Grafico plano ZX.pdf')
plt.show()
