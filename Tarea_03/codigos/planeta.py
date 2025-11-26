from scipy.integrate import RK45
import numpy as np


class Planeta(object):
    """
    Esta clase corresponde a la descripción de un planeta. Recibe una condición inicial
    que describe la posición y velocidad del planeta en ese instante. Y recibe un alpha
    como opcional, que corresponde a la correción del potencial gravitatorio sobre
    el planeta.
    """
    def __init__(self, condicion_inicial, alpha=0):
        """
        __init__ es un método especial que se usa para inicializar las
        instancias de una clase.

        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        """
        self.y_actual = condicion_inicial
        self.y_previo = [0, 0, 0, 0]
        self.t_actual = 0.
        self.alpha = alpha

    def ecuacion_de_movimiento(self, t, y):
        """
        Implementa la ecuación de movimiento, como sistema de ecuaciones de
        primer orden.
        """
        x, y, vx, vy = y
        fx = (-x)/(x**2+y**2)**(3/2) + self.alpha*(2*x)/(x**2+y**2)**(2)
        fy = (-y)/(x**2+y**2)**(3/2) + self.alpha*(2*y)/(x**2+y**2)**(2)
        return np.array([vx, vy, fx, fy])

    def avanza_rk4(self, dt):
        """
        Toma la condición actual del planeta y avanza su posición y velocidad
        en un intervalo de tiempo dt usando el método de rk4. El método no
        retorna nada, pero modifica los valores de self.y_actual, self.y_prev y
        de self.t_actual.
        """
        t_n = self.t_actual
        y_n = self.y_actual
        # Constantes k
        k1 = self.ecuacion_de_movimiento(t_n, y_n)
        k2 = self.ecuacion_de_movimiento(t_n + (1/2)*dt, y_n + (1/2)*k1*dt)
        k3 = self.ecuacion_de_movimiento(t_n + (1/2)*dt, y_n + (1/2)*k2*dt)
        k4 = self.ecuacion_de_movimiento(t_n + dt, y_n + k3*dt)
        self.y_previo = self.y_actual
        self.y_actual = y_n + (1/6)*dt * (k1 + 2*k2 + 2*k3 + k4)
        self.t_actual += dt

    def avanza_verlet(self, dt):
        """
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de verlet. El método no
        retorna nada, pero modifica los valores de self.y_actual, self.y_prev
        y de self.t_actual.
        """
        r_n = np.array(self.y_actual[:2])
        r_n_menos_1 = np.array(self.y_previo[:2])
        fuerza = self.ecuacion_de_movimiento(self.t_actual, self.y_actual)[2:]
        fuerza = np.array(fuerza)
        r_n_mas_1 = 2 * r_n - r_n_menos_1 + dt ** 2 * fuerza  # verlet posicion
        v_n_mas_1 = (r_n_mas_1-r_n)/dt  # verlet velocidad
        y_sgte = [r_n_mas_1[0], r_n_mas_1[1], v_n_mas_1[0], v_n_mas_1[1]]
        self.y_previo = self.y_actual
        self.t_actual += dt
        self.y_actual = y_sgte

    def avanza_beeman(self, dt):
        """
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de beeman. El método no
        retorna nada, pero modifica los valores de self.y_actual, self.y_prev
        y de self.t_actual.
        """
        t0 = self.t_actual
        r_n = np.array(self.y_actual[:2])
        v_n = np.array(self.y_actual[2:])
        fuerza_n = self.ecuacion_de_movimiento(t0, self.y_actual)[2:]
        fuerza_n = np.array(fuerza_n)
        fuerza_n_menos_1 = self.ecuacion_de_movimiento(t0, self.y_previo)[2:]
        r_n_mas_1 = r_n+dt*v_n+(2/3*fuerza_n-1/6*fuerza_n_menos_1)*dt**2
        v_n_mas_1 = v_n+1/2*(3*fuerza_n-fuerza_n_menos_1)*dt
        y_sgte = [r_n_mas_1[0], r_n_mas_1[1], v_n_mas_1[0], v_n_mas_1[1]]
        self.y_previo = self.y_actual
        self.t_actual += dt
        self.y_actual = y_sgte

    def energia_total(self):
        """
        Calcula la enérgía total del sistema en las condiciones actuales.
        """
        x, y, vx, vy = self.y_actual
        radio = np.sqrt(x**2+y**2)
        energia_cinetica = 1/2*(vx**2+vy**2)
        energia_potencial = -1/radio+self.alpha/radio**2
        energia = energia_cinetica+energia_potencial
        return energia
