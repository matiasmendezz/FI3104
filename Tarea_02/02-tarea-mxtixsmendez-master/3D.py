#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')

t = np.arange(0, 15 * np.pi, 0.01)

a = 10.
b = 8.
c = 6

wx = 3.2454
wy = 5.4234
wz = 9.1423

u = a * np.sin(wx * t)
v = b * np.sin(wy * t)
w = c * np.sin(wz * t)

ax.plot(u, v, w)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
