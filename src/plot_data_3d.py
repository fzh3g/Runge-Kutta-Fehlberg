#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import read_data


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
datafile = "../data/lorenz.dat"
datacols = [2, 3, 4]
data = read_data.read_cols(datafile, cols=datacols, header=1)

ax.plot(data[0], data[1], data[2], marker='.', markersize=2.0,
        linestyle='None', color='r')
ax.set_title('Lorenz Attractor', fontsize=16)
ax.dist = 12

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# bbox = (-2.0, 2.0)
# xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
# ax.set_zlim3d(zmin, zmax)
# ax.set_xlim3d(xmin, xmax)
# ax.set_ylim3d(ymin, ymax)

plt.savefig('../img/lorenz.png')
plt.show()
