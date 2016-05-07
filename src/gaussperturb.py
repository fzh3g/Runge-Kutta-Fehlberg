#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import read_data

fig1 = plt.figure(figsize=(12, 8))
ax1 = fig1.add_subplot(221)
ax2 = fig1.add_subplot(222)
ax3 = fig1.add_subplot(223)
ax4 = fig1.add_subplot(224)

datafile = "../data/gauss_perturb.dat"
datacols = [0, 2, 3, 4, 5, 6, 7]
data = read_data.read_cols(datafile, cols=datacols, header=1)

gauss_time = data[0]
gauss_a = data[1]
gauss_e = data[2]
gauss_omega = data[3]
gauss_m = data[4]
gauss_x = data[5]
gauss_y = data[6]

ax1.plot(gauss_time, gauss_a, 'bo', markersize=3)
ax1.plot(gauss_time, gauss_a, 'g')
ax1.set_xlabel('t')
ax1.set_ylabel('a')
ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

ax2.plot(gauss_time, gauss_e, 'bo', markersize=3)
ax2.plot(gauss_time, gauss_e, 'g')
ax2.set_xlabel('t')
ax2.set_ylabel('e')
ax2.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))

ax3.plot(gauss_time, gauss_omega, 'bo', markersize=3)
ax3.plot(gauss_time, gauss_omega, 'g')
ax3.set_xlabel('t')
ax3.set_ylabel('$\omega$')
ax3.ticklabel_format(style='sci', scilimits=(0, 0))

ax4.plot(gauss_time, gauss_m, 'bo', markersize=3)
ax4.plot(gauss_time, gauss_m, 'g')
ax4.set_xlabel('t')
ax4.set_ylabel('M')
ax4.ticklabel_format(style='sci', scilimits=(0, 0))

fig1.savefig('gauss_perturb.png')

fig2 = plt.figure()
ax4 = fig2.add_subplot(111)
ax4.plot(gauss_x, gauss_y, 'g.', markersize=3)
# ax4.plot(gauss_x, gauss_y, 'g')
ax4.set_xlabel('x')
ax4.set_ylabel('y')
plt.axis('equal')

plt.show()
