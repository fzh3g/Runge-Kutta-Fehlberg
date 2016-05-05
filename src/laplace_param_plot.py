#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import read_data

fig = plt.figure()
ax = fig.add_subplot(111)

data0 = read_data.read_cols('/home/faxiang/coding/matlab/laplace_param_test_0.dat')
data1 = read_data.read_cols('/home/faxiang/coding/matlab/laplace_param_test_1.dat')
data2 = read_data.read_cols('/home/faxiang/coding/matlab/laplace_param_test_2.dat')

ax.plot(data0[0], data0[1], 'r-', label='j=1,k=1/2')
ax.plot(data1[0], data1[1], 'r--', label='j=1,k=3/2')
ax.plot(data2[0], data2[1], 'r-.', label='j=1,k=5/2')

ax.set_xlabel(r'$\alpha{}$', fontsize=14)
ax.set_ylabel(r'$b_{k}(\alpha{})$', fontsize=14)
ax.set_ylim(0, 8)
ax.set_title("Laplace's Parameter", fontsize=16)

plt.legend(loc=2)
plt.savefig('laplace_param.png')
plt.show()
