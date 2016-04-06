#!/usr/bin/env python3

"""
Plot Poincare Section.
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif')


def poincare_section_sets(datafile, colx, coly, colvx):
    datf = open(datafile, 'r')
    data = np.array([0, 0, 0], dtype=np.longdouble)
    data1 = np.array([0, 0, 0], dtype=np.longdouble)
    firstline = datf.readline()
    firstline = firstline.strip()
    columns = firstline.split()
    cols = [colx, coly, colvx]
    secpoint = [[], []]
    for i in range(3):
        data[i] = np.longdouble(columns[cols[i]])
    for line in datf:
        for i in range(3):
            data1[i] = data[i]
            line = line.strip()
            columns = line.split()
        for i in range(3):
            data[i] = np.longdouble(columns[cols[i]])
        if (data[1] > 0 and data1[1] < 0):
            secpoint[0].append((data1[0] + data[0]) / 2)
            secpoint[1].append((data1[2] + data[2]) / 2)
        elif (data[1] == 0 and data1[1] < 0):
            secpoint[0].append(data[0])
            secpoint[1].append(data[2])
    datf.close()
    return secpoint

if __name__ == '__main__':
    datafile = '../data/poincare_section.dat'
    xyvxcols = [2, 3, 4]
    poincare_sec = poincare_section_sets(datafile, xyvxcols[0], xyvxcols[1],
                                         xyvxcols[2])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(poincare_sec[0], poincare_sec[1], marker='.', markersize=3.0,
            color='r', linestyle='None')
    ax.plot(0.499, 0, 'o')
    ax.annotate('$L_{4}$', xy=(0.499, 0), xytext=(0.49, -0.1),
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="angle,angleA=0,angleB=90,rad=10")
                )
    ax.set_xlabel('$x$', fontsize=14)
    ax.set_ylabel('$dx/dt$', fontsize=14)
    ax.set_xlim(0.42, 0.55)
    ax.set_ylim(-0.4, 0.4)
    # plt.axis('equal')
    ax.set_title('Poincaré Section', fontsize=16)
    fig.savefig('../img/poincare_section.png')
    plt.show()
