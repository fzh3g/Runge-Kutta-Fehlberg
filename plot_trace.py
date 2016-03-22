#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl


def main():
    f = open('twobody_output.dat', 'r')
    f.readline()
    x = []
    y = []
    for line in f:
        line = line.strip()
        columns = line.split()
        x.append(np.longdouble(columns[2]))
        y.append(np.longdouble(columns[3]))
    f.close()
    pl.figure(figsize=(8, 6))
    pl.plot(x, y, color='b')
    pl.xlabel("X")
    pl.ylabel("Y")
    pl.ticklabel_format(style='sci', scilimits=(0, 0))
    pl.axis('equal')
    pl.title("Orbit Trace")
    pl.savefig("orbit_trace.eps")
    pl.show()


if __name__ == '__main__':
    main()
