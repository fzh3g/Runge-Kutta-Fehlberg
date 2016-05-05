#!/usr/bin/env python3

import numpy as np


def read_cols(filename, cols=[0, 1], header=0):
    """Read two columns of data in a file"""
    f = open(filename, 'r')
    if header > 0:
        for i in range(header):
            f.readline()
    numberofcol = len(cols)
    data = [[] for i in range(numberofcol)]
    for line in f:
        line = line.strip()
        columns = line.split()
        for i in range(numberofcol):
            data[i].append(np.longdouble(columns[cols[i]]))
    f.close()
    return data

