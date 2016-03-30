#!/usr/bin/env python3

"""
Usage: simple_plot.py [OPTIONS] DATAFILE

  Read two columns of data from DATAFILE and plot.

Options:
  --columns INTEGER...     Two columns of the data file to plot.
  --del-header INTEGER     Number of header lines to delete.
  --labels <TEXT TEXT>...  Labels of X and Y axes.
  --xlim <FLOAT FLOAT>...  Limit of X axis.
  --ylim <FLOAT FLOAT>...  Limit of Y axis.
  --title TEXT             Title of the figure.
  --figname TEXT           Name of the figure to save (default: simple_plot).
  --figtype TEXT           Format of the figure to save (default: eps).
  --sci                    Apply scientific notation to the axes.
  --equal                  Same aspect for X and Y axes.
  --show                   Show the figure.

Examples:
  $ ./simple_plot.py twobody_output.dat --show --columns 2 3 --equal \
    --del-header 1 --title 'Orbit Trace' --figname 'orbit_trace'
"""

from matplotlib import rc
import numpy as np
import matplotlib.pyplot as pl
import click


def read_two_col(filename, col0, col1, header):
    """Read two columns of data in a file"""
    f = open(filename, 'r')
    if header > 0:
        for i in range(header):
            f.readline()
    x = []
    y = []
    for line in f:
        line = line.strip()
        columns = line.split()
        x.append(np.longdouble(columns[col0]))
        y.append(np.longdouble(columns[col1]))
    f.close()
    return x, y


@click.group()
def simple_plot():
    pass


@simple_plot.command()
@click.option('--columns', nargs=2, default=(0, 1),
              help='Two columns of the data file to plot.')
@click.option('--del-header', default=0,
              help='Number of header lines to delete.')
@click.option('--labels', nargs=2, default=('$x$', '$y$'),
              help='Labels of X and Y axes.')
@click.option('--xlim', nargs=2, type=(float, float), default=(0, 0),
              help='Limit of X axis.')
@click.option('--ylim', nargs=2, type=(float, float), default=(0, 0),
              help='Limit of Y axis.')
@click.option('--title', default='', help='Title of the figure.')
@click.option('--figname', default='simple_plot',
              help='Name of the figure to save (default: simple_plot).')
@click.option('--figtype', default='eps',
              help='Format of the figure to save (default: eps).')
@click.option('--sci', is_flag=True,
              help='Apply scientific notation to the axes.')
@click.option('--equal', is_flag=True,
              help='Same aspect for X and Y axes.')
@click.option('--show', is_flag=True,
              help='Show the figure.')
@click.argument('datafile')
def plot_fig(datafile, columns, del_header, labels, xlim, ylim,
             title, figname, sci, equal, show, figtype):
    """Read two columns of data from DATAFILE and plot."""
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)
    x, y = read_two_col(datafile, columns[0], columns[1], del_header)
    pl.rc('font', family='serif')
    pl.figure(figsize=(8, 6))
    pl.plot(x, y, marker='.', markersize=2.0, color='r',
            linestyle='None')
    pl.xlabel(labels[0], fontsize=16)
    pl.ylabel(labels[1], fontsize=16)
    if not (xlim[0] == xlim[1]):
        pl.xlim(xlim)
    if not (ylim[0] == ylim[1]):
        pl.ylim(ylim)
    if sci:
        pl.ticklabel_format(style='sci', scilimits=(0, 0))
    if equal:
        pl.axis('equal')
    if len(title) > 0:
        pl.title(title, fontsize=16)
    fig_name = figname + '.' + figtype
    pl.savefig(fig_name)
    if show:
        pl.show()


if __name__ == '__main__':
    plot_fig()
