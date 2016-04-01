#!/usr/bin/env python

"""
Script to plot Hill surface.
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def Hill_Surf(n, Miu, Cj):
    """Implicit equation of Hill surface."""
    def hill_surf(x, y, z):
        r1 = ((x + Miu)**2 + y**2 + z**2)**0.5
        r2 = ((x + Miu - 1)**2 + y**2 + z**2)**0.5
        return 0.5 * n**2 * (x**2 + y**2) + (1 - Miu)/r1 + Miu/r2 - 0.5 * Cj

    return hill_surf


def Hill_Surf_Cj_xy(n, Miu):
    """Cj function when z==0."""
    def func(x, y):
        r1 = ((x + Miu)**2 + y**2)**0.5
        r2 = ((x + Miu - 1)**2 + y**2)**0.5
        return n**2 * (x**2 + y**2) + 2 * (1 - Miu) / r1 + 2 * Miu / r2

    return func


def plot_implicit(ax, fn, bbox=(-2.5, 2.5)):
    """ create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval"""
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    A = np.linspace(xmin, xmax, 100)  # resolution of the contour
    B = np.linspace(xmin, xmax, 15)   # number of slices
    A1, A2 = np.meshgrid(A, A)  # grid on which the contour is plotted
    for z in B:  # plot contours in the XY plane
        X, Y = A1, A2
        Z = fn(X, Y, z)
        ax.contour(X, Y, Z+z, [z], zdir='z')
        # [z] defines the only level to plot for this contour for
        # this value of z

    for y in B:  # plot contours in the XZ plane
        X, Z = A1, A2
        Y = fn(X, y, Z)
        ax.contour(X, Y+y, Z, [y], zdir='y')

    for x in B:  # plot contours in the YZ plane
        Y, Z = A1, A2
        X = fn(x, Y, Z)
        ax.contour(X+x, Y, Z, [x], zdir='x')

    ax.tick_params(labelsize=6)  # smaller label size

    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends
    # the plot limits to encompass all values in the contour.
    ax.set_zlim3d(zmin, zmax)
    ax.set_xlim3d(xmin, xmax)
    ax.set_ylim3d(ymin, ymax)


def plot_contour(ax, f, bbox=(-2.5, 2.5), levels=[0]):
    """Plot contour."""
    A = np.linspace(bbox[0], bbox[1], 300)  # resolution of the contour
    A1, A2 = np.meshgrid(A, A)  # grid on which the contour is plotted
    z = f(A1, A2)
    if len(levels) > 1:
        ax.contour(A1, A2, z, levels, linewidths=0.5, colors='k')
        CF = ax.contourf(A1, A2, z, levels, cmap=plt.cm.jet)
        cbar = plt.colorbar(CF, shrink=0.8)
        # cbar.set_label('$C_{J}$', size=14)
        cbar.ax.tick_params(labelsize=8)
    else:
        ax.contour(A1, A2, z, levels)
    ax.tick_params(labelsize=8)


# parameters of Hill surface
Hill_Surf_Miu = 0.1
Hill_Surf_n = 1.0
Hill_Surf_Cj = np.array([3.5970, 3.4667, 3.0996])

# use serif font
plt.rc('font', family='serif')

# plot 3D Hill surface
fig1 = plt.figure(figsize=(8, 10))
ax1 = []
ax1.append(fig1.add_subplot(321, projection='3d'))
ax1.append(fig1.add_subplot(322))
ax1.append(fig1.add_subplot(323, projection='3d'))
ax1.append(fig1.add_subplot(324))
ax1.append(fig1.add_subplot(325, projection='3d'))
ax1.append(fig1.add_subplot(326))

for i in range(3):
    plot_implicit(ax1[2 * i], Hill_Surf(Hill_Surf_n, Hill_Surf_Miu,
                                        Hill_Surf_Cj[i]), bbox=(-2.0, 2.0))
    ax1[2 * i].set_title('$C_{L_{%d}}$=%.4f' % (i + 1, Hill_Surf_Cj[i]),
                         fontsize=10)
    plot_contour(ax1[2 * i + 1], Hill_Surf_Cj_xy(Hill_Surf_n, Hill_Surf_Miu),
                 bbox=(-2.0, 2.0), levels=[Hill_Surf_Cj[i]])

fig1.suptitle('Hill Surface, $\mu{}=0.1$', fontsize=16)
fig1.savefig('hill_surf.png')


# plot C_{J} contour
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
steps = 12
levels = np.linspace(2.9, 4.0, steps)

plot_contour(ax2, Hill_Surf_Cj_xy(Hill_Surf_n, Hill_Surf_Miu),
             bbox=(-2.0, 2.0), levels=levels)
ax2.set_xlabel('$x$', fontsize=14)
ax2.set_ylabel('$y$', fontsize=14)
ax2.set_title('$C_{J}$ Contour, $\mu{}=0.1$', fontsize=16)
fig2.savefig('cj_contour.png')


# plot C_{J} surface
fig3 = plt.figure()
ax3 = fig3.add_subplot(111, projection='3d')
bbox = (-2.0, 2.0)
A = np.linspace(bbox[0], bbox[1], 200)
X, Y = np.meshgrid(A, A)
Z = np.log10((Hill_Surf_Cj_xy(Hill_Surf_n, Hill_Surf_Miu))(X, Y))
surf = ax3.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=plt.cm.coolwarm,
                        linewidth=0, antialiased=False)
cbar = plt.colorbar(surf, shrink=0.5, aspect=8)
cbar.ax.tick_params(labelsize=8)
ax3.set_xlabel('$x$')
ax3.set_ylabel('$y$')
ax3.set_zlabel('$log(C_{J})$')
ax3.tick_params(labelsize=8)
ax3.set_title('$C_{J}$ Surface, $\mu{}=0.1, z=0$', fontsize=16)
ax3.dist = 12
fig3.savefig('cj_surface.png')

# show figures
plt.show()
