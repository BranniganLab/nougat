#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:30:26 2023

@author: js2746
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter
from utils import *
import os


def make_animated_heatmap(data, coordsys, dim1vals, dim2vals, Vmax, Vmin, colorbar):
    fig = plt.figure()
    if coordsys == "polar":
        ax = plt.subplot(projection="polar")
    else:
        ax = plt.subplot()
    plt.axis('off')
    label = ax.text(1.5 * np.pi, 200, 'Frame 0', ha='center', va='center')

    def init():
        c = plt.pcolormesh(dim2vals, dim1vals, data[:, :, 0], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
        plt.colorbar(c)
        return ax

    def animate(frame):
        c = plt.pcolormesh(dim2vals, dim1vals, data[:, :, frame], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
        label.set_text('Frame ' + str(frame))
        return c

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=np.shape(data)[2], repeat=False)
    anim.save("test.gif", progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'), dpi=500, writer=PillowWriter(fps=3))
    plt.show()


os.chdir("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO/lgPO_polar_5_10_0_-1_1")
data = np.load("npy/lgPO.zone.C1A.C1B.polar.meancurvature.npy")
dims = bin_prep("lgPO", "C1A.C1B", "polar", "OFF")
_, _, _, _, _, dim1vals, dim2vals = dims
make_animated_heatmap(data, "polar", dim1vals, dim2vals, .3, -.3, True)
os.chdir("/home/js2746/PolarHeightBinning/plotting")
