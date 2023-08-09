#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:30:26 2023

@author: js2746
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from utils import *
import os
from PIL import Image


def make_animated_heatmap(data, coordsys, dims, Vmax, Vmin, colorbar):
    _, _, _, _, _, dim1vals, dim2vals = dims
    fig = plt.figure()
    if coordsys == "polar":
        ax = plt.subplot(projection="polar")
    else:
        ax = plt.subplot()
    plt.axis('off')
    label = ax.text(1.5 * np.pi, 200, 'Frame 0', ha='center', va='center')
    c = ax.pcolormesh(dim2vals, dim1vals, data[:, :, 0], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    plt.colorbar(c)

    def animate(frame):
        c.set_array(data[:, :, frame])
        label.set_text('Frame ' + str(frame))
        return c

    anim = animation.FuncAnimation(fig, animate, frames=np.shape(data)[2], interval=333, repeat=False)
    anim.save("/home/js2746/PolarHeightBinning/plotting/test.gif", progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
    plt.close()


def make_animated_heatmap_with_avgovertheta(heatmap_data, coordsys, dims, Vmax, Vmin, colorbar):
    Nrbins, _, _, _, _, dim1vals, dim2vals = dims
    fig = plt.figure(layout="constrained")
    fig.set_figwidth(15)
    ax1 = plt.subplot(132, projection="polar")
    ax2 = plt.subplot(133)
    ax3 = plt.subplot(131)
    label1 = ax1.text(1.5 * np.pi, 250, 'Frame 0', ha='center', va='center')
    ax1.axis('off')
    ax3.axis('off')
    c = ax1.pcolormesh(dim2vals, dim1vals, heatmap_data[:, :, 0], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    plt.colorbar(c)
    avgovertheta_data = np.nanmean(heatmap_data, axis=1)
    avg_avgovertheta = np.nanmean(avgovertheta_data, axis=1)
    std_data = 2 * np.nanstd(avgovertheta_data, axis=1)
    X = np.linspace(0, Nrbins, Nrbins)
    e, = ax2.plot(X, avgovertheta_data[:, 0])
    avgavg = ax2.plot(X, avg_avgovertheta)
    std = ax2.fill_between(X, (avg_avgovertheta - std_data), (avg_avgovertheta + std_data), alpha=.1)
    ax2.set_ylim(Vmin, Vmax)
    imgpath = '/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_42us/movie_stills/lgPO.0000'
    img = np.asarray(Image.open(imgpath + "0.ppm"))
    f = ax3.imshow(img)

    def animate(frame):
        c.set_array(heatmap_data[:, :, frame])
        label1.set_text('Frame ' + str(frame))
        e.set_data(X, avgovertheta_data[:, frame])
        img = np.asarray(Image.open(imgpath + str(frame) + ".ppm"))
        f.set_array(img)
        return c, e, f

    anim = animation.FuncAnimation(fig, animate, frames=10, interval=333, repeat=False)
    anim.save("/home/js2746/PolarHeightBinning/plotting/test.gif", progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
    plt.close()


os.chdir("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_30us/lgPO_polar_5_10_0_-1_1")
hmap_data = np.load("npy/lgPO.zone.C1A.C1B.polar.meancurvature.npy")
dims = bin_prep("lgPO", "C1A.C1B", "polar", "OFF")
make_animated_heatmap_with_avgovertheta(hmap_data, "polar", dims, .3, -.3, True)
os.chdir("/home/js2746/PolarHeightBinning/plotting")
