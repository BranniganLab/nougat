#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:30:26 2023.

@author: js2746
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from utils import bin_prep, gifformat
import os
from PIL import Image


def make_animated_heatmap(data, coordsys, dims, Vmax, Vmin, colorbar):
    _, _, _, _, _, dim1vals, dim2vals = dims
    fig = plt.figure()
    if coordsys == "polar":
        ax = plt.subplot(projection="polar")
        label = ax.text(1.5 * np.pi, 200, 'Frame 0', ha='center', va='center')
        c = ax.pcolormesh(dim2vals, dim1vals, data[:, :, 0], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    else:
        ax = plt.subplot()
        label = ax.text(250, -1, 'Frame 0', ha='center', va='center')
        c = ax.pcolormesh(dim1vals, dim2vals, data[:, :, 0], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    plt.axis('off')
    plt.colorbar(c)

    def animate(frame):
        c.set_array(data[:, :, frame])
        label.set_text('Frame ' + str(frame))
        return c

    anim = animation.FuncAnimation(fig, animate, frames=np.shape(data)[2], interval=100, repeat=False)
    anim.save("/home/js2746/PolarHeightBinning/plotting/test.gif", progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
    plt.close()


def make_animated_heatmap_with_avgovertheta(heatmap_data, coordsys, dims, Vmax, Vmin, colorbar):
    Nrbins, dr, _, _, _, dim1vals, dim2vals = dims
    fig = plt.figure(layout="constrained")
    fig.set_figwidth(15)

    # average over theta on right
    avgovertheta_data = np.nanmean(heatmap_data, axis=1)
    if Vmax == "auto":
        Vmax = np.nanmax(np.abs(avgovertheta_data))
        Vmin = -1.0 * Vmax
    avgoverthetaovertime = np.nanmean(avgovertheta_data, axis=1)
    std_data = 2 * np.nanstd(avgovertheta_data, axis=1)
    xmin = (dr / 2.0) / 10
    xmax = ((Nrbins * dr) - xmin) / 10
    X = np.linspace(xmin, xmax, Nrbins)

    ax2 = plt.subplot(133)
    e, = ax2.plot(X, avgovertheta_data[:, 0])  # the comma after e is super important!
    ax2.plot(X, avgoverthetaovertime)
    ax2.fill_between(X, (avgoverthetaovertime - std_data), (avgoverthetaovertime + std_data), alpha=.1)
    ax2.set_ylim(Vmin, Vmax)
    ax2.set_xlim(0, xmax)

    # polar 2d heatmap in middle
    ax1 = plt.subplot(132, projection="polar")
    label1 = ax1.text(1.5 * np.pi, 250, 'Frame 0', ha='center', va='center')
    ax1.axis('off')
    c = ax1.pcolormesh(dim2vals, dim1vals, heatmap_data[:, :, 0], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    plt.colorbar(c)

    # VMD still on left
    imgpath = '/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_42us/movie_stills/lgPO.'
    img = np.asarray(Image.open(imgpath + "00000.ppm"))
    ax3 = plt.subplot(131)
    ax3.axis('off')
    f = ax3.imshow(img)

    def animate(frame):
        c.set_array(heatmap_data[:, :, frame])
        label1.set_text('Frame ' + str(frame))
        e.set_data(X, avgovertheta_data[:, frame])
        img = np.asarray(Image.open(imgpath + gifformat(frame, 5) + ".ppm"))
        f.set_array(img)
        return c, e, f

    anim = animation.FuncAnimation(fig, animate, frames=np.shape(heatmap_data)[2], interval=100, repeat=False)
    # anim = animation.FuncAnimation(fig, animate, frames=3, interval=1000, repeat=False)
    anim.save("lgPO_Hplus2_over_traj.gif", progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
    plt.close()


os.chdir("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_42us/lgPO_polar_40_12_0_-1_1")
hmap_data = np.load("npy/lgPO.H_plus2.npy")
dims = bin_prep("lgPO", "C1A.C1B", "polar", "OFF")
make_animated_heatmap_with_avgovertheta(hmap_data, "polar", dims, .0001, 0, True)
# make_animated_heatmap(hmap_data, "cart", dims, .025, -.025, True)
os.chdir("/home/js2746/PolarHeightBinning/plotting")
