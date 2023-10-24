#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 11:30:26 2023.

@author: js2746
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from utils import bin_prep, gifformat, read_log, calc_avg_over_time
import os
from PIL import Image


def make_animated_heatmap(data, coordsys, dims, Vmax, Vmin, colorbar):
    _, _, _, _, _, dim1vals, dim2vals = dims
    fig = plt.figure()
    if coordsys == "polar":
        ax = plt.subplot(projection="polar")
        label = ax.text(1.5 * np.pi, 200, 'Frame 0', ha='center', va='center')
        c = ax.pcolormesh(dim2vals, dim1vals, data[0, :, :], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    else:
        ax = plt.subplot()
        label = ax.text(250, -1, 'Frame 0', ha='center', va='center')
        c = ax.pcolormesh(dim1vals, dim2vals, data[0, :, :], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    plt.axis('off')
    plt.colorbar(c)

    def animate(frame):
        c.set_array(data[frame, :, :])
        label.set_text('Frame ' + str(frame))
        return c

    anim = animation.FuncAnimation(fig, animate, frames=np.shape(data)[0], interval=100, repeat=False)
    anim.save("/home/js2746/PolarHeightBinning/plotting/test.gif", progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
    plt.close()


def make_animated_heatmap_with_avgovertheta(heatmap_data, coordsys, dims, Vmax, Vmin, colorbar, bin_info, imgpath):
    dim1vals, dim2vals = dims
    Nrbins = bin_info["N1"]
    dr = bin_info["d1"]
    fig = plt.figure(layout="constrained")
    fig.set_figwidth(15)

    # average over theta on right
    avgovertheta_data = np.nanmean(heatmap_data, axis=2)
    if Vmax == "auto":
        Vmax = np.nanmax(np.abs(avgovertheta_data))
        Vmin = -1.0 * Vmax
    avgoverthetaovertime = np.nanmean(avgovertheta_data, axis=0)
    std_data = 2 * np.nanstd(avgovertheta_data, axis=0)
    xmin = (dr / 2.0)
    xmax = ((Nrbins * dr) - xmin)
    xmin = xmin / 10
    xmax = xmax / 10
    X = np.linspace(xmin, xmax, Nrbins)

    ax2 = plt.subplot(133)
    e, = ax2.plot(X, avgovertheta_data[0, :])  # the comma after e is super important!
    ax2.plot(X, avgoverthetaovertime)
    ax2.fill_between(X, (avgoverthetaovertime - std_data), (avgoverthetaovertime + std_data), alpha=.1)
    ax2.set_ylim(Vmin, Vmax)
    ax2.set_xlim(0, xmax)

    # polar 2d heatmap in middle
    ax1 = plt.subplot(132, projection="polar")
    label1 = ax1.text(1.5 * np.pi, 250, 'Frame 0', ha='center', va='center')
    ax1.axis('off')
    c = ax1.pcolormesh(dim2vals, dim1vals, heatmap_data[0, :, :], cmap="RdBu_r", vmax=Vmax, vmin=Vmin, zorder=0)
    plt.colorbar(c)

    # VMD still on left
    img = np.asarray(Image.open(imgpath + "00000.ppm"))
    ax3 = plt.subplot(131)
    ax3.axis('off')
    f = ax3.imshow(img)

    def animate(frame):
        c.set_array(heatmap_data[frame, :, :])
        label1.set_text('Frame ' + str(frame))
        e.set_data(X, avgovertheta_data[frame, :])
        img = np.asarray(Image.open(imgpath + gifformat(frame, 5) + ".ppm"))
        f.set_array(img)
        return c, e, f

    anim = animation.FuncAnimation(fig, animate, frames=np.shape(heatmap_data)[0] - 1, interval=100, repeat=False)
    # anim = animation.FuncAnimation(fig, animate, frames=3, interval=1000, repeat=False)
    anim.save("COMtiltspin_z2t1order_over_traj.gif", progress_callback=lambda i, n: print(f'Saving frame {i} of {n}'))
    plt.close()


def smooth(hmap, window_size):
    """
    Performs moving window average on 2nd axis (time).

    Parameters
    ----------
    hmap : numpy ndarray
        An array containing data that needs to be smoothed over time.
    window_size : int
        The window size for smoothing.

    Returns
    -------
    numpy ndarray
        The input array, smoothed over axis 2.

    """
    if window_size % 1 != 0:
        print("give me integer")
        return
    elif window_size < 1:
        print("give me positive integer")
        return
    elif window_size % 2 == 0:
        print("give me odd integer")
        return
    elif window_size == 1:
        print("window_size is set to 1, so no smoothing was done")
        return hmap

    middle = window_size // 2
    hmap_copy = hmap.copy()
    for frame in range(np.shape(hmap)[0] - 2 * middle):
        hmap_copy[frame + middle, :, :] = calc_avg_over_time(hmap[frame: frame + window_size, :, :])

    return hmap_copy


# os.chdir("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_42us/lgPO_polar_40_12_0_-1_1")
os.chdir("/home/js2746/Bending/PC/whole_mols/5x29/COMtiltspin/COMtiltspin_polar_10_10_0_-1_2")
bin_info_dict = read_log("COMtiltspin", "polar")
hmap_data = np.load("npy/COMtiltspin.POPC.tail1.ztwo.polar.order.npy")
print(np.shape(hmap_data))
# hmap_data = smooth(hmap_data, 11)
dims = bin_prep(bin_info_dict["bin_info"], "polar")
# imgpath = '/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_42us/movie_stills/lgPO.'
imgpath = '/home/js2746/Bending/PC/whole_mols/5x29/COMtiltspin/movie/COMtiltspin.'
make_animated_heatmap_with_avgovertheta(hmap_data, "polar", dims, "auto", 0, True, bin_info_dict["bin_info"], imgpath)
# make_animated_heatmap(hmap_data, "cart", dims, .025, -.025, True)
os.chdir("/home/js2746/PolarHeightBinning/plotting")
