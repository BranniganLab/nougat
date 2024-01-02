#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 13:44:55 2023.

@author: js2746
"""
import numpy as np
import matplotlib.pyplot as plt


color_dict = {
    "position": "red",
    "COMtilt": "blue"
}


def plot_tilt_colvar(list_of_paths, list_of_systems):
    fig, ax = plt.subplots()
    for file, sys in zip(list_of_paths, list_of_systems):
        data = np.genfromtxt(file)
        x = data[:, 0]
        y = data[:, 1]
        print(np.mean(y))
        tilt_angle = np.arccos(y) * (180 / np.pi)
        ax.plot(x, tilt_angle, label=sys)
    ax.set_ylabel("Tilt (degrees)")
    ax.set_xlabel("Frames")
    fig.legend()
    plt.savefig("/home/js2746/5x29_stiffness/COMtilt75k/tilt_colvar.pdf", dpi=700)
    plt.clf()
    plt.close()


position = "/home/js2746/5x29_stiffness/COMtilt75k/position.dat"
e75k = "/home/js2746/5x29_stiffness/COMtilt75k/tilt75k.dat"
e2k = "/home/js2746/5x29_stiffness/COMtilt75k/tilt2k.dat"
e20k = "/home/js2746/5x29_stiffness/COMtilt75k/tilt20k.dat"
e1MM = "/home/js2746/5x29_stiffness/COMtilt75k/tilt1MM.dat"
e1BB = "/home/js2746/5x29_stiffness/COMtilt1BB/tilt1BB.dat"
plot_tilt_colvar([position, e1BB, e1MM, e75k, e2k], ["position restraint", "1BB tilt restraint", "1MM tilt restraint", "75k tilt restraint", "2k tilt restraint"])
