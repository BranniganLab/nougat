import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os
from utils import *


def calculate_order(sys_name, names_dict, coordsys, inclusion, polar, dims, scale_dict):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    for species in names_dict['species_list']:
        for tail in names_dict[species]:
            zone = np.genfromtxt('tcl_output/' + sys_name + '.' + species + '.' + tail + '.zone.' + coordsys + '.order.dat', missing_values='nan', filling_values=np.nan)
            ztwo = np.genfromtxt('tcl_output/' + sys_name + '.' + species + '.' + tail + '.ztwo.' + coordsys + '.order.dat', missing_values='nan', filling_values=np.nan)

            # create a new array that has each frame in a different array level
            order_up = np.zeros((N1_bins, N2_bins, Nframes))
            order_down = np.zeros((N1_bins, N2_bins, Nframes))
            for frm in range(Nframes):
                order_up[:, :, frm] = zone[frm * N1_bins:(frm + 1) * N1_bins, 2:]
                order_down[:, :, frm] = ztwo[frm * N1_bins:(frm + 1) * N1_bins, 2:]

            order_up = mostly_empty(order_up)
            order_down = mostly_empty(order_down)

            avgouter = calc_avg_over_time(order_up)
            avginner = calc_avg_over_time(order_down)

            # make plots!
            plot_maker(dim1vals, dim2vals, avgouter, sys_name, species + '.' + tail + '.zone', scale_dict["order_max"], scale_dict["order_min"], inclusion, "avgOrder", False, coordsys)
            plot_maker(dim1vals, dim2vals, avginner, sys_name, species + '.' + tail + '.ztwo', scale_dict["order_max"], scale_dict["order_min"], inclusion, "avgOrder", False, coordsys)

            # save as file for debugging / analysis
            np.save('npy/' + sys_name + '.' + species + '.' + tail + '.zone.' + coordsys + '.order.npy', order_up)
            np.save('npy/' + sys_name + '.' + species + '.' + tail + '.ztwo.' + coordsys + '.order.npy', order_down)
            np.savetxt('dat/' + sys_name + '.' + species + '.' + tail + '.zone.' + coordsys + '.avgOrder.dat', avgouter, delimiter=',', fmt='%10.5f')
            np.savetxt('dat/' + sys_name + '.' + species + '.' + tail + '.ztwo.' + coordsys + '.avgOrder.dat', avginner, delimiter=',', fmt='%10.5f')

            print(sys_name + ' ' + species + " " + tail + " order done!")

        # if len(names_dict[species]) > 1:
        #  calculate_total_order(sys_name, species, names_dict, coordsys, inclusion, polar, dims)
        # elif len(names_dict[species]) < 1:
        #  print("Something is wrong with your tails list!")

    # if len(names_dict['species_list']) > 1:
    #  calculate_total_order(sys_name, "all", names_dict, coordsys, inclusion, polar, dims)
    # elif len(names_dict['species_list']) < 1:
    #  print("Something is wrong with your species_list!")


def calculate_total_order(sys_name, species, names_dict, coordsys, inclusion, polar, dims):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    if species != "all":
        for leaflet in ['zone', 'ztwo']:
            tot_order = np.zeros((N1_bins, N2_bins, Nframes))

            for tail in names_dict[species]:
                order_per_tail = np.load('npy/' + sys_name + '.' + species + '.' + tail + '.' + leaflet + '.' + coordsys + '.order.npy')
                tot_order = tot_order + order_per_tail
                # weight average by density!
