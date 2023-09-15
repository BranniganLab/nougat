"""Functions related to calculating order parameters."""
import numpy as np
from utils import *


def calculate_order(sys_name, system_dict, coordsys, inclusion, polar, dims, scale_dict):
    N1_bins = system_dict['bin_info']['N1']
    N2_bins = system_dict['bin_info']['N2']
    Nframes = system_dict['bin_info']['nframes']
    dim1vals, dim2vals = dims

    for species in system_dict['species']:
        for tail in range(system_dict['ntails'][species]):
            for leaflet in ["zone", "ztwo"]:
                order_file = np.genfromtxt('tcl_output/' + sys_name + '.' + species + '.tail' + str(tail) + '.' + leaflet + '.' + coordsys + '.order.dat', missing_values='nan', filling_values=np.nan)

                # create a new array that has each frame in a different array level
                order_array = np.zeros((N1_bins, N2_bins, Nframes))
                for frm in range(Nframes):
                    order_array[:, :, frm] = order_file[frm * N1_bins:(frm + 1) * N1_bins, 2:]

                order_array_pruned = mostly_empty(order_array)

                avgorder = calc_avg_over_time(order_array_pruned)

                # make plots!
                plot_maker(dim1vals, dim2vals, avgorder, sys_name, species + '.tail' + str(tail) + '.zone''.' + leaflet, scale_dict["order_max"], scale_dict["order_min"], inclusion, "avgOrder", False, coordsys, scale_dict)

                # save as file for debugging / analysis
                np.save('npy/' + sys_name + '.' + species + '.tail' + str(tail) + '.' + leaflet + '.' + coordsys + '.order.npy', order_array_pruned)
                np.save('npy/' + sys_name + '.' + species + '.tail' + str(tail) + '.' + leaflet + '.' + coordsys + '.avgorder.npy', avgorder)
                if coordsys == "polar":
                    avg_over_theta('npy/' + sys_name + '.' + species + '.tail' + str(tail) + '.' + leaflet + '.' + coordsys + '.avgorder')
                np.savetxt('dat/' + sys_name + '.' + species + '.tail' + str(tail) + '.' + leaflet + '.' + coordsys + '.avgOrder.dat', avgorder, delimiter=',', fmt='%10.5f')

            print(sys_name + ' ' + species + " tail" + str(tail) + " order done!")

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
