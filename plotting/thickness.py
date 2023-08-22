"""Functions related to calculating thickness."""
import numpy as np
from utils import *


def calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims, scale_dict):
    dim1vals, dim2vals = dims
    zone = np.load('npy/' + sys_name + '.zone.' + bead + '.' + coordsys + '.height.npy')
    ztwo = np.load('npy/' + sys_name + '.ztwo.' + bead + '.' + coordsys + '.height.npy')
    zzero = np.load('npy/' + sys_name + '.zzero.' + bead + '.' + coordsys + '.height.npy')

    for leaflet in ["zone", "ztwo", "whole"]:
        if leaflet == "zone":
            thickness = zone - zzero
        elif leaflet == "ztwo":
            thickness = zzero - ztwo
        elif leaflet == "whole":
            thickness = zone - ztwo

        avgthickness = calc_avg_over_time(thickness)

        # need to fix this to get t0 from an empty membrane
        # avgt0 = measure_t0(zone, ztwo, coordsys)
        # normthickness = avgthickness/avgt0

        # make plots!
        if leaflet == "whole":
            plot_maker(dim1vals, dim2vals, avgthickness, sys_name, leaflet, scale_dict["thick_whole_max"], scale_dict["thick_whole_min"], inclusion, "avgThickness", bead, coordsys, scale_dict)
        else:
            plot_maker(dim1vals, dim2vals, avgthickness, sys_name, leaflet, scale_dict["thick_max"], scale_dict["thick_min"], inclusion, "avgThickness", bead, coordsys, scale_dict)

        # save as file for debugging / analysis!
        np.save('npy/' + sys_name + '.' + leaflet + '.' + bead + '.' + coordsys + '.thickness.npy', thickness)
        np.save('npy/' + sys_name + '.' + leaflet + '.' + bead + '.' + coordsys + '.avgthickness.npy', avgthickness)
        if coordsys == "polar":
            avg_over_theta('npy/' + sys_name + '.' + leaflet + '.' + bead + '.' + coordsys + '.avgthickness')
        np.savetxt('dat/' + sys_name + '.' + leaflet + '.' + bead + '.' + coordsys + '.avgthickness.dat', avgthickness, delimiter=',', fmt='%10.5f')

    print(sys_name + ' ' + bead + " thickness done!")
