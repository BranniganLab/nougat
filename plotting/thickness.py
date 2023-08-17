import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os
from utils import *


def calculate_thickness(sys_name, bead, coordsys, inclusion, polar, dims, scale_dict):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims
    zone = np.load('npy/' + sys_name + '.zone.' + bead + '.' + coordsys + '.height.npy')
    ztwo = np.load('npy/' + sys_name + '.ztwo.' + bead + '.' + coordsys + '.height.npy')
    zzero = np.load('npy/' + sys_name + '.zzero.' + bead + '.' + coordsys + '.height.npy')

    outer_leaflet = zone - zzero
    inner_leaflet = zzero - ztwo
    whole_membrane = zone - ztwo

    # need to fix this to get t0 from an empty membrane
    # avgt0 = measure_t0(zone, ztwo, coordsys)

    avgouter = calc_avg_over_time(outer_leaflet)
    avginner = calc_avg_over_time(inner_leaflet)
    avgwhole = calc_avg_over_time(whole_membrane)

    # normouter = avgouter/avgt0
    # norminner = avginner/avgt0
    # normwhole = avgwhole/avgt0

    # make plots!
    plot_maker(dim1vals, dim2vals, avgouter, sys_name, 'zone', scale_dict["thick_max"], scale_dict["thick_min"], inclusion, "avgThickness", bead, coordsys, scale_dict["colorbar"])
    plot_maker(dim1vals, dim2vals, avginner, sys_name, 'ztwo', scale_dict["thick_max"], scale_dict["thick_min"], inclusion, "avgThickness", bead, coordsys, scale_dict["colorbar"])
    plot_maker(dim1vals, dim2vals, avgwhole, sys_name, 'whole', scale_dict["thick_whole_max"], scale_dict["thick_whole_min"], inclusion, "avgThickness", bead, coordsys, scale_dict["colorbar"])
    # plot_maker(dim1vals, dim2vals, normouter, sys_name, 'zone', 1.2, 0, inclusion, "normThickness", bead, coordsys)
    # plot_maker(dim1vals, dim2vals, norminner, sys_name, 'ztwo', 1.2, 0, inclusion, "normThickness", bead, coordsys)
    # plot_maker(dim1vals, dim2vals, normwhole, sys_name, 'whole', 1.2, 0, inclusion, "normThickness", bead, coordsys)

    # save as file for debugging / analysis!
    np.save('npy/' + sys_name + '.zone.' + bead + '.' + coordsys + '.thickness.npy', outer_leaflet)
    np.save('npy/' + sys_name + '.ztwo.' + bead + '.' + coordsys + '.thickness.npy', inner_leaflet)
    np.save('npy/' + sys_name + '.whole.' + bead + '.' + coordsys + '.thickness.npy', whole_membrane)
    np.savetxt('dat/' + sys_name + '.zone.' + bead + '.' + coordsys + '.avgthickness.dat', avgouter, delimiter=',', fmt='%10.5f')
    np.savetxt('dat/' + sys_name + '.ztwo.' + bead + '.' + coordsys + '.avgthickness.dat', avginner, delimiter=',', fmt='%10.5f')
    np.savetxt('dat/' + sys_name + '.whole.' + bead + '.' + coordsys + '.avgthickness.dat', avgwhole, delimiter=',', fmt='%10.5f')
    # np.savetxt('dat/'+sys_name+'.zone.'+bead+'.'+coordsys+'.normthickness.dat', normouter,delimiter = ',',fmt='%10.5f')
    # np.savetxt('dat/'+sys_name+'.ztwo.'+bead+'.'+coordsys+'.normthickness.dat', norminner,delimiter = ',',fmt='%10.5f')
    # np.savetxt('dat/'+sys_name+'.whole.'+bead+'.'+coordsys+'.normthickness.dat', normwhole,delimiter = ',',fmt='%10.5f')

    print(sys_name + ' ' + bead + " thickness done!")
