import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os
from utils import *


def calculate_density(sys_name, names_dict, coordsys, inclusion, polar, dims, scale_dict):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims
    areas = np.load('npy/' + sys_name + "." + coordsys + ".areas.npy")

    for species in names_dict['species_list']:
        for leaflet in ["zone", "ztwo"]:
            data = np.genfromtxt('tcl_output/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.density.dat', missing_values='nan', filling_values="0")

            # create a new array that has each frame in a different array level
            density_array = np.zeros((N1_bins, N2_bins, Nframes))

            for frm in range(Nframes):
                density_array[:, :, frm] = data[frm * N1_bins:(frm + 1) * N1_bins, 2:]

            avgdensity = calc_avg_over_time(density_array)

            # normalize
            normfactor = names_dict["density_norm"]
            avgdensity = avgdensity * normfactor / areas

            # make plots!
            if leaflet == "zone":
                plot_maker(dim1vals, dim2vals, avgdensity, sys_name, species + '.outer', scale_dict["density_max"], scale_dict["density_min"], inclusion, "avgDensity", False, coordsys, scale_dict)
            elif leaflet == "ztwo":
                plot_maker(dim1vals, dim2vals, avgdensity, sys_name, species + '.inner', scale_dict["density_max"], scale_dict["density_min"], inclusion, "avgDensity", False, coordsys, scale_dict)

            # save as file for debugging / analysis
            np.save('npy/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.density.npy', density_array)
            np.save('npy/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.avgdensity.npy', avgdensity)
            if polar is True:
                avg_over_theta('npy/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.avgdensity')
            np.savetxt('dat/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.avgdensity.dat', avgdensity, delimiter=',', fmt='%10.5f')

        print(sys_name + ' ' + species + " density done!")

    # if len(species_list) > 1:
    #  calculate_total_density(sys_name, names_dict['species_list'], coordsys, inclusion, polar, dims)
    # elif len(species_list) < 1:
    #  print("Something is wrong with species_list!")


def calculate_total_density(sys_name, species_list, coordsys, inclusion, polar, dims):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    for leaflet in ['zone', 'ztwo']:
        tot_density = np.zeros((N1_bins, N2_bins, Nframes))

        for species in species_list:
            dens_per_species = np.load('npy/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.density.npy')
            tot_density = tot_density + dens_per_species
            # normalize?

        tot_avg_density = np.mean(tot_density, axis=2)

        # make plots!
        plot_maker(dim1vals, dim2vals, tot_avg_density, sys_name, species + '.' + leaflet, density_max, density_min, inclusion, "totAvgDensity", False, polar, scale_dict["colorbar"])

        # save as file for debugging / analysis
        np.save('npy/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.totalDensity.npy', tot_density)
        np.savetxt('dat/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.totalAvgDensity.dat', tot_avg_density, delimiter=',', fmt='%10.5f')

    print("Total density done!")
