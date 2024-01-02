"""Functions related to calculating density."""
import numpy as np
from pathlib import Path
from utils import *


def calculate_density(sys_name, names_dict, coordsys, inclusion, dims, scale_dict, cwd):
    N1_bins = names_dict['bin_info']['N1']
    N2_bins = names_dict['bin_info']['N2']
    Nframes = names_dict['bin_info']['nframes']
    dim1vals, dim2vals = dims
    areas = np.load(cwd.joinpath("trajectory", "density", "areas.npy"))

    for species in names_dict['species']:
        for folder in ["trajectory", "average"]:
            species_dir = cwd.joinpath(folder, "density", species)
            species_dir.mkdir(parents = True, exist_ok = True)
        for leaflet in ["zone", "ztwo"]:
            data = np.genfromtxt(cwd.joinpath("tcl_output", "density", species, leaflet + ".dat"), missing_values='nan', filling_values="0")

            # create a new array that has each frame in a different array level
            density_array = np.zeros((Nframes, N1_bins, N2_bins))

            for frm in range(Nframes):
                density_array[frm, :, :] = data[frm * N1_bins:(frm + 1) * N1_bins, 2:]

            avgdensity = calc_avg_over_time(density_array)

            # normalize
            normfactor = names_dict["density_norm"][species]
            avgdensity = avgdensity * normfactor / areas

            # save as file for debugging / analysis
            np.save(cwd.joinpath("trajectory", "density", species, leaflet + ".npy"), density_array)
            np.save(cwd.joinpath("average", "density", species, leaflet + ".npy"), avgdensity)
            if coordsys == "polar":
                avg_over_theta(cwd.joinpath("average", "density", species, leaflet))
            np.savetxt(cwd.joinpath("average", "density", species, leaflet + ".dat"), avgdensity, delimiter=',', fmt='%10.5f')

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

        tot_avg_density = np.mean(tot_density, axis=0)

        # make plots!
        plot_maker(dim1vals, dim2vals, tot_avg_density, sys_name, species + '.' + leaflet, density_max, density_min, inclusion, "totAvgDensity", False, polar, scale_dict["colorbar"])

        # save as file for debugging / analysis
        np.save('npy/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.totalDensity.npy', tot_density)
        np.savetxt('dat/' + sys_name + '.' + species + '.' + leaflet + '.' + coordsys + '.totalAvgDensity.dat', tot_avg_density, delimiter=',', fmt='%10.5f')

    print("Total density done!")
