"""Functions related to calculating membrane quantities."""
import numpy as np
from pathlib import Path
from utils import *


def calculate_density(polar, system_dict, cwd):
    """
    Calculate density enrichment for each species in the system and save to files.

    Parameters
    ----------
    polar : bool
        Whether or not to use polar coordinates instead of cartesian.
    system_dict : dict
        Dictionary containing list of species, bin sizes, and normalization \
            constants.
    cwd : Pathlib Path object
        Path to current working directory.

    Returns
    -------
    None.

    """
    areas = np.load(cwd.joinpath("trajectory", "density", "areas.npy"))

    for species in system_dict['species']:

        # make species-specific folders for density
        for folder in ["trajectory", "average", "figures"]:
            species_dir = cwd.joinpath(folder, "density", species)
            species_dir.mkdir(parents=True, exist_ok=True)

        for leaflet in ["zone", "ztwo"]:
            counts_array = parse_dat_file(cwd.joinpath("tcl_output", "density", species, leaflet + ".dat"), system_dict["bin_info"], "density")

            density_array = counts_array / areas

            # normalize
            normfactor = system_dict["density_norm"][species]
            density_enrichment = density_array * normfactor

            avgdensity_enrichment = calc_avg_over_time(density_enrichment)

            # save as file for debugging / analysis
            np.save(cwd.joinpath("trajectory", "density", species, leaflet + ".npy"), density_enrichment)
            np.save(cwd.joinpath("average", "density", species, leaflet + ".npy"), avgdensity_enrichment)
            if polar:
                avg_over_theta(cwd.joinpath("average", "density", species, leaflet))
            np.savetxt(cwd.joinpath("average", "density", species, leaflet + ".dat"), avgdensity_enrichment, delimiter=',', fmt='%10.5f')

        print(species + " density done!")


def calculate_thickness(polar, cwd):
    """
    Calculate leaflet thickness from surface heights and save to files.

    Parameters
    ----------
    polar : bool
        Whether or not to use polar coordinates.
    cwd : pathlib Path object
        Current working directory.

    Returns
    -------
    None.

    """
    zone = np.load(cwd.joinpath("trajectory", "height", "zone.npy"))
    ztwo = np.load(cwd.joinpath("trajectory", "height", "ztwo.npy"))
    zzero = np.load(cwd.joinpath("trajectory", "height", "zzero.npy"))

    for field in ["zone", "ztwo", "whole"]:
        if field == "zone":
            thickness = zone - zzero
        elif field == "ztwo":
            thickness = zzero - ztwo
        elif field == "whole":
            thickness = zone - ztwo

        avgthickness = calc_avg_over_time(thickness)

        # save as file for debugging / analysis!
        np.save(cwd.joinpath("trajectory", "thickness", field + ".npy"), thickness)
        np.save(cwd.joinpath("average", "thickness", field + ".npy"), avgthickness)
        if polar:
            avg_over_theta(cwd.joinpath("average", "thickness", field))
        np.savetxt(cwd.joinpath("average", "thickness", field + ".dat"), avgthickness, delimiter=',', fmt='%10.5f')

    print("Thickness done!")


def calculate_order(polar, system_dict, cwd):
    """
    Parse order from .dat file and calculate time average before saving to file.

    Parameters
    ----------
    polar : bool
        Whether to use polar coordinates instead of cartesian.
    system_dict : dict
        Dict containing bin sizes and number of frames in trajectory.
    cwd : Pathlib Path object
        Path to current working directory.

    Returns
    -------
    None.

    """
    for species in system_dict['species']:
        for tail in range(system_dict['ntails'][species]):
            for leaflet in ["zone", "ztwo"]:

                # make species- and tail-specific folders for outputs
                for folder in ["trajectory", "average", "figures"]:
                    dirname = cwd.joinpath(folder, "order", species, "tail" + str(tail))
                    dirname.mkdir(parents=True, exist_ok=True)

                order_array = parse_dat_file(cwd.joinpath("tcl_output", "order", species, "tail" + str(tail), leaflet + ".dat"), system_dict["bin_info"], "order")

                order_array_pruned = mostly_empty(order_array)

                avgorder = calc_avg_over_time(order_array_pruned)

                # save as file for debugging / analysis
                np.save(cwd.joinpath("trajectory", "order", species, "tail" + str(tail), leaflet + ".npy"), order_array_pruned)
                np.save(cwd.joinpath("average", "order", species, "tail" + str(tail), leaflet + ".npy"), avgorder)
                if polar:
                    avg_over_theta(cwd.joinpath("average", "order", species, "tail" + str(tail), leaflet))
                np.savetxt(cwd.joinpath("average", "order", species, "tail" + str(tail), leaflet + ".dat"), avgorder, delimiter=',', fmt='%10.5f')

            print(species + " tail" + str(tail) + " order done!")


# still need to convert to pathlib, remove coordsys!
def calc_elastic_terms(path, polar, bin_info):
    """
    Calculate all the additional terms that appear in a hamiltonian or are \
        generally of interest.

    Parameters
    ----------
    path : string
        should point to the folder housing your nougat outputs for the given \
            system
    polar : bool
        if True, use polar (cylindrical) coordinates

    Returns
    -------
    None.

    """
    # load height and curvature data
    z_1 = np.load(path.joinpath("trajectory", "height", "zone.npy"))
    z_2 = np.load(path.joinpath("trajectory", "height", "ztwo.npy"))
    z_0 = np.load(path.joinpath("trajectory", "height", "zzero.npy"))
    z_plus = np.load(path.joinpath("trajectory", "height", "zplus.npy"))
    H_1 = np.load(path.joinpath("trajectory", "curvature", "mean", "zone.npy"))
    H_2 = np.load(path.joinpath("trajectory", "curvature", "mean", "ztwo.npy"))
    K_1 = np.load(path.joinpath("trajectory", "curvature", "gaussian", "zone.npy"))
    K_2 = np.load(path.joinpath("trajectory", "curvature", "gaussian", "ztwo.npy"))

    # measure terms of interest
    # removed z_minus terms until we have a better way of computing t0
    epsilon = z_0 - z_plus
    epsilon2 = epsilon**2
    H_plus = (H_1 + H_2) / 2
    K_plus = (K_1 + K_2) / 2
    K_minus = (K_1 - K_2) / 2
    H_plus2 = H_plus**2
    H_minus = (H_1 - H_2) / 2
    H_minus2 = H_minus**2
    epsilon_H = epsilon * H_plus

    # save useful trajectories
    np.save(path.joinpath("trajectory", "thickness", "epsilon.npy"), epsilon)
    np.save(path.joinpath("trajectory", "curvature", "mean", "Hplus.npy"), H_plus)
    np.save(path.joinpath("trajectory", "thickness", "epsilon2.npy"), epsilon2)
    np.save(path.joinpath("trajectory", "curvature", "mean", "Hplus2.npy"), H_plus2)

    # calculate averages
    avg_epsilon = calc_avg_over_time(epsilon)
    avg_epsilon2 = calc_avg_over_time(epsilon2)
    avg_H_plus = calc_avg_over_time(H_plus)
    avg_H_plus2 = calc_avg_over_time(H_plus2)
    avg_H_minus = calc_avg_over_time(H_minus)
    avg_H_minus2 = calc_avg_over_time(H_minus2)
    avg_epsilon_H = calc_avg_over_time(epsilon_H)
    avg_K_plus = calc_avg_over_time(K_plus)
    avg_K_minus = calc_avg_over_time(K_minus)

    # calculate correlations
    corr_eps_Hplus = calc_avg_over_time(epsilon * H_plus) - (avg_epsilon * avg_H_plus)
    corr_mag_eps_Hplus = calc_avg_over_time(np.sqrt(epsilon2) * np.sqrt(H_plus2)) - (np.sqrt(avg_epsilon2) * np.sqrt(avg_H_plus2))
    corr_eps_Kplus = calc_avg_over_time(epsilon * K_plus) - (avg_epsilon * avg_K_plus)

    # get proper plot dimensions
    dims = bin_prep(bin_info, polar)
    dim1vals, dim2vals = dims

    # make pretty pictures and save data
    data_list = [avg_K_plus, avg_K_minus, corr_eps_Kplus, corr_mag_eps_Hplus,
                 corr_eps_Hplus, avg_epsilon, avg_epsilon2, avg_H_plus,
                 avg_H_plus2, avg_H_minus, avg_H_minus2, avg_epsilon_H]
    name_list = ["avg_K_plus", "avg_K_minus", "corr_eps_Kplus",
                 "corr_mag_eps_Hplus", "corr_eps_Hplus", "avg_epsilon",
                 "avg_epsilon2", "avg_H_plus", "avg_H_plus2", "avg_H_minus",
                 "avg_H_minus2", "avg_epsilon_H"]
    for data, name in zip(data_list, name_list):
        # plot_maker(dim1vals, dim2vals, data, system, 'comb', .1, -.1, False, name, False, coordsys, scale_dict)
        np.save(path.joinpath("average", "misc", name + '.npy'), data)
        if polar:
            avg_over_theta(path.joinpath("average", "misc", name))
