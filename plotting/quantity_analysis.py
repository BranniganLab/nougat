"""Functions related to calculating membrane quantities."""
import numpy as np
from pathlib import Path
from utils import *


def calculate_density(system_dict, polar, cwd):
    """
    Calculate density enrichment for each species in the system and save to files.

    Parameters
    ----------
    system_dict : dict
        Dictionary containing list of species, bin sizes, and normalization \
            constants.
    polar : bool
        Whether or not to use polar coordinates instead of cartesian.
    cwd : Pathlib Path object
        Path to current working directory.

    Returns
    -------
    None.

    """
    areas = np.load(cwd.joinpath("trajectory", "density", "areas.npy"))

    for species in system_dict['species']:

        # make species-specific folders for density
        for folder in ["trajectory", "average"]:
            species_dir = cwd.joinpath(folder, "density", species)
            species_dir.mkdir(parents=True, exist_ok=True)

        for leaflet in ["zone", "ztwo"]:
            counts_array = parse_dat_file(cwd.joinpath("tcl_output", "density", species, leaflet + ".dat"), system_dict["bin_info"])

            density_array = counts_array / areas

            avgdensity = calc_avg_over_time(density_array)

            # normalize
            normfactor = system_dict["density_norm"][species]
            avgdensity = avgdensity * normfactor

            # save as file for debugging / analysis
            np.save(cwd.joinpath("trajectory", "density", species, leaflet + ".npy"), density_array)
            np.save(cwd.joinpath("average", "density", species, leaflet + ".npy"), avgdensity)
            if polar:
                avg_over_theta(cwd.joinpath("average", "density", species, leaflet))
            np.savetxt(cwd.joinpath("average", "density", species, leaflet + ".dat"), avgdensity, delimiter=',', fmt='%10.5f')

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
    dim1vals, dim2vals = dims
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


def calculate_order(system_dict, polar, cwd):
    """
    Parse order from .dat file and calculate time average before saving to file.

    Parameters
    ----------
    system_dict : dict
        Dict containing bin sizes and number of frames in trajectory.
    polar : bool
        Whether to use polar coordinates instead of cartesian.
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
                for folder in ["trajectory", "average"]:
                    dirname = cwd.joinpath(folder, "order", species, "tail" + str(tail))
                    dirname.mkdir(parents=True, exist_ok=True)

                order_array = parse_dat_file(cwd.joinpath("tcl_output", "order", species, "tail" + str(tail), leaflet + ".dat"), system_dict["bin_info"])

                order_array_pruned = mostly_empty(order_array)

                avgorder = calc_avg_over_time(order_array_pruned)

                # save as file for debugging / analysis
                np.save(cwd.joinpath("trajectory", "order", species, "tail" + str(tail), leaflet + ".npy"), order_array_pruned)
                np.save(cwd.joinpath("average", "order", species, "tail" + str(tail), leaflet + ".npy"), avgorder)
                if polar:
                    avg_over_theta(cwd.joinpath("average", "order", species, "tail" + str(tail), leaflet))
                np.savetxt(cwd.joinpath("average", "order", species, "tail" + str(tail), leaflet + ".dat"), avgorder, delimiter=',', fmt='%10.5f')

            print(species + " tail" + str(tail) + " order done!")
