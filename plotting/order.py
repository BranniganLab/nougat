"""Functions related to calculating order parameters."""
import numpy as np
from pathlib import Path
from utils import *


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
