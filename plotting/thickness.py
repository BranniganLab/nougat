"""Functions related to calculating thickness."""
import numpy as np
from pathlib import Path
from utils import *


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
