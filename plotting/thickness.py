"""Functions related to calculating thickness."""
import numpy as np
from pathlib import Path
from utils import *


def calculate_thickness(sys_name, coordsys, inclusion, dims, scale_dict, cwd):
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

        # need to fix this to get t0 from an empty membrane
        # avgt0 = measure_t0(zone, ztwo, coordsys)
        # normthickness = avgthickness/avgt0

        # make plots!
        # if leaflet == "whole":
        # plot_maker(dim1vals, dim2vals, avgthickness, sys_name, leaflet, scale_dict["thick_whole_max"], scale_dict["thick_whole_min"], inclusion, "avgThickness", bead, coordsys, scale_dict)
        # else:
        # plot_maker(dim1vals, dim2vals, avgthickness, sys_name, leaflet, scale_dict["thick_max"], scale_dict["thick_min"], inclusion, "avgThickness", bead, coordsys, scale_dict)

        # save as file for debugging / analysis!
        np.save(cwd.joinpath("trajectory", "thickness", field + ".npy"), thickness)
        np.save(cwd.joinpath("average", "thickness", field + ".npy"), avgthickness)
        if coordsys == "polar":
            avg_over_theta(cwd.joinpath("average", "thickness", field))
        np.savetxt(cwd.joinpath("average", "thickness", field + ".dat"), avgthickness, delimiter=',', fmt='%10.5f')

    print(sys_name + " thickness done!")
