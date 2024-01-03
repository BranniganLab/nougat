"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from height import *
from curvature import *
from quantity_analysis import *
from tilt import *
from utils import *


def run_nougat(polar, inclusion_drawn):
    """
    Run nougat's averaging and image processing routines.

    Parameters
    ----------
    polar : boolean
        True for polar coordinate system; False for cartesian
    inclusion_drawn : boolean
        This feature currently isn't implemented, but would be how you acccount
        include an inclusion in your graphics

    Returns
    -------
    None.

    """
    cwd = Path.cwd()

    # make necessary folders
    create_outfile_directories(cwd)

    # define inclusion if present
    if inclusion_drawn is True:
        inclusion = add_inclusion(name, field_list)  # this proc doesn't exist yet!
    else:
        inclusion = False

    # figure out all the important info about the system you'll need
    system_dict = read_log()

    # prep heatmap plot dimensions
    hmap_dims = bin_prep(system_dict['bin_info'], polar)

    # read in height files and calculate surface heights
    # parse_height returns system_dict bc it adds nframes to the dictionary
    system_dict = parse_height(system_dict, polar, cwd)

    calculate_thickness(polar, cwd)
    calculate_curvature(polar, system_dict, cwd)

    save_areas(system_dict["bin_info"], 0, polar)

    # $$$$$$$$$$ UNTESTED FEATURES IN DEVELOPMENT BELOW $$$$$$$$$$$

    calculate_density(polar, system_dict, cwd)

    calculate_order(polar, system_dict, cwd)

    # calculate_tilt(sys_name, system_dict, coordsys, inclusion, cwd)

    # calc_elastic_terms(".", coordsys, config_dict, system_dict['bin_info'])

    if polar:
        coordsys = "polar"
    else:
        coordsys = "cart"

    for species in system_dict['species']:
        for field in ["zone", "ztwo", "zzero", "zplus"]:
            for quantity in ['height', 'curvature/gaussian', 'curvature/mean']:
                hmap_data = np.genfromtxt(cwd.joinpath("average", quantity, field + ".dat"), delimiter=",")
                fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, quantity, polar)
                plt.savefig(cwd.joinpath("figures", quantity, field + ".pdf"))
                plt.close()
        for field in ["zone", "ztwo", "whole"]:
            hmap_data = np.genfromtxt(cwd.joinpath("average", "thickness", field + ".dat"), delimiter=",")
            fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, "thickness", polar)
            plt.savefig(cwd.joinpath("figures", "thickness", field + ".pdf"))
            plt.close()
        for field in ["zone", "ztwo"]:
            hmap_data = np.genfromtxt(cwd.joinpath("average", "density", species, field + ".dat"), delimiter=",")
            fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, "density", polar)
            plt.savefig(cwd.joinpath("figures", "density", species, field + ".pdf"))
            plt.close()
            for tail in range(system_dict['ntails'][species]):
                hmap_data = np.genfromtxt(cwd.joinpath("average", "order", species, "tail" + str(tail), field + ".dat"), delimiter=",")
                fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, "order", polar)
                plt.savefig(cwd.joinpath("figures", "order", species, "tail" + str(tail), field + ".pdf"))
                plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Produce plots based on output from nougat.tcl")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl in polar coordinates")
    # parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")
    args = parser.parse_args()

    run_nougat(args.polar, False)

    print("Thank you for using nougat!")
