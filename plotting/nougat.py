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

    print(system_dict)

    for species in system_dict['species']:
        for field in ["zone", "ztwo", "zzero", "zplus"]:
            for quantity in ['avgheight', 'avgKcurvature', 'avgcurvature']:
                print(sys_name, field, coordsys, quantity)
                hmap_data = np.genfromtxt("./dat/" + sys_name + "." + field + "." + coordsys + "." + quantity + ".dat", delimiter=",")
                fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, quantity, coordsys)
                save_figure(sys_name, field, coordsys, quantity)
        for field in ["zone", "ztwo", "whole"]:
            hmap_data = np.genfromtxt("./dat/" + sys_name + "." + field + "." + coordsys + ".avgthickness.dat", delimiter=",")
            fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, "avgthickness", coordsys)
            save_figure(sys_name, field, coordsys, "avg_thickness")
        for field in ["zone", "ztwo"]:
            hmap_data = np.genfromtxt(cwd.joinpath("average", "density", species, field + ".dat"), delimiter=",")
            fig, ax = plot_maker(hmap_dims, hmap_data, field, inclusion, "avgdensity", polar)
            plt.save(cwd.joinpath("figures", "density", species, field + ".pdf"))
            for tail in range(system_dict['ntails'][species]):
                hmap_data = np.genfromtxt("./dat/" + sys_name + "." + species + ".tail" + str(tail) + "." + field + "." + coordsys + ".avgOrder.dat", delimiter=",")
                fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, "avgorder", coordsys)
                plt.savefig('pdf/' + sys_name + "_" + field + "_" + coordsys + "tail" + str(tail) + "_avgorder.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Produce plots based on output from nougat.tcl")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl in polar coordinates")
    # parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")
    args = parser.parse_args()

    run_nougat(args.polar, False)

    print("Thank you for using nougat!")
