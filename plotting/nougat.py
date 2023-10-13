"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import argparse
import os
import numpy as np
from height import *
from thickness import *
from curvature import *
from density import *
from tilt import *
from order import *
from utils import *


def run_nougat(sys_name, polar, inclusion_drawn, config_dict):
    """
    Run nougat's averaging and image processing routines.

    Parameters
    ----------
    sys_name : string
        The name assigned to the nougat.tcl run
    polar : boolean
        True for polar coordinate system; False for cartesian
    inclusion_drawn : boolean
        This feature currently isn't implemented, but would be how you acccount
        include an inclusion in your graphics
    config_dict : dict
        Dictionary containing useful information about the system

    Returns
    -------
    None.

    """
    cwd = os.getcwd()

    for filetype in ["npy", "dat", "pdf"]:
        dirname = os.path.join(cwd, filetype)
        try:
            os.mkdir(dirname)
        except OSError:
            continue

    if inclusion_drawn is True:
        inclusion = add_inclusion(name, field_list)  # this proc doesn't exist yet!
    else:
        inclusion = False

    if polar is True:
        coordsys = 'polar'
    elif polar is False:
        coordsys = 'cart'

    field_list = ["zone", "ztwo", "zzero"]

    # figure out all the important info about the system you'll need
    system_dict = read_log(sys_name, coordsys)
    save_areas(system_dict['bin_info']['N1'], system_dict['bin_info']['d1'], system_dict['bin_info']['N2'], system_dict['bin_info']['d2'], 0, coordsys, sys_name)

    # prep heatmap plot dimensions
    hmap_dims = bin_prep(system_dict['bin_info'], coordsys)

    # analyze height
    system_dict = analyze_height(sys_name, system_dict, coordsys, inclusion, polar, field_list)

    for bead in system_dict['headnames'].values():
        calculate_thickness(sys_name, bead, coordsys, inclusion, polar, hmap_dims, config_dict)
        calculate_curvature(sys_name, bead, coordsys, inclusion, polar, hmap_dims, field_list, config_dict, system_dict)

    calculate_density(sys_name, system_dict, coordsys, inclusion, polar, hmap_dims, config_dict)
    calculate_order(sys_name, system_dict, coordsys, inclusion, polar, hmap_dims, config_dict)
    # calculate_tilt(sys_name, system_dict, coordsys, inclusion, polar, hmap_dims, config_dict)

    calc_elastic_terms(sys_name, ".", coordsys, config_dict, system_dict['bin_info'])

    for species in system_dict['headnames'].keys():
        bead = system_dict['headnames'][species]
        for field in ["zone", "ztwo", "zzero", "zplus"]:
            for quantity in ['avgheight', 'avgKcurvature', 'avgcurvature']:
                hmap_data = np.genfromtxt("./dat/" + sys_name + "." + field + "." + bead + "." + coordsys + "." + quantity + ".dat")
                fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, quantity, bead, coordsys)
                save_figure(bead, sys_name, field, coordsys, quantity)
        for field in ["zone", "ztwo", "whole"]:
            hmap_data = np.genfromtxt("./dat/" + sys_name + "." + field + "." + bead + "." + coordsys + ".avgthickness.dat")
            fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, "avgthickness", bead, coordsys)
            save_figure(bead, sys_name, field, coordsys, "avg_thickness")
        for field in ["zone", "ztwo"]:
            hmap_data = np.genfromtxt("./dat/" + sys_name + "." + species + "." + field + "." + coordsys + ".avgdensity.dat")
            fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, "avgdensity", bead, coordsys)
            save_figure(bead, sys_name, field, coordsys, "avgdensity")
            for tail in range(system_dict['ntails'][species]):
                hmap_data = np.genfromtxt("./dat/" + sys_name + "." + species + ".tail" + str(tail) + "." + field + "." + coordsys + ".avgOrder.dat")
                fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, "avgorder", bead, coordsys)
                save_figure(bead, sys_name, field, coordsys, "avgorder")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Produce plots based on output from nougat.tcl")
    parser.add_argument("sys_name", help="what system did you name this?")
    parser.add_argument("config", help="what config file should nougat use?")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl in polar coordinates")
    # parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")
    args = parser.parse_args()

    config_dict = read_config(args.config)

    run_nougat(args.sys_name, args.polar, False, config_dict)

    print("Thank you for using nougat!")
