"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import argparse
from pathlib import Path

import numpy as np
from height import *
from thickness import *
from curvature import *
from density import *
from tilt import *
from order import *
from utils import *


def run_nougat(polar, inclusion_drawn, config_dict):
    """
    Run nougat's averaging and image processing routines.

    Parameters
    ----------
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

    # make necessary folders
    cwd = Path.cwd()
    for filetype in ["trajectory", "average"]:
        for quantity in ["height", "density", "curvature", "thickness", "order", "tilt", "misc"]:
            if quantity == "curvature":
                for curv in ["mean", "gaussian", "normal_vectors"]:
                    dirname = cwd.joinpath(filetype, quantity, curv)
                    dirname.mkdir(parents=True, exist_ok=True)
            else:
                dirname = cwd.joinpath(filetype, quantity)
                dirname.mkdir(parents=True, exist_ok=True)

    # define inclusion if present
    if inclusion_drawn is True:
        inclusion = add_inclusion(name, field_list)  # this proc doesn't exist yet!
    else:
        inclusion = False

    # TO DO: eliminate coordsys
    if polar is True:
        coordsys = 'polar'
    elif polar is False:
        coordsys = 'cart'

    # figure out all the important info about the system you'll need
    system_dict = read_log()
    sys_name = system_dict["sysname"]

    # prep heatmap plot dimensions
    hmap_dims = bin_prep(system_dict['bin_info'], polar)

    # analyze height
    system_dict = analyze_height(sys_name, system_dict, coordsys, inclusion, cwd)

    calculate_thickness(sys_name, coordsys, inclusion, hmap_dims, config_dict, cwd)
    calculate_curvature(sys_name, coordsys, system_dict, cwd)

    save_areas(system_dict["bin_info"], 0, polar, sys_name, cwd)

    calculate_density(sys_name, system_dict, coordsys, inclusion, hmap_dims, config_dict, cwd)
    calculate_order(sys_name, system_dict, coordsys, inclusion, hmap_dims, config_dict, cwd)
    # calculate_tilt(sys_name, system_dict, coordsys, inclusion, hmap_dims, config_dict, cwd)

    calc_elastic_terms(sys_name, ".", coordsys, config_dict, system_dict['bin_info'])


"""
    for species in system_dict['headnames'].keys():
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
            hmap_data = np.genfromtxt("./dat/" + sys_name + "." + species + "." + field + "." + coordsys + ".avgdensity.dat", delimiter=",")
            fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, "avgdensity", coordsys)
            save_figure(sys_name, field, coordsys, "avgdensity")
            for tail in range(system_dict['ntails'][species]):
                hmap_data = np.genfromtxt("./dat/" + sys_name + "." + species + ".tail" + str(tail) + "." + field + "." + coordsys + ".avgOrder.dat", delimiter=",")
                fig, ax = plot_maker(hmap_dims, hmap_data, sys_name, field, config_dict, inclusion, "avgorder", coordsys)
                plt.savefig('pdf/' + sys_name + "_" + field + "_" + coordsys + "tail" + str(tail) + "_avgorder.pdf")
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Produce plots based on output from nougat.tcl")
    parser.add_argument("config", help="what config file should nougat use?")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl in polar coordinates")
    # parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")
    args = parser.parse_args()

    config_dict = read_config(args.config)

    run_nougat(args.polar, False, config_dict)

    print("Thank you for using nougat!")
