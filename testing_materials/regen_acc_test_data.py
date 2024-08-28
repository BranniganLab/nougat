#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 12:16:04 2024

@author: js2746
"""
import shutil
from nougat import run_nougat


def overwrite_ref_data_directories(path):
    """
    Create the preliminary directory hierarchy for nougat.py reference outputs.

    Parameters
    ----------
    path  :  Path object
        the path to the desired outfiles directory.

    Returns
    -------
    None.

    """
    quantities = ["height", "curvature", "thickness"]
    for filetype in ["trajectory", "average", "figures"]:
        for quantity in quantities:
            if quantity == "curvature":
                for curv in ["mean", "gaussian", "normal_vectors"]:
                    shutil.rmtree(path, ignore_errors=True)
                    dirname = path.joinpath(filetype, quantity, curv)
                    dirname.mkdir(parents=True, exist_ok=True)
            else:
                shutil.rmtree(path, ignore_errors=True)
                dirname = path.joinpath(filetype, quantity)
                dirname.mkdir(parents=True, exist_ok=True)


def regenerate_ref_data(path, polar):
    """
    Run nougat on the test system and print the correct fields to file with \
    naming convention that matches acceptance testing.

    Parameters
    ----------
    path : Path or str
        The path to the directory above the tcl_outputs directory for the test\
        system.
    polar : bool
        If True, use polar coordinates. If False, use Cartesian.

    Returns
    -------
    None.

    """
    m = run_nougat(path, polar)

    height = m.children['z']
    height_zero = m.children['z_zero']

    thickness = m.children['t']

    mean_curv = m.children['H']
    mean_zero = m.children['h_zero']

    gauss_curv = m.children['K']
    gauss_plus = m.children['k_plus']
    gauss_zero = m.children['k_zero']

    height.outer.print_to_file(path.joinpath('zone'))