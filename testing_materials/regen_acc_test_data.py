#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 12:16:04 2024

@author: js2746
"""
from nougat import run_nougat
from pathlib import Path


def create_ref_data_directories(path):
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
                    dirname = path.joinpath(filetype, quantity, curv)
                    dirname.mkdir(parents=True, exist_ok=True)
            else:
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

    for filetype in ['trajectory', 'average']:
        height.outer.save_to_file(path.joinpath(filetype, 'height'), filetype, name="zone")
        height.inner.save_to_file(path.joinpath(filetype, 'height'), filetype, name="ztwo")
        height.plus.save_to_file(path.joinpath(filetype, 'height'), filetype, name="zplus")
        height_zero.save_to_file(path.joinpath(filetype, 'height'), filetype, name="zzero")

        thickness.outer.save_to_file(path.joinpath(filetype, 'thickness'), filetype, name="zone")
        thickness.inner.save_to_file(path.joinpath(filetype, 'thickness'), filetype, name="ztwo")

        mean_curv.outer.save_to_file(path.joinpath(filetype, 'curvature', 'mean'), filetype, name="zone")
        mean_curv.inner.save_to_file(path.joinpath(filetype, 'curvature', 'mean'), filetype, name="ztwo")
        mean_curv.plus.save_to_file(path.joinpath(filetype, 'curvature', 'mean'), filetype, name="zplus")
        mean_zero.save_to_file(path.joinpath(filetype, 'curvature', 'mean'), filetype, name="zzero")

        gauss_curv.outer.save_to_file(path.joinpath(filetype, 'curvature', 'gaussian'), filetype, name="zone")
        gauss_curv.inner.save_to_file(path.joinpath(filetype, 'curvature', 'gaussian'), filetype, name="ztwo")
        gauss_plus.save_to_file(path.joinpath(filetype, 'curvature', 'mean'), filetype, name="zplus")
        gauss_zero.save_to_file(path.joinpath(filetype, 'curvature', 'mean'), filetype, name="zzero")


if __name__ == '__main__':
    