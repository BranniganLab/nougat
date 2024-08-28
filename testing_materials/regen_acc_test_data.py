#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 12:16:04 2024.

@author: js2746
"""
from pathlib import Path
import argparse
from shutil import rmtree
import os
import sys

sys.path.append(os.path.abspath('../plotting/'))
from nougat import run_nougat


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
        dirname = path.joinpath(filetype, 'height')
        dirname.mkdir(parents=True, exist_ok=False)
        height.outer.save_to_file(dirname, filetype, name="zone")
        height.inner.save_to_file(dirname, filetype, name="ztwo")
        height.plus.save_to_file(dirname, filetype, name="zplus")
        height_zero.save_to_file(dirname, filetype, name="zzero")

        dirname = path.joinpath(filetype, 'thickness')
        dirname.mkdir(parents=True, exist_ok=False)
        thickness.outer.save_to_file(dirname, filetype, name="zone")
        thickness.inner.save_to_file(dirname, filetype, name="ztwo")

        dirname = path.joinpath(filetype, 'curvature', 'mean')
        dirname.mkdir(parents=True, exist_ok=False)
        mean_curv.outer.save_to_file(dirname, filetype, name="zone")
        mean_curv.inner.save_to_file(dirname, filetype, name="ztwo")
        mean_curv.plus.save_to_file(dirname, filetype, name="zplus")
        mean_zero.save_to_file(dirname, filetype, name="zzero")

        dirname = path.joinpath(filetype, 'curvature', 'gaussian')
        dirname.mkdir(parents=True, exist_ok=False)
        gauss_curv.outer.save_to_file(dirname, filetype, name="zone")
        gauss_curv.inner.save_to_file(dirname, filetype, name="ztwo")
        gauss_plus.save_to_file(dirname, filetype, name="zplus")
        gauss_zero.save_to_file(dirname, filetype, name="zzero")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Regenerate reference data files for acceptance testing.")
    parser.add_argument("path", default=".", help="the path to the directory above your tcl_outputs folder")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl with polar coordinates")
    parser.add_argument("-o", "--overwrite", action="store_true", help="add this flag if you want to overwrite existing python outputs")

    args = parser.parse_args()
    path = Path(args.path)

    if args.overwrite:
        rmtree(path.joinpath('trajectory'))
        rmtree(path.joinpath('average'))
    regenerate_ref_data(path, args.polar)
