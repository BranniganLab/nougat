#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:12:59 2024

@author: js2746
"""

from pathlib import Path
import argparse
import sys
import os
import matplotlib.pyplot as plt
from nougat import Membrane
from nougat.curvature import calculate_curvature


def make_todo_list(quantities):
    """
    Make the todo_list for the Membrane object.

    Parameters
    ----------
    quantities : str
        A string specifying which quantities should be analyzed by nougat.py.

    Raises
    ------
    ValueError
        The user must specify a valid nougat quantity.

    Returns
    -------
    todo_list  :  list
        A list of nougat Field_sets to be computed.

    """
    if quantities is None:
        todo_list = ["height", "curvature", "thickness", "order"]
    else:
        todo_list = []
        for letter in quantities:
            if letter == "h":
                todo_list.append("height")
            elif letter == "c":
                todo_list.append("curvature")
            elif letter == "t":
                todo_list.append("thickness")
            elif letter == "o":
                todo_list.append("order")
            elif letter == "n":
                return NotImplemented
                # todo_list.append("tilt")
            elif letter == "d":
                return NotImplemented
                # todo_list.append("density")
            else:
                raise ValueError("Must specify a valid nougat quantity")
    return todo_list


def run_nougat(path, polar, quantities):
    """
    Run nougat's averaging and image processing routines.

    Parameters
    ----------
    path : Path or str
        The path to your nougat results.
    polar : boolean
        True for cylindrical coordinate system, False for Cartesian.
    quantities : str
        A string that specifies which quantities to carry out analysis on. If\
        None, assume all quantities should be analyzed.

    Returns
    -------
    Membrane object containing the fields specified in $quantities.

    List of all default Fields and Field_sets
    ----------
    height  :  Field_set
        The height (in z) of the outer and inner leaflets, as well as the\
        symmetric/antisymmetric variables "z plus" and "z minus". Z plus is\
        the bilayer midplane. Z minus is the bilayer thickness. Heights may be\
        relative to some reference point on a protein or could be the absolute\
        height of the membrane in VMD.
    zzero  :  Field
        The height (in z) of the interface between the outer and inner\
        leaflet. "Z zero" is commonly thought to be equal to the bilayer\
        midplane but this is not always the case, especially around inclusions.
    thickness  :  Field_set
        The thickness of the outer and inner leaflets, as well as the\
        symmetric and anti-symmetric variables "t plus" and "t minus." T plus\
        is the average leaflet thickness and t minus is the leaflet thickness\
        asymmetry.
    epsilon  :  Field
        An alternative measurement of leaflet thickness asymmetry, defined in\
        [Watson & Brown, PRL, 2012].
    mean_curvature  :  Field_set
        The mean curvature (H) of the membrane outer and inner leaflets, as\
        well as the symmetric/anti-symmetric variables "H plus" and "H minus".
    gaussian_curvature  :  Field_set
        The Gaussian curvature (K) of the membrane outer and inner leaflets,\
        as well as the symmetric/anti-symmetric variables "K plus" and "K\
        minus".
    """
    if isinstance(path, str):
        path = Path(path)
    elif not isinstance(path, Path):
        raise Exception("path must be a Path object or a string.")

    todo_list = make_todo_list(quantities)

    m = Membrane(polar)

    if "height" in todo_list:
        zone = m.create_Field(path, "z_one", "height", "zone")
        ztwo = m.create_Field(path, "z_two", "height", "ztwo")
        zzero = m.create_Field(path, "z_zero", "height", "zzero")
        height = m.create_Field_set(zone, ztwo, "z")
    if "thickness" in todo_list:
        tone = m.create_Field(zone - zzero, "t_one")
        ttwo = m.create_Field(zzero - ztwo, "t_two")
        thickness = m.create_Field_set(tone, ttwo, "t")
    if "curvature" in todo_list:
        hone, kone, _ = calculate_curvature(zone, m.polar, m.grid_dims)
        htwo, ktwo, _ = calculate_curvature(ztwo, m.polar, m.grid_dims)
        hzero, kzero, _ = calculate_curvature(zzero, m.polar, m.grid_dims)

        hone = m.create_Field(hone, "h_one")
        kone = m.create_Field(kone, "k_one")
        htwo = m.create_Field(htwo, "h_two")
        ktwo = m.create_Field(ktwo, "k_two")

        mean_curv = m.create_Field_set(hone, htwo, "H")
        gauss_curv = m.create_Field_set(kone, ktwo, "K")

        hzero = m.create_Field(hzero, "h_zero")
        kzero = m.create_Field(kzero, "k_zero")

        # gaussian curvature of height.plus is not the same as gauss_curv.plus!
        _, kplus, _ = calculate_curvature(height.plus, m.polar, m.grid_dims)
        kplus = m.create_Field(kplus, "k_plus")

        # gaussian curvature of height.minus is not the same as gauss_curv.minus!
        _, kminus, _ = calculate_curvature(height.minus, m.polar, m.grid_dims)
        kminus = m.create_Field(kminus, "k_minus")

    return m


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze output from nougat.tcl")
    parser.add_argument("path", default=".", help="the path to your nougat outputs folder")
    parser.add_argument("-p", "--polar", action="store_true", help="add this flag if you ran nougat.tcl with polar coordinates")
    parser.add_argument("-q", "--quantities", help="Specify the quantities you want to calculate: height=h, thickness=t, curvature=c, order=o")
    parser.add_argument("-d", "--dump", action="store_true", help="Print all fields to file")
    parser.add_argument("-i", "--inclusion", action="store_true", help="add this flag if you ran nougat.tcl with Protein_Position turned on")

    args = parser.parse_args()
    path = Path(args.path)

    m = run_nougat(path, args.polar, args.quantities)

    if args.inclusion:
        m.add_protein_helices(path)

    print(m.helix_locations['zone'])

    fig, ax = m.plot2d(getattr(m.children['z'], 'outer'), 15, -15, helix_surface='zone')

    if args.dump:
        m.dump(path)

    print("Thank you for using nougat!")
