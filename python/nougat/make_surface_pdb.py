#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 09:28:42 2026.

@author: js2746
"""
from nougat.utils import compute_bin_centers


def make_pdb(filename, list_of_surfaces, list_of_names, bin_info, box_z=' 200.000'):
    """
    Make a pdb file from average height surfaces.

    Parameters
    ----------
    filename : str or path
        The name of the pdb file you wish to create.
    list_of_surfaces : list of numpy ndarrays
        List of 2D arrays containing the average height surface(s) you wish to
        represent in your pdb.
    list_of_names : list of str
        The names of each surface. Limit 4 characters per name.
    bin_info : named tuple
        The lattice information (bin widths, number, coordinate system, etc).
    box_z : str, OPTIONAL
        The height in z of the box. Default is 200 A.

    Returns
    -------
    None.

    """
    index = 0
    box_x = format_coordinate(bin_info.N1 * bin_info.d1)
    box_y = format_coordinate(bin_info.N2 * bin_info.d2)
    with open(filename, "w") as pdb:
        print(f'CRYST1 {box_x} {box_y} {box_z}  90.00  90.00  90.00 P 1           1', file=pdb)
        for surface, name in zip(list_of_surfaces, list_of_names):
            index = print_surface_to_pdb(surface, bin_info, pdb, index, name)
        print("END", file=pdb)


def print_surface_to_pdb(data, bin_info, f, index_num, field_name):
    """
    Add HETATOM records to a .pdb file for each bin in a nougat height Field.

    Parameters
    ----------
    data : numpy ndarray
        2D array containing heights from a nougat Field (e.g. z1, z2, etc).
    bin_info : named tuple
        The bin widths and numbers.
    f : file
        An open .pdb file you wish to add these lines to.
    index_num : int
        The index number for the first row you wish to print.
    field_name : str
        The name to be recorded in segname.

    Returns
    -------
    index_num : int
        The index number that should be provided to the next call to this fx.

    """
    resid_num = 1
    x_centers, y_centers = compute_bin_centers(bin_info)
    for d1bin in range(bin_info.N1):
        for d2bin in range(bin_info.N2):
            if str(data[d1bin][d2bin]) != "nan":
                x = x_centers[d1bin][d2bin]
                y = y_centers[d1bin][d2bin]
                print(
                    'HETATM'
                    f'{pad_str_with_spaces(index_num, 5)} '         # index
                    'SURF '                                         # name
                    f'{pad_str_with_spaces(field_name, 4, False)}'  # resname
                    'S'                                             # chain
                    f'{pad_str_with_spaces(resid_num, 5)}    '      # resid
                    f'{format_coordinate(x)}'                       # x
                    f'{format_coordinate(y)}'                       # y
                    f'{format_coordinate(data[d1bin][d2bin])}'      # z
                    f'{pad_str_with_spaces(d1bin, 3)}.00'           # occupancy
                    f'{pad_str_with_spaces(d2bin, 3)}.00'           # beta
                    f'      {field_name[:4]}  ',                    # segname
                    file=f,
                )
                index_num += 1
                resid_num += 1
    return index_num


def pad_str_with_spaces(inp, desired_len, left_pad=True):
    """
    Format a string to be of length desired_len by cutting or adding spaces.

    Parameters
    ----------
    inp : str, int, float
        The input you wish to turn into a string of a certain length.
    desired_len : int
        How many characters will your output string have?
    left_pad : bool, optional
        If True, add spaces to the left side of string. If False, add spaces to
        the right side of string. The default is True.

    Returns
    -------
    output_string : str
        The properly formatted output string.

    """
    input_string = str(inp)[:desired_len]
    padding = (desired_len - len(input_string)) * ' '
    if left_pad:
        output_string = padding + input_string
    else:
        output_string = input_string + padding
    return output_string


def format_coordinate(value):
    """
    Round an x/y coordinate and/or pad it with blank spaces.

    Creates the correct number of chars to fit in a pdb. Needs to handle the
    LHS and RHS of the decimal separately.

    Parameters
    ----------
    value : float
        A number.

    Returns
    -------
    float
        The same number, rounded and with blank spaces added to make it fit in
        a pdb file coordinate column.

    """
    rounded = round(value, 3)
    leftside, rightside = str(rounded).split('.')
    leftside = pad_str_with_spaces(leftside, 4)
    rightside = pad_str_with_spaces(rightside, 3, left_pad=False)
    return leftside + '.' + rightside
