#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 09:28:42 2026.

@author: js2746
"""
import numpy as np


def print_surface_to_pdb(data, d1, d2, f, index_num, field_name, coordsys):
    """
    Add HETATOM records to a .pdb file for each bin in a nougat height Field.

    Parameters
    ----------
    data : numpy ndarray
        2D array containing heights from a nougat Field (e.g. z1, z2, etc).
    d1 : float
        The bin width in x or r.
    d2 : float
        The bin width in y or theta.
    f : file
        An open .pdb file you wish to add these lines to.
    index_num : int
        The index number for the first row you wish to print.
    field_name : str
        The name to be recorded in segname.
    coordsys : str
        The coordinate system nougat used for your analysis. Must be "polar" or
        "cartesian".

    Returns
    -------
    index_num : int
        The index number that should be provided to the next call to this fx.

    """
    resid_num = 1
    chain = 'S'
    name = 'SURF '
    resname = pad_str_with_spaces(field_name, 4, left_pad=False)
    segname = '      ' + field_name[:4] + '  '

    row, col = data.shape
    for d1bin in range(row):
        for d2bin in range(col):
            if str(data[d1bin][d2bin]) != "nan":
                # calculate bin center in x and y
                if coordsys == "polar":
                    x = (d1 * d1bin + .5 * d1) * (np.cos(d2bin * d2 + 0.5 * d2))
                    y = (d1 * d1bin + .5 * d1) * (np.sin(d2bin * d2 + 0.5 * d2))
                else:
                    L1 = d1 * row
                    L2 = d2 * col
                    x = (d1 * d1bin + .5 * d1) - L1 / 2
                    y = (d2 * d2bin + .5 * d2) - L2 / 2

                # compile mis en place
                index = pad_str_with_spaces(index_num, 5) + ' '
                resid = pad_str_with_spaces(resid_num, 5) + '    '
                x = format_coordinate(x)
                y = format_coordinate(y)
                z = format_coordinate(data[d1bin][d2bin])
                occupancy = pad_str_with_spaces(d1bin, 3) + '.00'
                beta = pad_str_with_spaces(d2bin, 3) + '.00'

                # print row of .pdb file
                print('HETATM' + index + name + resname + chain + resid + x + y + z + occupancy + beta + segname, file=f)

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
