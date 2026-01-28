#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  5 09:28:42 2026.

@author: js2746
"""
from scipy.spatial import Delaunay  # pylint: disable-msg=E0611
import numpy as np
from nougat.utils import compute_bin_centers


def make_pdb(filename, list_of_surfaces, list_of_names, bin_info, box_dims=(200, 200, 200)):
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
    box_dims : tuple, OPTIONAL
        The dimensions of the box in Angstroms (x, y, z). Default is 200x200x200.

    Returns
    -------
    None.

    """
    index = 1
    with open(filename, "w", encoding='utf-8') as pdb:
        print(
            f"CRYST1  {box_dims[0]}.000  "
            f"{box_dims[1]}.000  "
            f"{box_dims[2]}.000  "
            "90.00  90.00  90.00 P 1           1",
            file=pdb,
        )
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
                    f'{pad_str_with_spaces(index_num, 5)} '             # index
                    'SURF '                                              # name
                    f'{pad_str_with_spaces(field_name, 3, False)}'    # resname
                    ' S'                                                # chain
                    f'{pad_str_with_spaces(resid_num, 4)}    '          # resid
                    f'{format_coordinate_for_pdb(x)}'                       # x
                    f'{format_coordinate_for_pdb(y)}'                       # y
                    f'{format_coordinate_for_pdb(data[d1bin][d2bin])}'      # z
                    f'{pad_str_with_spaces(d1bin, 3)}.00'           # occupancy
                    f'{pad_str_with_spaces(d2bin, 3)}.00'                # beta
                    f'      {field_name[:4]} C',          # segname and element
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


def format_coordinate_for_pdb(value):
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


def make_triangle_coordinates_file(xy, z, path):
    """
    Save a file that has triangle coordinate points on each row.

    Every three rows constitutes one triangle. Use nougat's drawTriangles proc
    to load into VMD.

    Parameters
    ----------
    xy : 2D numpy ndarray
        An array with two columns and as many rows as there are triangle points.
        Column 0 contains the x-coordinate and column 1 contains the y-coordinate
        for each point. This is the format required by scipy.spatial.Delaunay.
    z : list
        List of z-coordinates; same length as number of rows in xy.
    path : pathlib Path or str
        Path (including name and suffix) to file that will be created.

    Returns
    -------
    None.

    """
    if len(z) != xy.shape[0]:
        raise IndexError("xy must have same number of entries as z")
    triangles = Delaunay(xy)
    with open(path, 'w', encoding='utf-8') as f:
        for simplex in triangles.simplices:
            for index in simplex:
                print(xy[index][0], xy[index][1], z[index], file=f)


def format_triangle_points_and_values(surface_values, x_coords, y_coords):
    """
    Format xy coordinates and z coordinates to be useable by Delaunay module.

    Parameters
    ----------
    surface_values : 2D numpy ndarray
        A 2D array of values (e.g. average height over time) that will form
        the z component of your triangles.
    x_coords : 2D numpy ndarray
        The x-coordinates for every bin in the lattice.
    y_coords : 2D numpy ndarray
        The y-coordinates for every bin in the lattice.

    Returns
    -------
    2D numpy ndarray
        An array with two columns and as many rows as there are triangle points.
        Column 0 contains the x-coordinate and column 1 contains the y-coordinate
        for each point. This is the format required by scipy.spatial.Delaunay.
    list
        List of z-coordinates; same length as number of rows in xy ndarray.

    """
    if x_coords.shape != y_coords.shape:
        raise IndexError("x_coords and y_coords must be same shape.")
    num_rows, num_cols = x_coords.shape
    points_list = []
    values_list = []
    for row_i in range(num_rows):
        for col_j in range(num_cols):
            if not np.isnan(surface_values[row_i, col_j]):
                point = [x_coords[row_i, col_j], y_coords[row_i, col_j]]
                points_list.append(point)
                values_list.append(surface_values[row_i, col_j])
    return np.array(points_list), values_list


def save_surface_triangle_coordinates(path, surface, bin_info):
    """
    Save triangle coordinates to file to be read-in to molvis software (e.g. VMD).

    Parameters
    ----------
    path : pathlib Path
        The full path (including name) of the file you wish to save.
    surface : 2D numpy ndarray
        A 2D array of values (e.g. average height over time) that will form
        the z component of your triangles.
    bin_info : namedtuple
        Contains information about number of bins, step size, and coordinate
        system.

    Returns
    -------
    None.

    """
    x_centers, y_centers = compute_bin_centers(bin_info)
    points, values = format_triangle_points_and_values(surface, x_centers, y_centers)
    make_triangle_coordinates_file(points, values, path)
