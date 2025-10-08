"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings


def strip_blank_lines(file):
    """
    Remove any lines that are empty when reading a file.

    Parameters
    ----------
    file : list
        A list of all lines in a file.

    Yields
    ------
    stripped_line : string
        Each line that has text in it.

    """
    for line in file:
        stripped_line = line.strip()
        if stripped_line:
            yield stripped_line


def read_config(path):
    """
    Read the nougat.py config file and save everything to a dict.

    Parameters
    ----------
    path : string
        Path to the nougat.py config file.

    Returns
    -------
    config_dict : dict
        Dict containing all config entries in key/val pairs.

    """
    config_dict = {}
    with open(path, "r+") as config_file:
        for line in strip_blank_lines(config_file):
            if line.startswith("#") is True:
                continue
            else:
                line = line.partition('#')[0]
                if line.rstrip():
                    key = line.partition('=')[0].strip()
                    value = line.partition('=')[2].strip()
                    config_dict.update({key: value})

    return config_dict


def gifformat(num, size):
    """
    Format number with proper amount of zeros in front.

    Parameters
    ----------
    num : float/int
        The number in need of formatting.
    size : int
        The number of spaces it needs to fill.

    Returns
    -------
    padded_val : string
        The number with the appropriate amount of zeros in front.
    """
    numzeros = size - len(str(num))
    padded_val = "0" * numzeros + str(num)
    return padded_val


def calc_avg_over_time(matrix_data):
    """
    Calculate mean over 3rd axis (time, in nougat terms) of an array. np.nanmean \
        is wrapped inside of a catch_warnings loop so that the user does not get \
        bombarded with a bunch of warnings every time a nan is encountered, \
        which is supposed to be the whole point of nanmean anyway.

    Parameters
    ----------
    matrix_data : numpy ndarray
        A 3d array.

    Returns
    -------
    avg : numpy ndarray
        A 2d array (matrix_data averaged over the 3rd dimension).

    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avg = np.nanmean(matrix_data, axis=0)
        return avg


def bin_prep(bin_info, polar):
    """
    Configure the arrays needed for plotting heatmaps.

    Parameters
    ----------
    bin_info : DICT
        Contains N1, N2, d1, and d2 information.
    polar : BOOL
        Whether or not to use polar coordinates.

    Returns
    -------
    list
        The two numpy ndarrays needed for plotting a heatmap.

    """
    dim1 = np.linspace(0, bin_info['N1'] * bin_info['d1'], bin_info['N1'] + 1)
    if polar:
        dim2 = np.linspace(0, 2 * np.pi, bin_info['N2'] + 1)
    else:
        dim2 = dim1
    dim1vals, dim2vals = np.meshgrid(dim1, dim2, indexing='ij')

    return [dim1vals, dim2vals]


def mostly_empty(data_array):
    """
    Replace bin values with np.nan if that bin has lipids in it less than 10% \
        of the trajectory frames.

    Parameters
    ----------
    data_array : numpy ndarray
        The input data in a 3D array with dimensions time, [x, r], [y, theta].

    Returns
    -------
    data_array : numpy ndarray
        A pruned array.

    """
    # if a bin only has lipids in it <10% of the time, it shouldn't be considered part of the membrane
    Nframes, N1_bins, N2_bins = np.shape(data_array)
    for row in range(N1_bins):
        for col in range(N2_bins):
            zerocount = np.count_nonzero(data_array[:, row, col])
            count = np.count_nonzero(np.isnan(data_array[:, row, col]))
            if (zerocount - count) / Nframes <= .1:
                data_array[:, row, col] = np.nan
    return data_array


def read_log(input_path):
    """
    Read log file output by nougat.tcl and save important info for later.

    Parameters
    ----------
    input_path  :  Path or str
        The Path object or string path to your nougat.log file.

    Returns
    -------
    system_dict : DICT
        A dictionary containing the list of lipid species, their respective \
            headnames, their respective density normalization factors, and \
            the bin sizes.

    """
    system_dict = {}

    # open log file
    with open(input_path, "r+") as log_file:
        lines = [line.rstrip('\n') for line in log_file]

        system_dict["sysname"] = lines[1]
        system_dict["coordsys"] = lines[2]

        species_list = []
        # get all lipid species names from species section
        species_start_line = lines.index("#SYSTEM CONTENTS")
        for species in lines[species_start_line + 1].split(' '):
            species_list.append(species)
        system_dict["species"] = species_list

        # get number of tails on each lipid
        tails_start_line = lines.index("#NUMBER OF TAILS") + 1
        system_dict["ntails"] = {}
        for line in range(len(system_dict["species"])):
            names_line = lines[tails_start_line].split(':')
            system_dict["ntails"][names_line[0]] = int(names_line[1])
            tails_start_line += 1

        # get density norm info from density section
        density_start_line = lines.index("#DENSITY NORMALIZATION") + 1
        system_dict["density_norm"] = {}
        for line in range(len(system_dict["species"])):
            names_line = lines[density_start_line].split(':')
            system_dict["density_norm"][names_line[0]] = float(names_line[1])
            density_start_line += 1

        # get bin size info from bin info section
        bin_start_line = lines.index("#BIN INFO") + 1
        N1, N2 = np.int64(lines[bin_start_line].split(' '))
        d1, d2 = np.float64(lines[bin_start_line + 1].split(' '))
        system_dict['bin_info'] = {"N1": N1, "N2": N2, "d1": d1, "d2": d2}

    return system_dict


def plot_maker(dims, data, polar, vmax, vmin, protein=False):
    """
    Create and save 2D heatmaps.

    Parameters
    ----------
    dims : list
        np.meshgrid output
    data : array
        the 2d array/matrix of values to be heatmapped
    polar : bool
        Whether or not to use polar coordinates
    vmax : float
        The maximum value for your colorbar
    vmin : float
        The minimum value for your colorbar
    protein : list or False
        if --inclusion turned on, list of helix coordinates; if no protein, False

    Returns
    -------
    matplotlib figure and axes objects.

    """
    dim1vals, dim2vals = dims

    fig = plt.figure()
    if polar:
        ax = plt.subplot(projection="polar")
    else:
        ax = plt.subplot()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    create_heatmap(polar, dim1vals, dim2vals, data, vmax, vmin, True)
    if protein:
        draw_protein(protein, polar)
    return fig, ax


def create_heatmap(polar, dim1vals, dim2vals, data, Vmax, Vmin, colorbar):
    """
    Create a 2d heatmap of your data.

    Parameters
    ----------
    polar : bool
        Whether or not to use polar coordinates
    dim1vals : list
        meshgrid output 1
    dim2vals : list
        meshgrid output 2
    data : numpy ndarray
        the 2d array of values to be heatmapped
    Vmax : float
        max value for colorbar
    Vmin : float
        min value for colorbar
    colorbar : bool
        draws colorbar legend if True

    Returns
    -------
    None.
    """
    if polar:
        if Vmax != "auto":
            c = plt.pcolormesh(dim2vals, dim1vals, data, cmap="RdBu_r", zorder=0, vmax=Vmax, vmin=Vmin)
        else:
            c = plt.pcolormesh(dim2vals, dim1vals, data, cmap="RdBu_r", zorder=0)
    else:
        if Vmax != "auto":
            c = plt.pcolormesh(dim1vals, dim2vals, data, cmap="RdBu_r", zorder=0, vmax=Vmax, vmin=Vmin)
        else:
            c = plt.pcolormesh(dim1vals, dim2vals, data, cmap="RdBu_r", zorder=0)

    if colorbar:
        plt.colorbar(c)

    plt.axis('off')


def draw_protein(protein, polar):
    """
    Draw protein alpha helix positions.

    Parameters
    ----------
    protein : list
        List of helix coordinates
    polar : bool
        Whether or not to use polar coordinates

    Returns
    -------
    None.

    """
    for i in range(0, 10, 2):
        protein[i + 1] = np.deg2rad(protein[i + 1])
        if polar is False:
            protein[i], protein[i + 1] = convert_to_cart(protein[i], protein[i + 1])
        plt.scatter(protein[i + 1], protein[i], c="black", linewidth=4, zorder=2)


def convert_to_cart(rval, thetaval):
    """
    Convert from polar to cartesian coordinates.

    Parameters
    ----------
    rval : float
        The r coordinate value.
    thetaval : float
        The theta coordinate value.

    Returns
    -------
    xval : float
        The x coordinate value.
    yval : float
        The y coordinate value.

    """
    xval = rval * np.cos(thetaval)
    yval = rval * np.sin(thetaval)
    return xval, yval


def coord_format(value):
    """
    Round an x/y coordinate and/or pad it with blank spaces so that it is the \
        correct number of chars to fit in a pdb.

    Parameters
    ----------
    value : float
        A number.

    Returns
    -------
    final_value : float
        The same number, rounded and with blank spaces added to make it fit in \
            a pdb file coordinate column.

    """
    rounded = round(value, 3)
    leftside, rightside = str(rounded).split('.')
    if len(rightside) < 3:
        rightside = rightside + (' ' * (3 - len(rightside)))
    if len(leftside) < 4:
        leftside = (' ' * (4 - len(leftside))) + leftside
    final_value = leftside + '.' + rightside
    return final_value


def bin_format(value):
    """
    Round a bin number and/or pad it with blank spaces so that it is the \
        correct number of chars to fit in a pdb.

    Parameters
    ----------
    value : int
        The number of the bin.

    Returns
    -------
    final_value : string
        The same number, with .00 appended to the end and the correct number \
            of blank spaces in front to fit in the pdb column.

    """
    strval = str(value)
    length = len(strval)
    final_value = (' ' * (3 - length)) + strval + '.00'
    return final_value


def measure_quant_in_empty_sys(path, system, polar, quantity):
    """
    Measure the average thickness of a membrane.

    Parameters
    ----------
    path : string
        path to the directory where your nougat outputs are
    system : string
        name of the system you gave nougat
    polar : bool
        Whether or not to use polar coordinates\
            if polar, will ignore small r bins (area too small)
    quantity : string
        The quantity you want to take the average of

    Returns
    -------
    avg : float
        the average quantity of the membrane

    """
    data = np.load(path + '/npy/' + system + '.' + quantity + '.npy')
    if polar:
        data = data[:, 4:, :]  # this could be smarter
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avg = np.nanmean(data)
    return avg
