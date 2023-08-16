"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import matplotlib.pyplot as plt
from matplotlib import animation
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


def find_first_val(in_list):
    """
    Find first non-nan value in a list.

    Parameters
    ----------
    in_list : list
        A list.

    Returns
    -------
    float
        The first non-nan value in the list. Returns nan if no such value exists.

    """
    for value in in_list:
        if np.isnan(value):
            continue
        else:
            return value
    return np.nan


def find_last_val(in_list):
    """
    Find last non-nan value in a list.

    Parameters
    ----------
    in_list : list
        A list.

    Returns
    -------
    float
        The last non-nan value in the list. Returns nan if no such value exists.

    """
    last_index = len(in_list) - 1
    while last_index >= 0:
        if not np.isnan(in_list[last_index]):
            return in_list[last_index]
        else:
            last_index -= 1
    return np.nan


def filename_generator(sys_name, lipid_name, field, beadname, coordsys, measure, dtype):
    if measure == "height" or measure == "curvature" or measure == "Kcurvature" or measure == "thickness":
        if dtype == "dat":
            if measure == "thickness":
                fullmeasure = "normthickness"
            else:
                fullmeasure = "avg" + measure
        elif dtype == "npy":
            if measure == "curvature":
                fullmeasure = "meancurvature"
            elif measure == "Kcurvature":
                fullmeasure = "gausscurvature"
            else:
                fullmeasure = measure
        filename = sys_name + "." + field + "." + beadname + "." + coordsys + "." + fullmeasure + "." + dtype
    elif measure == "density":
        if dtype == "dat":
            filename = sys_name + "." + lipid_name + "." + field + "." + coordsys + ".avg" + measure + "." + dtype
        elif dtype == "npy":
            filename = sys_name + "." + lipid_name + "." + field + "." + coordsys + "." + measure + "." + dtype
    elif measure == "tail1" or measure == "tail0":
        if dtype == "dat":
            filename = sys_name + "." + lipid_name + "." + measure + "." + field + "." + coordsys + ".avgOrder." + dtype
        elif dtype == "npy":
            filename = sys_name + "." + lipid_name + "." + measure + "." + field + "." + coordsys + ".order." + dtype
    if measure not in ["height", "curvature", "Kcurvature", "thickness", "density", "tail1", "tail0"]:
        raise RuntimeWarning('You used \"' + measure + '\" but that is not an allowed measurement')
    return filename


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
        avg = np.nanmean(matrix_data, axis=2)
        return avg


def bin_prep(sys_name, beadnames, coordsys, density):

    sample_data = np.genfromtxt('tcl_output/' + sys_name + '.zone.' + beadnames + '.' + coordsys + '.height.dat', missing_values='nan', filling_values=np.nan)

    N1_bins, d1, N2_bins, d2, Nframes, min_val = dimensions_analyzer(sample_data, coordsys)

    # prep plot dimensions
    dim1 = sample_data[0:N1_bins, 0]
    dim1 = np.append(dim1, sample_data[N1_bins - 1, 1])
    if coordsys == "polar":
        dim2 = np.linspace(0, 2 * np.pi, N2_bins + 1)
    elif coordsys == "cart":
        dim2 = np.linspace(0, N2_bins + 1, N2_bins + 1)
    dim1vals, dim2vals = np.meshgrid(dim1, dim2, indexing='ij')

    if density == "ON":
        # save an array that represents the area per bin for normalizing density later
        save_areas(N1_bins, d1, N2_bins, d2, min_val, coordsys, sys_name)

    return [N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals]


def save_areas(N1_bins, d1, N2_bins, d2, min_val, coordsys, sys_name):

    areas = np.ones([N1_bins, N2_bins])

    areas = areas * d1 * d2
    if coordsys == "polar":
        for row in range(N1_bins):
            dist_to_center = min_val + row * d1 + d1 / 2.0
            areas[row, :] = areas[row, :] * dist_to_center
    np.save('npy/' + sys_name + "." + coordsys + ".areas.npy", areas)


def mostly_empty(data_array):
    # if a bin only has lipids in it <10% of the time, it shouldn't be considered part of the membrane
    N1_bins, N2_bins, Nframes = np.shape(data_array)
    for row in range(N1_bins):
        for col in range(N2_bins):
            zerocount = np.count_nonzero(data_array[row, col, :])
            count = np.count_nonzero(np.isnan(data_array[row, col, :]))
            if (zerocount - count) / Nframes <= .1:
                data_array[row, col, :] = np.nan
    return data_array


def read_log(sys_name, coordsys):
    # this proc is not robust to multiple lipid species in the same system!!
    # specifically, species with differing headnames
    names_dict = {}
    names_dict['beads_list'] = []
    # open log file
    with open("tcl_output/" + sys_name + "." + coordsys + ".log", "r+") as log_file:
        lines = [line.rstrip('\n') for line in log_file]

        # get contents of line 1 and save as species_list
        names_dict['species_list'] = lines[1].split(' ')

        # get contents of line 2 and save as beads_list
        headnames = lines[2].split(' ')
        filename = headnames[0]
        for indx in range(1, len(headnames)):
            filename = filename + "." + headnames[indx]
        if filename not in names_dict['beads_list']:
            names_dict['beads_list'].append(filename)

        # get density norm factor
        start_line = lines.index("#DENSITY NORMALIZATION") + 1
        names_dict['density_norm'] = float(lines[start_line].split(":")[1])

    return names_dict


def plot_maker(dim1vals, dim2vals, data, name, field, Vmax, Vmin, protein, dataname, bead, coordsys, config_dict):
    """
    Create and save 2D heatmaps.

    Parameters
    ----------
    dim1vals : list
        meshgrid output 1
    dim2vals : list
        meshgrid output 2
    data : array
        the 2d array/matrix of values to be heatmapped
    name : string
        the system name you gave nougat.py
    field : string
        usually describes which membrane field (z1, z2, etc) to be heatmapped
    Vmax : float
        max value for colorbar
    Vmin : float
        min value for colorbar
    protein : list or False
        if --inclusion turned on, list of helix coordinates; if no protein, False
    dataname : string
        the type of measurement (thickness, height, curvature, etc)
    bead : string or False
        if bead specified, name of bead; else False
    coordsys : string
        "polar" or "cart"
    config_dict : dict
        Dict containing config info

    Returns
    -------
    None.

    """
    fig = plt.figure()
    if coordsys == "polar":
        ax = plt.subplot(projection="polar")
    else:
        ax = plt.subplot()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    # fig.set_size_inches(6,6)
    create_heatmap(coordsys, dim1vals, dim2vals, data, Vmax, Vmin, config_dict['colorbar'])
    if protein:
        draw_protein(protein, coordsys)
    save_figure(bead, name, field, coordsys, dataname)
    plt.clf()
    plt.close()


def create_heatmap(coordsys, dim1vals, dim2vals, data, Vmax, Vmin, colorbar):
    """
    Create a 2d heatmap of your data.

    Parameters
    ----------
    coordsys : string
        "cart" or "polar"
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
    if coordsys == "polar":
        if Vmax != "auto":
            c = plt.pcolormesh(dim2vals, dim1vals, data, cmap="RdBu_r", zorder=0, vmax=Vmax, vmin=Vmin)
        else:
            c = plt.pcolormesh(dim2vals, dim1vals, data, cmap="RdBu_r", zorder=0)
    elif coordsys == "cart":
        if Vmax != "auto":
            c = plt.pcolormesh(dim1vals, dim2vals, data, cmap="RdBu_r", zorder=0, vmax=Vmax, vmin=Vmin)
        else:
            c = plt.pcolormesh(dim1vals, dim2vals, data, cmap="RdBu_r", zorder=0)
    else:
        print("something's wrong with coordsys")

    if colorbar:
        plt.colorbar(c)

    plt.axis('off')


def draw_protein(protein, coordsys):
    """
    Draw protein alpha helix positions.

    Parameters
    ----------
    protein : list or False
        if --inclusion turned on, list of helix coordinates; if no protein, False
    coordsys : string
        "cart" or "polar"

    Returns
    -------
    None.

    """
    for i in range(0, 10, 2):
        protein[i + 1] = np.deg2rad(protein[i + 1])
        if coordsys == "cart":
            protein[i], protein[i + 1] = convert_to_cart(protein[i], protein[i + 1])
        plt.scatter(protein[i + 1], protein[i], c="black", linewidth=4, zorder=2)

    # This circle is a custom E protein thing and should be removed
    # circle1 = plt.Circle((0, 0), 28.116, transform=ax.transData._b, color='black', linestyle='dashed', linewidth=4, fill=False)
    # if field == "zone":
    #    ax.add_artist(circle1)


def save_figure(bead, name, field, coordsys, dataname):
    """
    Save the current figure.

    Parameters
    ----------
    bead : string or False
        if bead specified, name of bead; else False
    name : string
        the system name you gave nougat.py
    field : string
        usually describes which membrane field (z1, z2, etc) to be heatmapped
    coordsys : string
        "cart" or "polar"
    dataname : string
        the type of measurement (thickness, height, curvature, etc)

    Returns
    -------
    None.

    """
    if bead is False:
        plt.savefig('pdf/' + name + "_" + field + "_" + coordsys + "_" + dataname + ".pdf", dpi=700)
    else:
        plt.savefig('pdf/' + name + "_" + bead + "_" + field + "_" + coordsys + "_" + dataname + ".pdf", dpi=700)


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


def dimensions_analyzer(data, coordsys):
    # figure out how many radial or x bins there are
    counter = 1
    flag = True
    match_value = data[0, 0]
    while (flag is True):
        try:
            if data[counter, 0] == match_value:
                flag = False
            else:
                counter = counter + 1
        # what if there is only 1 frame? Will raise IndexError
        except IndexError:
            flag = False
    N1_bins = counter

    # figure out how many azimuthal or y bins there are
    N2_bins = len(data[0, :]) - 2

    # figure out how many frames there are in the traj
    Nframes = int(len(data[:, 0]) / N1_bins)

    # error check
    if (len(data[:, 0]) % N1_bins) != 0:
        raise Exception("There is something wrong with the Nframes calculation")

    if coordsys == "polar":
        d1 = data[0, 1] - data[0, 0]
        d2 = (np.pi * 2) / N2_bins
    elif coordsys == "cart":
        # compute average d1, assume d2 is the same
        d1list = []
        for row in range(Nframes):
            d1list.append(data[row * N1_bins, 1])
        d1 = np.mean(d1list)
        d2 = d1

    return N1_bins, d1, N2_bins, d2, Nframes, match_value


def calc_elastic_terms(system, path, coordsys, scale_dict):
    """
    Calculate all the additional terms that appear in a hamiltonian or are \
        generally of interest.

    Parameters
    ----------
    system : string
        the same name you gave nougat.tcl and nougat.py
    path : string
        should point to the folder housing your nougat outputs for the given \
            system
    coordsys : string
        "polar" or "cart"
    scale_dict : dict
        contains scale bounds from the nougat config file

    Returns
    -------
    None.

    """
    # load height and curvature data
    z_1 = np.load(path + '/npy/' + system + '.zone.C1A.C1B.' + coordsys + '.height.npy')
    z_2 = np.load(path + '/npy/' + system + '.ztwo.C1A.C1B.' + coordsys + '.height.npy')
    z_0 = np.load(path + '/npy/' + system + '.zzero.C1A.C1B.' + coordsys + '.height.npy')
    z_plus = np.load(path + '/npy/' + system + '.zplus.C1A.C1B.' + coordsys + '.height.npy')
    H_1 = np.load(path + '/npy/' + system + '.zone.C1A.C1B.' + coordsys + '.meancurvature.npy')
    H_2 = np.load(path + '/npy/' + system + '.ztwo.C1A.C1B.' + coordsys + '.meancurvature.npy')
    K_1 = np.load(path + '/npy/' + system + '.zone.C1A.C1B.' + coordsys + '.gausscurvature.npy')
    K_2 = np.load(path + '/npy/' + system + '.ztwo.C1A.C1B.' + coordsys + '.gausscurvature.npy')

    # measure terms of interest
    # removed z_minus terms until we have a better way of computing t0
    epsilon = z_plus - z_0
    epsilon2 = epsilon**2
    H_plus = (H_1 + H_2) / 2
    K_plus = (K_1 + K_2) / 2
    K_minus = (K_1 - K_2) / 2
    H_plus2 = H_plus**2
    H_minus = (H_1 - H_2) / 2
    H_minus2 = H_minus**2
    epsilon_H = epsilon * H_plus
    total_t = z_1 - z_2

    # save useful trajectories
    np.save(path + '/npy/' + system + '.epsilon.npy', epsilon)
    np.save(path + '/npy/' + system + '.H_plus.npy', H_plus)
    np.save(path + '/npy/' + system + '.epsilon2.npy', epsilon2)
    np.save(path + '/npy/' + system + '.H_plus2.npy', H_plus2)
    np.save(path + '/npy/' + system + '.total_t.npy', total_t)

    # calculate averages
    avg_epsilon = calc_avg_over_time(epsilon)
    avg_epsilon2 = calc_avg_over_time(epsilon2)
    avg_H_plus = calc_avg_over_time(H_plus)
    avg_H_plus2 = calc_avg_over_time(H_plus2)
    avg_H_minus = calc_avg_over_time(H_minus)
    avg_H_minus2 = calc_avg_over_time(H_minus2)
    avg_epsilon_H = calc_avg_over_time(epsilon_H)
    avg_total_t = calc_avg_over_time(total_t)
    avg_K_plus = calc_avg_over_time(K_plus)
    avg_K_minus = calc_avg_over_time(K_minus)

    # calculate correlations
    corr_eps_Hplus = calc_avg_over_time(epsilon * H_plus) - (avg_epsilon * avg_H_plus)
    corr_mag_eps_Hplus = calc_avg_over_time(np.sqrt(epsilon2) * np.sqrt(H_plus2)) - (np.sqrt(avg_epsilon2) * np.sqrt(avg_H_plus2))
    corr_eps_Kplus = calc_avg_over_time(epsilon * K_plus) - (avg_epsilon * avg_K_plus)

    # get proper plot dimensions
    dims = bin_prep(system, "C1A.C1B", coordsys, "OFF")
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    # measure average thickness
    avgt0 = measure_t0(path, system, coordsys)
    np.save(path + '/npy/' + system + '.avg_t0.npy', avgt0)

    # make pretty pictures and save data
    data_list = [avg_K_plus, avg_K_minus, corr_eps_Kplus, corr_mag_eps_Hplus,
                 corr_eps_Hplus, avg_epsilon, avg_epsilon2, avg_H_plus,
                 avg_H_plus2, avg_H_minus, avg_H_minus2, avg_epsilon_H,
                 avg_total_t]
    name_list = ["avg_K_plus", "avg_K_minus", "corr_eps_Kplus",
                 "corr_mag_eps_Hplus", "corr_eps_Hplus", "avg_epsilon",
                 "avg_epsilon2", "avg_H_plus", "avg_H_plus2", "avg_H_minus",
                 "avg_H_minus2", "avg_epsilon_H", "avg_total_t"]
    for data, name in zip(data_list, name_list):
        plot_maker(dim1vals, dim2vals, data, system, 'comb', .1, -.1, False, name, False, coordsys, scale_dict)
        np.save(path + '/npy/' + system + '.' + name + '.npy', data)
        if coordsys == "polar":
            avg_over_theta(path + '/npy/' + system + '.' + name)


def avg_over_theta(path):
    """
    Compute average of the quantity in question over the theta dimension.

    Parameters
    ----------
    path : string
        The directory in which nougat npy outputs are located

    Returns
    -------
    None.

    """
    data = np.load(path + ".npy")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avg_vals = np.nanmean(data, axis=1)
        std = np.nanstd(data, axis=1)
    np.save(path + ".avg_over_theta.npy", avg_vals)
    np.save(path + ".avg_over_theta.std.npy", std)


def bad_measure_t0(zone, ztwo, coordsys):
    """
    Compute average bulk membrane thickness by measuring thickness at the box \
    border. This is not a good way of doing things and will be deprecated in \
    favor of measuring average thickness of an empty membrane.

    Parameters
    ----------
    zone : numpy array
        data for outer leaflet height
    ztwo : numpy array
        data for inner leaflet height
    coordsys : string
        "polar" or "cart"

    Returns
    -------
    avgt0 : float
        the average thickness of the membrane at the borders of the box

    """
    thickness = zone - ztwo

    avgthickness = calc_avg_over_time(thickness)

    if coordsys == "cart":
        leftcol = np.mean(avgthickness[:, 0])
        rightcol = np.mean(avgthickness[:, -1])
        toprow = np.mean(avgthickness[0, :])
        botrow = np.mean(avgthickness[-1, :])
        avgt0 = (leftcol + rightcol + toprow + botrow) / 4.0
    elif coordsys == "polar":
        avgt0 = np.mean(avgthickness[-1:])

    avgt0 = avgt0 / 2.0

    return avgt0


def measure_quant_in_empty_sys(path, system, coordsys, quantity):
    """
    Measure the average thickness of a membrane.

    Parameters
    ----------
    path : string
        path to the directory where your nougat outputs are
    system : string
        name of the system you gave nougat
    coordsys : string
        "polar" or "cart"; if polar, will ignore small r bins (area too small)
    quantity : string
        The quantity you want to take the average of

    Returns
    -------
    avg : float
        the average quantity of the membrane

    """
    data = np.load(path + '/npy/' + system + '.' + quantity + '.npy')
    if coordsys == "polar":
        data = data[4:, :, :]  # this could be smarter
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avg = np.nanmean(data)
    return avg
