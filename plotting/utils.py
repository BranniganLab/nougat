"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os


def strip_blank_lines(f):
    for l in f:
        line = l.strip()
        if line:
            yield line


def read_config(path):
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


def find_first_val(l):
    for value in l:
        if np.isnan(value):
            continue
        else:
            return value
    return np.nan


def find_last_val(l):
    last_index = len(l) - 1
    while last_index >= 0:
        if not np.isnan(l[last_index]):
            return l[last_index]
        else:
            last_index -= 1
    return "this didnt work"


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


def mostly_empty(data_array, N1_bins, N2_bins, Nframes):
    # if a bin only has lipids in it <10% of the time, it shouldn't be considered part of the membrane
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
        names_dict['density_norm'] = lines[start_line].split(":")[1]

    return names_dict


def plot_maker(dim1vals, dim2vals, data, name, field, Vmax, Vmin, protein, dataname, bead, coordsys):
    """
    Make 2D heatmaps.

    Parameters
    ----------
    dim1vals : list
        meshgrid output 2
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
        if protein present, list of helix coordinates; if no protein, False
    dataname : string
        the type of measurement (thickness, height, curvature, etc)
    bead : string or False
        if bead specified, name of bead; else False
    coordsys : string
        "polar" or "cart"

    Returns
    -------
    None.

    """
    fig = plt.figure()

    if coordsys == "polar":
        ax = plt.subplot(projection="polar")
        if Vmax != "auto":
            c = plt.pcolormesh(dim2vals, dim1vals, data, cmap="RdBu_r", zorder=0, vmax=Vmax, vmin=Vmin)
        else:
            c = plt.pcolormesh(dim2vals, dim1vals, data, cmap="RdBu_r", zorder=0)
    elif coordsys == "cart":
        ax = plt.subplot()
        if Vmax != "auto":
            c = plt.pcolormesh(dim1vals, dim2vals, data, cmap="RdBu_r", zorder=0, vmax=Vmax, vmin=Vmin)
        else:
            c = plt.pcolormesh(dim1vals, dim2vals, data, cmap="RdBu_r", zorder=0)
    else:
        print("something's wrong with coordsys")

    # cbar = plt.colorbar(c)

    if protein is not False:
        print(protein)
        for i in range(0, 10, 2):
            protein[i + 1] = np.deg2rad(protein[i + 1])
            if coordsys == "cart":
                protein[i], protein[i + 1] = convert_to_cart(protein[i], protein[i + 1])
            plt.scatter(protein[i + 1], protein[i], c="black", linewidth=4, zorder=2)
        circle1 = plt.Circle((0, 0), 28.116, transform=ax.transData._b, color='black', linestyle='dashed', linewidth=4, fill=False)
        if field == "zone":
            ax.add_artist(circle1)

    plt.axis('off')
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # fig.set_size_inches(6,6)

    if bead is False:
        plt.savefig('pdf/' + name + "_" + field + "_" + coordsys + "_" + dataname + ".pdf", dpi=700)
    else:
        plt.savefig('pdf/' + name + "_" + bead + "_" + field + "_" + coordsys + "_" + dataname + ".pdf", dpi=700)
    plt.clf()
    plt.close()


def convert_to_cart(rval, thetaval):
    xval = rval * np.cos(thetaval)
    yval = rval * np.sin(thetaval)
    return xval, yval


def coord_format(value):
    rounded = round(value, 3)
    leftside, rightside = str(rounded).split('.')
    if len(rightside) < 3:
        rightside = rightside + (' ' * (3 - len(rightside)))
    if len(leftside) < 4:
        leftside = (' ' * (4 - len(leftside))) + leftside
    final_value = leftside + '.' + rightside
    return final_value


def bin_format(value):
    strval = str(value)
    length = len(strval)
    final_value = (' ' * (3 - length)) + strval + '.00'
    return final_value


def dimensions_analyzer(data, coordsys):
    # figure out how many radial or x bins there are
    counter = 1
    flag = True
    match_value = data[0, 0]
    while (flag == True):
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


def calc_elastic_terms(system, path, coordsys):
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

    # make pretty pictures and save data
    data_list = [avg_K_plus, avg_K_minus, corr_eps_Kplus, corr_mag_eps_Hplus, corr_eps_Hplus, avg_epsilon, avg_epsilon2, avg_H_plus, avg_H_plus2, avg_H_minus, avg_H_minus2, avg_epsilon_H, avg_total_t]
    name_list = ["avg_K_plus", "avg_K_minus", "corr_eps_Kplus", "corr_mag_eps_Hplus", "corr_eps_Hplus", "avg_epsilon", "avg_epsilon2", "avg_H_plus", "avg_H_plus2", "avg_H_minus", "avg_H_minus2", "avg_epsilon_H", "avg_total_t"]
    for data, name in zip(data_list, name_list):
        plot_maker(dim1vals, dim2vals, data, system, 'comb', .1, -.1, False, name, False, coordsys)
        np.save(path + '/npy/' + system + '.' + name + '.npy', data)


def measure_t0(zone, ztwo, coordsys):
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
