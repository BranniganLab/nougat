"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

from itertools import product
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import os
from utils import *


def prep_config(category_names, lol_of_categories):
    """
    Create a text file that the user can manually add paths to.

    Parameters
    ----------
    category_names : LIST
        a list of strings containing category names.
    lol_of_categories : LIST
        list of lists (lol) containing the entries for each category.

    Returns
    -------
    None.

    """
    # error check
    if len(category_names) != len(lol_of_categories):
        print("different number of categories than lists!")
        return

    with open('comp_config.txt', 'w') as f:

        # write comments at top
        f.write(
            "### Please enter paths to nougat-generated output folder after \
                each entry below\n")
        f.write("### If no path exists, write NULL\n")
        f.write("### Entry format is " + ';'.join(category_names) + ":\
                [path]\n")
        f.write("\n")

        # create entries for user to fill in
        for item in product(*lol_of_categories):
            f.write(";".join(item) + ':\n')


def strip_blank_lines(file_lines):
    """
    Generate non-blank lines.

    Parameters
    ----------
    file_lines : LIST
        the list of lines.

    Yields
    ------
    non_blank : STRING
        any non-blank line.

    """
    for line in file_lines:
        non_blank = line.strip()
        if non_blank:
            yield non_blank


def read_config(config_file, category_names):
    """
    Read in the config file provided by user and return a dict.

    Parameters
    ----------
    config_file : STRING
        config filename (must be in working directory)
    category_names : LIST
        a list of strings containing category names

    Returns
    -------
    config_dict : DICT
        dict with paths and other user-specified values

    """
    config_dict = {}
    config_dict['TOC'] = {}
    counter = 0

    with open(config_file, "r+") as f:
        # skip blank lines
        for line in strip_blank_lines(f):

            # ignore lines starting with comment
            if line.startswith("#") is True:
                continue

            else:
                # ignore in-line comment if present
                line = line.partition('#')[0]

                # double-check that is probably unnecessary
                if line.rstrip():

                    # initialize dict key to contain nested dict
                    config_dict[str(counter)] = {}

                    # split the line into 2 parts
                    values = line.split(":")
                    path = values[1]
                    cat_vals = values[0].split(";")

                    # save to nested dict
                    for key, val in zip(category_names, cat_vals):
                        config_dict[str(counter)][str(key)] = val
                        if key in config_dict['TOC']:
                            if val not in config_dict['TOC'][str(key)]:
                                config_dict['TOC'][str(key)].append(str(val))
                        else:
                            config_dict['TOC'][str(key)] = [str(val)]

                    config_dict[str(counter)]["path"] = path

                    counter = counter + 1
    return config_dict


def generate_combinations(config_dict):
    """
    Generate list of the different paths needed to make each figure, as well \
        as the name of the combination type.

    Parameters
    ----------
    config_dict : DICT
        the dictionary made by read_config

    Yields
    ------
    templist : TUPLE
        list of paths corresponding to the correct combination of variables.
    name : STRING
        the name corresponding to that combination of variables.

    """
    key_list = config_dict.keys()
    category_list = list(config_dict['TOC'].keys())
    if len(category_list) == 3:
        for category_combo in combinations(category_list, 2):
            for category1 in config_dict['TOC'][category_combo[0]]:
                for category2 in config_dict['TOC'][category_combo[1]]:
                    templist = []
                    for key in key_list:
                        if key == 'TOC':
                            continue
                        else:
                            if ((config_dict[key][category_combo[0]] == category1) and (config_dict[key][category_combo[1]] == category2)):
                                if config_dict[key]['path'] != "NULL":
                                    templist.append(config_dict[key]['path'])
                    name = category_combo[0] + "_" + category1 + \
                        "_" + category_combo[1] + "_" + category2
                    yield (templist, name)
    elif len(category_list) == 2:
        for category in category_list:
            for item in config_dict['TOC'][category]:
                templist = []
                for key in key_list:
                    if key == 'TOC':
                        continue
                    else:
                        if config_dict[key][category] == item:
                            if config_dict[key]['path'] != "NULL":
                                templist.append(config_dict[key]['path'])
                name = category + "_" + item
                yield (templist, name)
    elif len(category_list) == 1:
        for key in key_list:
            templist = []
            if key == 'TOC':
                continue
            else:
                if config_dict[key]['path'] != "NULL":
                    templist.append(config_dict[key]['path'])
        name = "all"
        yield (templist, name)
    else:
        print("I only made this to work with 3 dimensions - sorry!")


def normalize_by_same_quantity_in_empty_membrane(path, quantity, sysname, species_name):
    """
    Calculate a more complicated quantity that requires custom work.

    Parameters
    ----------
    path : STRING
        The path to the files you are plotting
    quantity : STRING
        The name of the quantity needing to be calculated.
    sysname : STRING
        The name given to nougat.py.

    Returns
    -------
    values : NDARRAY
        The values corresponding to the y axis on the figure you are plotting.

    """
    if species_name in ["lgDT", "lgDL", "lgDP", "lgDB", "lgDX"]:
        empty_sims_path = "/home/js2746/Bending/PC/whole_mols/empty/" + species_name + "/" + species_name + "_polar_5_10_0_-1_1"
    elif species_name in ["lgDY", "lgDO", "lgDG"]:
        empty_sims_path = "/home/js2746/Bending/PC/whole_mols/empty/" + species_name + "/" + species_name + "_polar_10_10_100_-1_1"
    if quantity == "avg_tilde_total_t":
        exp_quantity = "total_t"
        rms = False
    elif quantity == "avg_tilde_epsilon2":
        exp_quantity = "epsilon2"
        rms = True
    elif quantity == "avg_tilde_H_plus2":
        exp_quantity = "H_plus2"
        rms = True
    else:
        print("I think you spelled something wrong")
    exp_value = np.load(path + "/npy/" + sysname + "." + exp_quantity + ".npy")
    bulk_avg = measure_quant_in_empty_sys(empty_sims_path, species_name, "polar", exp_quantity)
    normed_values = calc_avg_over_time(exp_value / bulk_avg)
    if rms is True:
        normed_values = np.sqrt(normed_values)
    np.save(path + "/npy/" + sysname + "." + quantity + ".npy", normed_values)
    avg_over_theta(path + "/npy/" + sysname + "." + quantity)
    return np.load(path + "/npy/" + sysname + "." + quantity + ".avg_over_theta.npy")


def calc_eps_t0(path, quantity, sysname, species_name):
    """
    Calculate epsilon over t0.

    Parameters
    ----------
    path : STRING
        The path to the files you are plotting
    quantity : STRING
        The name of the quantity needing to be calculated.
    sysname : STRING
        The name given to nougat.py.

    Returns
    -------
    values : NDARRAY
        The values corresponding to the y axis on the figure you are plotting.

    """
    if species_name in ["lgDT", "lgDL", "lgDP", "lgDB", "lgDX"]:
        empty_sims_path = "/home/js2746/Bending/PC/whole_mols/empty/" + species_name + "/" + species_name + "_polar_5_10_0_-1_1"
    elif species_name in ["lgDY", "lgDO", "lgDG"]:
        empty_sims_path = "/home/js2746/Bending/PC/whole_mols/empty/" + species_name + "/" + species_name + "_polar_10_10_100_-1_1"
    exp_value = np.load(path + "/npy/" + sysname + ".epsilon.npy")
    bulk_avg = measure_quant_in_empty_sys(empty_sims_path, species_name, "polar", "total_t")
    normed_values = calc_avg_over_time(exp_value / bulk_avg)
    np.save(path + "/npy/" + sysname + "." + quantity + ".npy", normed_values)
    avg_over_theta(path + "/npy/" + sysname + "." + quantity)
    return np.load(path + "/npy/" + sysname + "." + quantity + ".avg_over_theta.npy")


def plot_combination(paths, name, quantity, stds, rmin):
    """
    Plot the quantity specified from each of the paths on the same figure.

    Parameters
    ----------
    paths : LIST
        a list of paths where the correct nougat results can be found
    name : STRING
        the name of the figure (output from generate_combinations)
    quantity : STRING
        the name of the values being plotted (e.g. height, thickness, etc.)
    stds : BOOLEAN
        Whether or not you want standard deviation shown on your plots
    rmin : FLOAT
        The r value below which no line should be plotted.

    Returns
    -------
    None.

    """
    fig, axs = plt.subplots(layout="constrained")
    fig.supxlabel(r'$r \;(\mathrm{nm})$')
    fig.supylabel(y_label_dict[quantity])
    for path in paths:
        # find the correct system name
        nougval = [i for i in path.split("/") if "polar" in i][0]
        sysname = nougval.split("_polar")[0]
        if "_" in sysname:
            species_name = sysname.split("_")[0]
        else:
            species_name = sysname

        if quantity in ["avg_tilde_total_t", "avg_tilde_epsilon2", "avg_tilde_H_plus2"]:
            y_vals = normalize_by_same_quantity_in_empty_membrane(path, quantity, sysname, species_name)
        elif quantity == "avg_epsilon_over_t0":
            y_vals = calc_eps_t0(path, quantity, sysname, species_name)
        else:
            y_vals = np.load(path + "/npy/" + sysname + "." + quantity + ".avg_over_theta.npy")

        # figure out what the x axis values should be
        tcl_output = np.genfromtxt(path + '/tcl_output/' + sysname + '.zone.C1A.C1B.polar.height.dat',
                                   missing_values='nan', filling_values=np.nan)
        Nr, dr, _, _, _, _ = dimensions_analyzer(tcl_output, "polar")
        xmin = dr / 2
        xmax = Nr * dr - xmin
        x_vals = np.linspace(xmin, xmax, Nr) / 10
        flag = True
        i = 0
        if rmin > x_vals[i]:
            while flag is True and i < len(x_vals):
                i += 1
                if rmin < x_vals[i]:
                    x_vals = x_vals[i:]
                    y_vals = y_vals[i:]
                    flag = False

        axs.plot(x_vals, y_vals, color=color_dict[species_name], linestyle=style_dict[species_name])
        if "tilde" in quantity:
            axs.axhline(1, color="gray", linestyle="")
        if stds is True:
            std_data = np.load(path + "/npy/" + sysname + "." + quantity + ".avg_over_theta.std.npy")
            std_data = std_data[i:] / np.sqrt(10)
            axs.fill_between(x_vals, (y_vals - std_data), (y_vals + std_data), alpha=.4, color=color_dict[species_name])
        axs.set_xlim(0, xmax / 10)
        if quantity + "_min" in scale_dict:
            axs.set_ylim(scale_dict[quantity + "_min"], scale_dict[quantity + "_max"])
    if stds is True:
        plt.savefig(name + "_" + quantity + "_with_stdv.pdf", dpi=700)
    else:
        plt.savefig(name + "_" + quantity + ".pdf", dpi=700)
    plt.clf()
    plt.close()


y_label_dict = {
    "avg_epsilon_over_t0": r'$\langle \epsilon / t_0 \rangle$',
    "avg_abs_epsilon": r'$\langle | \epsilon | \rangle\; (\mathrm{\dot A})$',
    "avg_abs_epsilon_over_t0": r'$\langle | \epsilon / t_0 | \rangle$',
    "avg_epsilon2_over_t02": r'$\langle ( \epsilon / t_0 )^2 \rangle$',
    "avg_epsilon_H_over_t0": r'$\langle \epsilon H^+ / t_0 \rangle\; (\mathrm{\dot A^{-1}})$',
    "avg_epsilon2": r'$\langle \epsilon^2 \rangle\; (\mathrm{\dot A^2})$',
    "avg_tilde_epsilon2": r'$\langle \sqrt{\tilde \epsilon ^ 2} \rangle$',
    "avg_H_plus2": r'$\langle ( H^+ )^2 \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_tilde_H_plus2": r'$\langle \sqrt{\tilde H_+ ^ 2} \rangle$',
    "avg_tilde_total_t": r'$\langle \tilde t \rangle$',
    "avg_epsilon": r'$\langle \epsilon \rangle\; (\mathrm{\dot A})$',
    "avg_total_t": r'$\langle t \rangle\; (\mathrm{\dot A})$',
    "corr_mag_epst0_Hplus": r'$\langle | \epsilon| | H^+ | / t_0 \rangle - \langle |\epsilon| / t_0 \rangle \langle |H^+ | \rangle\; (\mathrm{\dot A^{-1}})$',
    "corr_epst0_Hplus": r'$ \langle \epsilon H^+ / t_0 \rangle - \langle \epsilon/ t_0 \rangle \langle H^+ \rangle \; ( \mathrm{\dot A^{-1}} )$',
    "avg_rms_epsilon_over_t0": r'$\langle \mathrm{rms}\;\epsilon / t_0 \rangle$',
    "avg_K_plus": r'$\langle K^+ \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_K_minus": r'$\langle K^- \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_H_plus": r'$\langle H^+ \rangle\; (\mathrm{\dot A^{-1}})$',
    "avg_H_minus": r'$\langle H^- \rangle\; (\mathrm{\dot A^{-1}})$',
    "avg_H_minus2": r'$\langle \left ( H^- \right )^2 \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_epsilon_H": r'$\langle  \epsilon H^+  \rangle$',
    "avg_z_minus": r'$\langle z^- \rangle\; (\mathrm{\dot A})$',
    "avg_z_minus2": r'$\langle \left ( z^- \right )^2 \rangle\; (\mathrm{\dot A^2})$',
    "avg_z_minus_H_minus": r'$\langle z^- H^- \rangle$',
    "avg_z_minus2_over_t02": r'$\langle \left ( z^- / t_0 \right )^2 \rangle$',
    "avg_z_minus_H_minus_over_t0": r'$\langle z^- H^- / t_0 \rangle\; (\mathrm{\dot A^{-1}})$',
    "corr_epst0_Kplus": r'$\langle \epsilon K^+ / t_0 \rangle - \langle \epsilon / t_0 \rangle \langle K^+ \rangle\; (\mathrm{\dot A^{-2}})$',
    "corr_eps_Kplus": r'$\langle \epsilon K^+ \rangle - \langle \epsilon \rangle \langle K^+ \rangle\; (\mathrm{\dot A^{-1}})$',
    "corr_mag_eps_Hplus": r'$\langle | \epsilon | | H^+ | \rangle - \langle | \epsilon | \rangle \langle | H^+ | \rangle\; (\mathrm{\dot A^{-2}})$',
    "corr_eps_Hplus": r'$\langle  \epsilon H^+ \rangle - \langle  \epsilon  \rangle \langle  H^+  \rangle\; (\mathrm{\dot A^{-2}})$'
}

scale_dict = {
    "avg_K_plus_min": -.0001,
    "avg_K_plus_max": .0001,
    "avg_tilde_epsilon2_min": 0,
    "avg_tilde_epsilon2_max": 18
}

color_dict = {
    "lgDT": "red",
    "lgDL": "orange",
    "lgDY": "orange",
    "lgDP": "green",
    "lgPO": "green",
    "lgDO": "green",
    "lgDB": "blue",
    "lgDG": "blue",
    "lgDX": "purple",
    "DTPC": "red",
    "DLPC": "orange",
    "DYPC": "orange",
    "DPPC": "green",
    "POPC": "green",
    "DOPC": "green",
    "DBPC": "blue",
    "DGPC": "blue",
    "DXPC": "purple"
}

style_dict = {
    "lgDT": "solid",
    "lgDL": "solid",
    "lgDY": "dashed",
    "lgDP": "solid",
    "lgPO": "dotted",
    "lgDO": "dashed",
    "lgDB": "solid",
    "lgDG": "dashed",
    "lgDX": "solid",
    "DTPC": "solid",
    "DLPC": "solid",
    "DYPC": "dashed",
    "DPPC": "solid",
    "POPC": "dotted",
    "DOPC": "dashed",
    "DBPC": "solid",
    "DGPC": "dashed",
    "DXPC": "solid"
}

if __name__ == "__main__":
    quant_list1 = ["avg_K_plus", "avg_K_minus", "corr_eps_Kplus",
                   "corr_mag_eps_Hplus", "corr_eps_Hplus", "avg_epsilon",
                   "avg_epsilon2", "avg_H_plus", "avg_H_plus2", "avg_H_minus",
                   "avg_H_minus2", "avg_epsilon_H", "avg_total_t", "avg_epsilon_over_t0", "avg_tilde_total_t", "avg_tilde_epsilon2", "avg_tilde_H_plus2"]
    # prep_config(["Lipid Tail Length", "Saturation"], [["2", "3", "4", "5", "6"], ["Saturated", "Mono-unsaturated"]])

    config_dict = read_config('comp_config.txt', ["length", "saturation"])
    cwd = os.getcwd()
    try:
        os.mkdir("avg_over_theta_comparisons")
    except:
        pass
    os.chdir("avg_over_theta_comparisons")
    for combination in generate_combinations(config_dict):
        for quantity in quant_list1:
            plot_combination(combination[0], combination[1], quantity, True, 2.75)
            plot_combination(combination[0], combination[1], quantity, False, 2.75)
    os.chdir(cwd)
