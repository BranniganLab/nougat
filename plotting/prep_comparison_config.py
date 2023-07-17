"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

from itertools import product
from itertools import combinations


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
    templist : LIST
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
            for item in config_dict['TOC'][category].keys():
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


def plot_combination(paths, name, quantity):
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

    Returns
    -------
    None.

    """
    return


if __name__ == "__main__":
    quant_list = ["avg_K_plus", "avg_K_minus", "corr_epst0_Kplus",
                  "corr_mag_epst0_Hplus", "corr_epst0_Hplus", "avg_epsilon",
                  "avg_rms_epsilon_over_t0", "avg_abs_epsilon",
                  "avg_abs_epsilon_over_t0", "avg_epsilon2", "avg_H_plus",
                  "avg_H_plus2", "avg_H_minus", "avg_H_minus2",
                  "avg_epsilon_H", "avg_total_t", "avg_z_minus",
                  "avg_z_minus2", "avg_z_minus_H_minus", "avg_epsilon_over_t0",
                  "avg_epsilon_H_over_t0", "avg_epsilon2_over_t02",
                  "avg_tilde_t", "avg_z_minus2_over_t02",
                  "avg_z_minus_H_minus_over_t0"]
    # prep_config(["Lipid Tail Length", "Saturation", "Structure"], [["2", "3", "4", "5", "6"],["Saturated", "Mono-unsaturated"],["capped","uncapped","protein-less"]])
    config_dict = read_config('comp_config_main.txt', [
                              "length", "saturation", "structure"])
    for combination in generate_combinations(config_dict):
        for quantity in quant_list:
            plot_combination(combination[0], combination[1], quantity)
