"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings


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


def create_outfile_directories(cwd):
    """
    Create the preliminary directory hierarchy for nougat.py outputs.

    Parameters
    ----------
    cwd  :  Path object
        the path to the current working directory.

    Returns
    -------
    None.

    """
    #quantities = ["height", "density", "curvature", "thickness", "order", "tilt", "misc"]
    quantities = ["height", "curvature", "thickness", "misc"]
    for filetype in ["trajectory", "average", "figures"]:
        for quantity in quantities:
            if quantity == "curvature":
                for curv in ["mean", "gaussian", "normal_vectors"]:
                    dirname = cwd.joinpath(filetype, quantity, curv)
                    dirname.mkdir(parents=True, exist_ok=True)
            else:
                dirname = cwd.joinpath(filetype, quantity)
                dirname.mkdir(parents=True, exist_ok=True)


def parse_dat_file(path, bin_info, quant):
    """
    Parse nougat .dat file and turn it into a numpy ndarray.

    Parameters
    ----------
    path : Pathlib Path object
        The path to the .dat file.
    bin_info : dict
        Dictionary containing bin sizes and number of frames in trajectory.
    quant : str
        The quantity being parsed.

    Returns
    -------
    data_array : numpy ndarray
        3D array with dimensions time, [x, r], [y, theta].

    """
    N1_bins = bin_info["N1"]
    N2_bins = bin_info["N2"]
    Nframes = bin_info["nframes"]

    if quant != "density":
        data = np.genfromtxt(path, missing_values='nan', filling_values=np.nan)
    else:
        data = np.genfromtxt(path, missing_values='nan', filling_values="0")

    # create a new array that has each frame in a different array level
    data_array = np.zeros((Nframes, N1_bins, N2_bins))

    # put each frame in its own level of the matrix
    for frm in range(Nframes):
        data_array[frm, :, :] = data[frm * N1_bins:(frm + 1) * N1_bins, 2:]

    return data_array


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


def filename_generator(sys_name, lipid_name, field, beadname, coordsys, measure, dtype):
    """
    Generate old filenames for legacy uses.

    Parameters
    ----------
    sys_name : str
        System name provided to nougat.tcl.
    lipid_name : str
        The lipid species present.
    field : str
        The surface being analyzed.
    beadname : str
        The beads specified as headnames.
    coordsys : str
        "polar" or "cart".
    measure : str
        That quantity being measured.
    dtype : str
        The file suffix (npy or dat, for example).

    Raises
    ------
    RuntimeWarning
        Raised if wrong str supplied as measure.

    Returns
    -------
    filename : str
        Complete legacy filename.

    """
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


def save_areas(bin_info, min_val, polar):
    """
    Calculate area of each bin and save to file.

    Parameters
    ----------
    bin_info : dict
        Dict containing bin sizes.
    min_val : float
        Distance from origin that first radial bin starts. Only used in polar \
            coordinate systems. Default is zero.
    polar : bool
        Whether to use polar coordinates or cartesian.

    Returns
    -------
    None.

    """
    N1_bins = bin_info["N1"]
    N2_bins = bin_info["N2"]
    d1 = bin_info["d1"]
    d2 = bin_info["d2"]

    areas = np.ones([N1_bins, N2_bins])

    areas = areas * d1 * d2

    if polar:
        for row in range(N1_bins):
            dist_to_center = min_val + row * d1 + d1 / 2.0
            areas[row, :] = areas[row, :] * dist_to_center

    np.save('trajectory/density/areas.npy', areas)


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


def read_log():
    """
    Read log file output by nougat.tcl and save important info for later.

    Parameters
    ----------
    None.

    Returns
    -------
    system_dict : DICT
        A dictionary containing the list of lipid species, their respective \
            headnames, their respective density normalization factors, and \
            the bin sizes.

    """
    system_dict = {}

    # open log file
    with open("tcl_output/nougat.log", "r+") as log_file:
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


def plot_all_quantities(polar, system_dict, cwd, inclusion):
    """
    Plot all heatmaps for height, curvature, thickness, density, and order.

    Parameters
    ----------
    polar : bool
        Whether or not to use polar coordinates instead of cartesian.
    system_dict : dict
        Dictionary containing list of species, bin sizes, and normalization \
            constants.
    cwd : Pathlib Path object
        Path to current working directory.
    inclusion : False or list
        If no inclusion, False. Else, contains list of coordinates to be plotted. \
        This feature not fully implemented.

    Returns
    -------
    None.

    """
    # prep heatmap plot dimensions
    hmap_dims = bin_prep(system_dict['bin_info'], polar)

    for species in system_dict['species']:
        for field in ["zone", "ztwo", "zzero", "zplus"]:
            for quantity in ['height', 'curvature/gaussian', 'curvature/mean']:
                hmap_data = np.genfromtxt(cwd.joinpath("average", quantity, field + ".dat"), delimiter=",")
                fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, quantity, polar)
                plt.savefig(cwd.joinpath("figures", quantity, field + ".pdf"))
                plt.close()
        for field in ["zone", "ztwo", "whole"]:
            hmap_data = np.genfromtxt(cwd.joinpath("average", "thickness", field + ".dat"), delimiter=",")
            fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, "thickness", polar)
            plt.savefig(cwd.joinpath("figures", "thickness", field + ".pdf"))
            plt.close()
        """
        for field in ["zone", "ztwo"]:
            hmap_data = np.genfromtxt(cwd.joinpath("average", "density", species, field + ".dat"), delimiter=",")
            fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, "density", polar)
            plt.savefig(cwd.joinpath("figures", "density", species, field + ".pdf"))
            plt.close()
            for tail in range(system_dict['ntails'][species]):
                hmap_data = np.genfromtxt(cwd.joinpath("average", "order", species, "tail" + str(tail), field + ".dat"), delimiter=",")
                fig, ax = plot_maker(hmap_dims, hmap_data, inclusion, "order", polar)
                plt.savefig(cwd.joinpath("figures", "order", species, "tail" + str(tail), field + ".pdf"))
                plt.close()
        """


def plot_maker(dims, data, protein, quant, polar):
    """
    Create and save 2D heatmaps.

    Parameters
    ----------
    dims : list
        np.meshgrid output
    data : array
        the 2d array/matrix of values to be heatmapped
    protein : list or False
        if --inclusion turned on, list of helix coordinates; if no protein, False
    quant : str
        The quantity being plotted.
    polar : bool
        Whether or not to use polar coordinates

    Returns
    -------
    None.

    """
    dim1vals, dim2vals = dims

    if quant == "density":
        Vmin = 0
        Vmax = 2
    elif quant == "height":
        Vmin = -60
        Vmax = 60
    else:
        Vmin, Vmax = "auto", "auto"

    fig = plt.figure()
    if polar:
        ax = plt.subplot(projection="polar")
    else:
        ax = plt.subplot()
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    create_heatmap(polar, dim1vals, dim2vals, data, Vmax, Vmin, True)
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
    protein : list or False
        if --inclusion turned on, list of helix coordinates; if no protein, False
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

    # This circle is a custom E protein thing and should be removed
    # circle1 = plt.Circle((0, 0), 28.116, transform=ax.transData._b, color='black', linestyle='dashed', linewidth=4, fill=False)
    # if field == "zone":
    #    ax.add_artist(circle1)


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


def dimensions_analyzer(data, polar):
    """
    Determine system dimensions (number of bins, size of bins, number of frames) \
        from .dat file.

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    polar : TYPE
        Whether to use polar coordinates or cartesian.

    Raises
    ------
    Exception
        If the number of rows is not divisible by the number of frames, that \
            means there is a problem in the construction or parsing of the file.

    Returns
    -------
    N1_bins : int
        DESCRIPTION.
    d1 : float
        DESCRIPTION.
    N2_bins : int
        DESCRIPTION.
    d2 : float
        DESCRIPTION.
    Nframes : int
        DESCRIPTION.
    match_value : float
        The first radial bin in polar coords should have match_value of 0 unless \
            you've specified a min value.

    """
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

    if polar is False:
        d1 = data[0, 1] - data[0, 0]
        d2 = (np.pi * 2) / N2_bins
    elif polar is True:
        # compute average d1, assume d2 is the same
        d1list = []
        for row in range(Nframes):
            d1list.append(data[row * N1_bins, 1])
        d1 = np.mean(d1list)
        d2 = d1

    return N1_bins, d1, N2_bins, d2, Nframes, match_value


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
    data = np.load(path.with_suffix(".npy"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avg_vals = np.nanmean(data, axis=1)
        std = np.nanstd(data, axis=1)
    np.save(path.with_suffix(".avg_over_theta.npy"), avg_vals)
    np.save(path.with_suffix(".avg_over_theta.std.npy"), std)


def bad_measure_t0(zone, ztwo, polar):
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
    polar : bool
        Whether or not to use polar coordinates

    Returns
    -------
    avgt0 : float
        the average thickness of the membrane at the borders of the box

    """
    thickness = zone - ztwo

    avgthickness = calc_avg_over_time(thickness)

    if polar is False:
        leftcol = np.mean(avgthickness[:, 0])
        rightcol = np.mean(avgthickness[:, -1])
        toprow = np.mean(avgthickness[0, :])
        botrow = np.mean(avgthickness[-1, :])
        avgt0 = (leftcol + rightcol + toprow + botrow) / 4.0
    elif polar is True:
        avgt0 = np.mean(avgthickness[-1:])

    avgt0 = avgt0 / 2.0

    return avgt0


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
