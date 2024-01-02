"""Functions relating to calculation of curvature."""
import numpy as np
from utils import *


def calculate_curvature(sys_name, coordsys, system_dict, cwd):
    """
    Calculate mean and Gaussian curvature, as well as normal vectors for each \
        surface.

    Parameters
    ----------
    sys_name : str
        System name you gave nougat.py.
    coordsys : str
        "polar" or "cart".
    system_dict : dict
        Dictionary containing key 'bin_info' that has bin dimensions.
    cwd : PathLib Path object
        Path to current working directory.

    Returns
    -------
    None.

    """
    for field in ["zone", "ztwo", "zzero", "zplus"]:
        field_height = np.load(cwd.joinpath("trajectory", "height", field + ".npy"))

        wrapped_height = make_pbc(field_height, coordsys, system_dict)

        diffs = take_finite_differences(wrapped_height, system_dict)
        r_vector = calculate_r(system_dict['bin_info'], coordsys)
        H = measure_mean_curvature(diffs, r_vector, coordsys)
        K = measure_gaussian_curvature(diffs, r_vector, coordsys)
        Nvecs = measure_normal_vectors(diffs, system_dict['bin_info'], r_vector, coordsys)

        # take the average curvatures over all frames
        avgH = calc_avg_over_time(H)
        avgK = calc_avg_over_time(K)

        # save as files for debugging / analysis
        np.savetxt(cwd.joinpath("average", "curvature", "mean", field + ".dat"), avgH, delimiter=',', fmt='%10.5f')
        np.savetxt(cwd.joinpath("average", "curvature", "gaussian", field + ".dat"), avgK, delimiter=',', fmt='%10.5f')
        np.save(cwd.joinpath("trajectory", "curvature", "mean", field + ".npy"), H)
        np.save(cwd.joinpath("trajectory", "curvature", "gaussian", field + ".npy"), K)
        np.save(cwd.joinpath("average", "curvature", "mean", field + ".npy"), avgH)
        np.save(cwd.joinpath("average", "curvature", "gaussian", field + ".npy"), avgK)
        np.save(cwd.joinpath("trajectory", "curvature", "normal_vectors", field + ".npy"), Nvecs)
        if coordsys == "polar":
            avg_over_theta(cwd.joinpath("average", "curvature", "mean", field))
            avg_over_theta(cwd.joinpath("average", "curvature", "gaussian", field))

        print(sys_name + " " + field + " curvatures done!")


def make_pbc(height, coordsys, system_dict):
    """
    Implement periodic boundary conditions by wrapping the height values.

    Parameters
    ----------
    height : numpy ndarray
        3D matrix containing height values.
    coordsys : str
        "polar" or "cart".
    system_dict : dict
        Dictionary containing bin information in key 'bin_info'.

    Returns
    -------
    wrapped_inputs : numpy ndarray
        Expanded height matrices that allow for calculation of curvature with \
            PBCs.

    """
    N1_bins = system_dict['bin_info']['N1']
    N2_bins = system_dict['bin_info']['N2']
    Nframes = np.shape(height)[0]

    # create arrays for storing curvature data
    wrapped_inputs = np.zeros((Nframes, N1_bins + 2, N2_bins + 2))

    wrapped_inputs[:, 1:(N1_bins + 1), 1:(N2_bins + 1)] = height

    # wrap in column directions
    wrapped_inputs[:, :, 0] = wrapped_inputs[:, :, N2_bins]
    wrapped_inputs[:, :, (N2_bins + 1)] = wrapped_inputs[:, :, 1]

    if coordsys == "polar":
        # set top and bottom row to nan
        wrapped_inputs[:, 0, :] = np.nan
        wrapped_inputs[:, N1_bins + 1, :] = np.nan
    elif coordsys == "cart":
        # if cartesian, wrap in both directions
        wrapped_inputs[:, 0, :] = wrapped_inputs[:, N1_bins, :]
        wrapped_inputs[:, (N1_bins + 1), :] = wrapped_inputs[:, 1, :]
        # and fill in the corners
        wrapped_inputs[:, 0, 0] = wrapped_inputs[:, N1_bins, N2_bins]
        wrapped_inputs[:, N1_bins + 1, N2_bins + 1] = wrapped_inputs[:, 1, 1]
        wrapped_inputs[:, 0, N2_bins + 1] = wrapped_inputs[:, N1_bins, 1]
        wrapped_inputs[:, N1_bins + 1, 0] = wrapped_inputs[:, 1, N2_bins]

    return wrapped_inputs


def calculate_r(bin_info, coordsys):
    """
    Calculate the r value for each bin, measuring at the center of the bin.

    Parameters
    ----------
    bin_info : dict
        Dictionary containing bin numbers and sizes.
    coordsys : str
        "polar" or "cart".

    Returns
    -------
    r : numpy ndarray
        A column vector listing all the r values for a given system.

    """
    if coordsys == "polar":
        dr = bin_info['d1']
        Nr = bin_info['N1']

        # calculate the r values for each bin
        r = np.linspace(dr / 2, (dr / 2) + ((Nr - 1) * dr), Nr)

        # turn row vector into column vector for correct multiplication later
        r = r[:, None]
    elif coordsys == "cart":
        r = np.nan

    return r


def take_finite_differences(heights, system_dict):
    """
    Take 1st and 2nd order finite differences of height field.

    Parameters
    ----------
    heights : numpy ndarray
        height field h(t,a,b) where a/b can be x/y or r/theta, depending on \
            the coordinate system used. t is frames in the trajectory.
    system_dict : dict
        dictionary with key 'bin_info' that holds number and size of bins in \
            each direction.

    Returns
    -------
    h_1 : numpy ndarray
        first order difference along dimension 1 (x/r).
    h_2 : numpy ndarray
        first order difference along dimension 2 (y/theta).
    h_11 : numpy ndarray
        second order difference along dimension 1 (x/r).
    h_22 : numpy ndarray
        second order difference along dimension 2 (y/theta).
    h_12 : numpy ndarray
        second order cross difference along dimensions 1 and 2.

    """
    N1_bins = system_dict['bin_info']['N1']
    d1 = system_dict['bin_info']['d1']
    N2_bins = system_dict['bin_info']['N2']
    d2 = system_dict['bin_info']['d2']

    shift_row_up = heights[:, 0:N1_bins, 1:N2_bins + 1]
    shift_row_down = heights[:, 2:N1_bins + 2, 1:N2_bins + 1]
    shift_col_left = heights[:, 1:N1_bins + 1, 0:N2_bins]
    shift_col_right = heights[:, 1:N1_bins + 1, 2:N2_bins + 2]
    no_shift = heights[:, 1:N1_bins + 1, 1:N2_bins + 1]
    shift_up_left = heights[:, 0:N1_bins, 0:N2_bins]
    # shift_up_right = heights[:, 0:N1_bins, 2:N2_bins + 2]
    shift_down_right = heights[:, 2:N1_bins + 2, 2:N2_bins + 2]
    # shift_down_left = heights[:, 2:N1_bins+2, 0:N2_bins]

    h_1 = (shift_row_down - shift_row_up) / (2 * d1)
    h_11 = (shift_row_up + shift_row_down - 2 * no_shift) / d1**2
    h_2 = (shift_col_right - shift_col_left) / (2 * d2)
    h_22 = (shift_col_left + shift_col_right - 2 * no_shift) / d2**2
    # h_12 = (shift_up_left + shift_down_right - shift_up_right - shift_down_left) / (4 * d1 * d2)
    h_12 = (shift_up_left + shift_down_right + 2 * no_shift - shift_col_left - shift_col_right - shift_row_up - shift_row_down) / (2 * d1 * d2)

    return [h_1, h_2, h_11, h_22, h_12]


def measure_mean_curvature(finite_differences, r, coordsys):
    """
    Measure mean curvature H of a membrane field.

    Parameters
    ----------
    finite_differences : list
        List containing the finite differences generated by take_finite_differences.
    r : numpy ndarray
        A column vector listing all the r values for a given system.
    coordsys : str
        "cart" or "polar".

    Returns
    -------
    H : numpy ndarray
        Mean curvature H of the field being measured.

    """
    h_1, _, h_11, h_22, _ = finite_differences

    if coordsys == "polar":
        # mean curvature: 1/2 * [h_rr + 1/r(h_r) + 1/r**2(h_thetatheta)]
        H = (h_11 + r**(-1) * h_1 + r**(-2) * h_22) / 2.0

    elif coordsys == "cart":
        # mean curvature: (Hxx + Hyy)/2
        H = (h_11 + h_22) / 2.0

    else:
        print("something is wrong with coordsys")

    return H


def measure_gaussian_curvature(finite_differences, r, coordsys):
    """
    Measure Gaussian curvature K of a membrane field.

    Parameters
    ----------
    finite_differences : list
        List containing the finite differences generated by take_finite_differences.
    r : numpy ndarray
        A column vector listing all the r values for a given system.
    coordsys : str
        "cart" or "polar".

    Returns
    -------
    K : numpy ndarray
        Gaussian curvature K of the field being measured.

    """
    h_1, h_2, h_11, h_22, h_12 = finite_differences

    if coordsys == "polar":
        # Gaussian curvature: 1/r(h_r*h_rr) + 2/r^3(h_rtheta*h_theta) - 1/r^4(h_theta^2) - 1/r^2(h_rtheta^2 - h_rr*h_thetatheta)
        K = r**(-1) * (h_1 * h_11) + 2 * r**(-3) * (h_12 * h_2) - r**(-4) * h_2**2 - r**(-2) * (h_12**2 - h_11 * h_22)

    elif coordsys == "cart":
        # mean curvature: HxxHyy - Hxy^2
        K = h_22 * h_11 - h_12**2

    else:
        print("something is wrong with coordsys")

    return K


def measure_normal_vectors(finite_differences, bin_info, r, coordsys):
    """
    Measure Gaussian curvature K of a membrane field.

    Parameters
    ----------
    finite_differences : list
        List containing the finite differences generated by take_finite_differences.
    bin_info : dict
        Dictionary containing bin numbers and sizes for each dimension.
    r : numpy ndarray
        A column vector listing all the r values for a given system.
    coordsys : str
        "cart" or "polar".

    Returns
    -------
    Nvecs : numpy ndarray
        Normal vectors of the field being measured.

    """
    h_1, h_2, _, _, _ = finite_differences
    N2 = bin_info['N2']
    d2 = bin_info['d2']

    if coordsys == "polar":
        # Determine theta value for each column
        theta = np.linspace(d2 / 2, d2 / 2 + d2 * (N2 - 1), N2)

        # normal vector: [(((1/r)*sin(theta)*h_theta) - (cos(theta)*h_r))/N, (((-1/r)*cos(theta)*h_theta) - (sin(theta)*h_r))/N, 1/N]
        # N = normalization constant = sqrt(r^(-2) * h_theta^2 + h_r^2 + 1)
        N = np.sqrt(r**(-2) * h_2**2 + h_1**2 + 1)
        Nx = (r**(-1) * np.sin(theta) * h_2 - np.cos(theta) * h_1) * N**(-1)
        Ny = (-1 * r**(-1) * np.cos(theta) * h_2 - np.sin(theta) * h_1) * N**(-1)
        Nz = N**(-1)

    elif coordsys == "cart":
        # normal vector: [-h_r/N, -h_y/N, 1/N]
        # N = normalization constant = sqrt(h_x^2 + h_y^2 + 1)
        N = np.sqrt(h_2**2 + h_1**2 + 1)
        Nx = -1 * h_2 * N**(-1)
        Ny = -1 * h_1 * N**(-1)
        Nz = N**(-1)
    else:
        print("something is wrong with coordsys")

    return [Nx, Ny, Nz]
