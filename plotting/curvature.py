"""Functions relating to calculation of curvature."""
import numpy as np
from utils import *


def init_curvature_data(height, polar, system_dict):
    N1_bins = system_dict['bin_info']['N1']
    N2_bins = system_dict['bin_info']['N2']
    Nframes = np.shape(height)[0]

    # create arrays for storing curvature data
    if polar is True:
        curvature_inputs = np.zeros((Nframes, N1_bins, N2_bins + 2))
        curvature_outputs = np.zeros((Nframes, N1_bins, N2_bins + 2))
        kgauss_outputs = np.zeros((Nframes, N1_bins, N2_bins + 2))
        normal_vector_outputs = np.zeros((Nframes, N1_bins, 3 * (N2_bins + 2)))
    elif polar is False:
        curvature_inputs = np.zeros((Nframes, N1_bins + 2, N2_bins + 2))
        curvature_outputs = np.zeros((Nframes, N1_bins + 2, N2_bins + 2))
        kgauss_outputs = np.zeros((Nframes, N1_bins + 2, N2_bins + 2))
        normal_vector_outputs = np.zeros((Nframes, N1_bins + 2, 3 * (N2_bins + 2)))

    if polar is True:
        # wrap the inputs in the theta direction for calculating curvature
        curvature_inputs[:, :, 1:(N2_bins + 1)] = height
        curvature_inputs[:, :, 0] = curvature_inputs[:, :, N2_bins]
        curvature_inputs[:, :, (N2_bins + 1)] = curvature_inputs[:, :, 1]
    elif polar is False:
        # if cartesian, wrap in both directions
        curvature_inputs[:, 1:(N1_bins + 1), 1:(N2_bins + 1)] = height
        curvature_inputs[:, :, 0] = curvature_inputs[:, :, N2_bins]
        curvature_inputs[:, :, (N2_bins + 1)] = curvature_inputs[:, :, 1]
        curvature_inputs[:, 0, :] = curvature_inputs[:, N1_bins, :]
        curvature_inputs[:, (N1_bins + 1), :] = curvature_inputs[:, 1, :]
        # and fill in the corners
        curvature_inputs[:, 0, 0] = curvature_inputs[:, N1_bins, N2_bins]
        curvature_inputs[:, N1_bins + 1, N2_bins + 1] = curvature_inputs[:, 1, 1]
        curvature_inputs[:, 0, N2_bins + 1] = curvature_inputs[:, N1_bins, 1]
        curvature_inputs[:, N1_bins + 1, 0] = curvature_inputs[:, 1, N2_bins]

    return curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs


def calculate_curvature(sys_name, bead, coordsys, inclusion, polar, dims, field_list, scale_dict, system_dict):
    N1_bins = system_dict['bin_info']['N1']
    N2_bins = system_dict['bin_info']['N2']
    dim1vals, dim2vals = dims

    leaflist = field_list.copy()
    leaflist.append("zplus")

    for field in leaflist:
        field_height = np.load('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.height.npy')

        curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs = init_curvature_data(field_height, polar, system_dict)

        # if a bin is empty, you can't (nicely) measure the curvature of its neighbors
        nan_test, knan_test = empty_neighbor_test(curvature_inputs)

        # measure the laplacian and gaussian curvatures
        if polar is True:
            curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_polar(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, system_dict)
            # diffs = take_finite_differences(curvature_inputs, system_dict)
            # H = measure_mean_curvature(diffs, system_dict['bin_info'], "polar")
        elif polar is False:
            curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_cart(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, system_dict)
            # diffs = take_finite_differences(curvature_inputs, system_dict)
            # H = measure_mean_curvature(diffs, system_dict['bin_info'], "polar")
        # print(np.allclose(curvature_outputs, H, rtol=0, atol=1e-11, equal_nan=True), field)

        # unwrap along dim2 direction
        meancurvature = curvature_outputs[:, :, 1:N2_bins + 1]
        kcurvature = kgauss_outputs[:, :, 1:N2_bins + 1]
        normal_vectors = normal_vector_outputs[:, :, 3:3 * (N2_bins + 1)]

        # if cartesian, unwrap along dim1 direction too
        if polar is False:
            meancurvature = meancurvature[:, :, 1:N1_bins + 1]
            kcurvature = kcurvature[:, :, 1:N1_bins + 1]
            normal_vectors = normal_vectors[:, :, 1:N1_bins + 1]

        # take the average curvatures over all frames
        avgcurvature = calc_avg_over_time(meancurvature)
        avgkcurvature = calc_avg_over_time(kcurvature)

        # make plots!
        plot_maker(dim1vals, dim2vals, avgkcurvature, sys_name, field, scale_dict["gauss_curv_max"], scale_dict["gauss_curv_min"], inclusion, "gausscurvature", bead, coordsys, scale_dict)
        plot_maker(dim1vals, dim2vals, avgcurvature, sys_name, field, scale_dict["mean_curv_max"], scale_dict["mean_curv_min"], inclusion, "curvature", bead, coordsys, scale_dict)

        # save as files for debugging / analysis
        np.savetxt('dat/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.avgcurvature.dat', avgcurvature, delimiter=',', fmt='%10.5f')
        np.savetxt('dat/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.avgKcurvature.dat', avgkcurvature, delimiter=',', fmt='%10.5f')
        np.save('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.meancurvature.npy', meancurvature)
        np.save('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.gausscurvature.npy', kcurvature)
        np.save('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.avgmeancurvature.npy', avgcurvature)
        np.save('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.avggausscurvature.npy', avgkcurvature)
        np.save('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.normal_vectors.npy', normal_vectors)
        if polar is True:
            avg_over_theta('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.avgmeancurvature')
            avg_over_theta('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.avggausscurvature')

        print(sys_name + ' ' + bead + ' ' + field + " curvatures done!")


def measure_curvature_cart(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, system_dict):
    N1_bins = system_dict['bin_info']['N1']
    d1 = system_dict['bin_info']['d1']
    N2_bins = system_dict['bin_info']['N2']
    d2 = system_dict['bin_info']['d2']
    Nframes = np.shape(curvature_inputs)[0]

    # mean curvature: Hxx + Hyy
    # gaussian curvature: HxxHyy - Hxy^2

    for frm in range(Nframes):
        for row in range(N1_bins + 2):
            for col in range(N2_bins + 2):
                if knan_test[frm, row, col] == False:

                    del2x = curvature_inputs[frm, row - 1, col] + curvature_inputs[frm, row + 1, col] - 2 * curvature_inputs[frm, row, col]
                    del2x = del2x / (d1**2)

                    del2y = curvature_inputs[frm, row, col - 1] + curvature_inputs[frm, row, col + 1] - 2 * curvature_inputs[frm, row, col]
                    del2y = del2y / (d2**2)

                    delxy = (curvature_inputs[frm, row + 1, col + 1] - curvature_inputs[frm, row + 1, col] - curvature_inputs[frm, row, col + 1] + 2 * curvature_inputs[frm, row, col] - curvature_inputs[frm, row - 1, col] - curvature_inputs[frm, row, col - 1] + curvature_inputs[frm, row - 1, col - 1])
                    delxy = delxy / (2 * d1 * d2)

                    # delxy = curvature_inputs[frm, row+1,col+1] - curvature_inputs[frm, row+1,col-1] - curvature_inputs[frm, row-1,col+1] + curvature_inputs[frm, row-1,col-1]
                    # delxy = delxy / (4*d1*d2)

                    delx = (curvature_inputs[frm, row + 1, col] - curvature_inputs[frm, row - 1, col]) / (2 * d1)

                    dely = (curvature_inputs[frm, row, col + 1] - curvature_inputs[frm, row, col - 1]) / (2 * d2)

                    normalization_factor = np.sqrt(1 + delx**2 + dely**2)
                    norm_vec_x = -1 * delx / normalization_factor
                    norm_vec_y = -1 * dely / normalization_factor
                    norm_vec_z = 1 / normalization_factor

                    curvature_outputs[frm, row, col] = (del2x + del2y) / 2.0
                    kgauss_outputs[frm, row, col] = del2x * del2y - delxy**2
                    normal_vector_outputs[frm, row, col * 3] = norm_vec_x
                    normal_vector_outputs[frm, row, col * 3 + 1] = norm_vec_y
                    normal_vector_outputs[frm, row, col * 3 + 2] = norm_vec_z

                elif nan_test[frm, row, col] == False:
                    del2x = curvature_inputs[frm, row - 1, col] + curvature_inputs[frm, row + 1, col] - 2 * curvature_inputs[frm, row, col]
                    del2x = del2x / d1**2

                    del2y = curvature_inputs[frm, row, col - 1] + curvature_inputs[frm, row, col + 1] - 2 * curvature_inputs[frm, row, col]
                    del2y = del2y / d2**2

                    delx = (curvature_inputs[frm, row + 1, col] - curvature_inputs[frm, row - 1, col]) / (2 * d1)

                    dely = (curvature_inputs[frm, row, col + 1] - curvature_inputs[frm, row, col - 1]) / (2 * d2)

                    normalization_factor = np.sqrt(1 + delx**2 + dely**2)
                    norm_vec_x = -1 * delx / normalization_factor
                    norm_vec_y = -1 * dely / normalization_factor
                    norm_vec_z = 1 / normalization_factor

                    curvature_outputs[frm, row, col] = (del2x + del2y) / 2.0
                    kgauss_outputs[frm, row, col] = np.nan
                    normal_vector_outputs[frm, row, col * 3] = norm_vec_x
                    normal_vector_outputs[frm, row, col * 3 + 1] = norm_vec_y
                    normal_vector_outputs[frm, row, col * 3 + 2] = norm_vec_z

                else:

                    curvature_outputs[frm, row, col] = np.nan
                    kgauss_outputs[frm, row, col] = np.nan
                    normal_vector_outputs[frm, row, col * 3] = np.nan
                    normal_vector_outputs[frm, row, col * 3 + 1] = np.nan
                    normal_vector_outputs[frm, row, col * 3 + 2] = np.nan

    return curvature_outputs, kgauss_outputs, normal_vector_outputs


def take_finite_differences(curvature_inputs, system_dict):
    """
    Take 1st and 2nd order finite differences of height field.

    Parameters
    ----------
    curvature_inputs : numpy ndarray
        height field h(a,b,t) where a/b can be x/y or r/theta, depending on \
            coordinate system.
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

    up = curvature_inputs[:, 0:N1_bins - 2, 1:N2_bins - 1]
    down = curvature_inputs[:, 2:N1_bins, 1:N2_bins - 1]
    left = curvature_inputs[:, 1:N1_bins - 1, 0:N2_bins - 2]
    right = curvature_inputs[:, 1:N1_bins - 1, 2:N2_bins]
    center = curvature_inputs[:, 1:N1_bins - 1, 1:N2_bins - 1]
    up_left = curvature_inputs[:, 0:N1_bins - 2, 0:N2_bins - 2]
    # up_right = curvature_inputs[:, 0:N1_bins - 2, 2:N2_bins]
    down_right = curvature_inputs[:, 2:N1_bins, 2:N2_bins]
    # down_left = curvature_inputs[:, 2:N1_bins, 0:N2_bins - 2]

    h_1 = (up - down) / (2 * d1)
    h_11 = (up + down - 2 * center) / d1**2
    h_2 = (left - right) / (2 * d2)
    h_22 = (left + right - 2 * center) / d2**2
    # h_12 = (up_left + down_right - up_right - down_left) / (4 * d1 * d2)
    h_12 = (up_left + down_right + 2 * center - left - right - up - down) / (2 * d1 * d2)

    return [h_1, h_2, h_11, h_22, h_12]


def measure_mean_curvature(finite_differences, bin_info, coordsys):
    """
    Measure mean curvature H of a membrane field.

    Parameters
    ----------
    finite_differences : list
        List containing the finite differences generated by take_finite_differences.
    bin_info : dict
        Dictionary containing bin numbers and sizes.
    coordsys : str
        "cart" or "polar".

    Returns
    -------
    H : numpy ndarray
        Mean curvature H of the field being measured.

    """
    h_1, _, h_11, h_22, _ = finite_differences
    dr = bin_info['d1']
    Nr = bin_info['N1']

    if coordsys == "polar":
        # mean curvature: 1/2 * [h_rr + 1/r(h_r) + 1/r**2(h_thetatheta)]
        r = np.linspace(dr + dr / 2, dr * (Nr - 2) + (dr / 2), Nr - 2)
        r = r[:, None]  # turns row vector into column vector
        print(np.shape(r**(-1)), np.shape(r**(-1) * h_11))

        H = 0.5 * (h_11 + r**(-1) * h_1 + r**(-2) * h_22)

    elif coordsys == "cart":
        # mean curvature: (Hxx + Hyy)/2
        H = 0.5 * (h_11 + h_22)

    return H


def measure_curvature_polar(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, system_dict):
    N1_bins = system_dict['bin_info']['N1']
    d1 = system_dict['bin_info']['d1']
    N2_bins = system_dict['bin_info']['N2']
    d2 = system_dict['bin_info']['d2']
    Nframes = np.shape(curvature_inputs)[0]

    # mean curvature: 1/2 * [h_rr + 1/r(h_r) + 1/r**2(h_thetatheta)]
    # gaussian curvature: 1/r(h_r*h_rr) + 2/r**3(h_rtheta*h_theta) - 1/r**4(h_theta**2) - 1/r**2(h_rtheta**2 - h_rr*h_thetatheta)

    for frm in range(Nframes):
        for row in range(N1_bins):
            for col in range(N2_bins + 2):
                if knan_test[frm, row, col] == False:

                    # calculate d2h/dr2
                    del2r = curvature_inputs[frm, row - 1, col] + curvature_inputs[frm, row + 1, col] - 2 * curvature_inputs[frm, row, col]
                    del2r = del2r / d1**2

                    # calculate dh/dr
                    delr = (curvature_inputs[frm, row + 1, col] - curvature_inputs[frm, row - 1, col]) / (2 * d1)

                    # calculate d2h/drdtheta
                    delrdeltheta = (curvature_inputs[frm, row + 1, col + 1] - curvature_inputs[frm, row + 1, col] - curvature_inputs[frm, row, col + 1] + 2 * curvature_inputs[frm, row, col] - curvature_inputs[frm, row - 1, col] - curvature_inputs[frm, row, col - 1] + curvature_inputs[frm, row - 1, col - 1])
                    delrdeltheta = delrdeltheta / (2 * d1 * d2)

                    # calculate dh/dtheta
                    deltheta = (curvature_inputs[frm, row, col + 1] - curvature_inputs[frm, row, col - 1]) / (2 * d2)

                    # calculate d2h/dtheta2
                    del2theta = curvature_inputs[frm, row, col - 1] + curvature_inputs[frm, row, col + 1] - 2 * curvature_inputs[frm, row, col]
                    del2theta = del2theta / d2**2

                    # calculate coefficients
                    r = (row * d1) + (d1 / 2)
                    c1 = 1 / r
                    c2 = 1 / r**2
                    c3 = 1 / r**3
                    c4 = 1 / r**4
                    theta = ((col - 1) * d2) + (d2 / 2)  # col-1 because this is wrapped in theta direction

                    # calculate normal vector x,y components
                    normalization_factor = np.sqrt(1 + c2 * deltheta**2 + delr**2)
                    norm_vec_x = (c1 * np.sin(theta) * deltheta) - (np.cos(theta) * delr) / normalization_factor
                    norm_vec_y = (-1 * c1 * np.cos(theta) * deltheta) - (np.sin(theta) * delr) / normalization_factor
                    norm_vec_z = 1 / normalization_factor

                    # calculate polar laplacian and gaussian curvature
                    curvature_outputs[frm, row, col] = (del2r + c1 * delr + c2 * del2theta) / 2.0
                    kgauss_outputs[frm, row, col] = c1 * delr * del2r + 2 * c3 * delrdeltheta * deltheta - c4 * deltheta**2 - c2 * (delrdeltheta**2 - del2r * del2theta)
                    normal_vector_outputs[frm, row, col * 3] = norm_vec_x
                    normal_vector_outputs[frm, row, col * 3 + 1] = norm_vec_y
                    normal_vector_outputs[frm, row, col * 3 + 2] = norm_vec_z

                elif nan_test[frm, row, col] == False:

                    # calculate d2h/dr2
                    del2r = curvature_inputs[frm, row - 1, col] + curvature_inputs[frm, row + 1, col] - 2 * curvature_inputs[frm, row, col]
                    del2r = del2r / d1**2

                    # calculate dh/dr
                    delr = (curvature_inputs[frm, row + 1, col] - curvature_inputs[frm, row - 1, col]) / (2 * d1)

                    # calculate dh/dtheta
                    deltheta = (curvature_inputs[frm, row, col + 1] - curvature_inputs[frm, row, col - 1]) / (2 * d2)

                    # calculate d2h/dtheta2
                    del2theta = curvature_inputs[frm, row, col - 1] + curvature_inputs[frm, row, col + 1] - 2 * curvature_inputs[frm, row, col]
                    del2theta = del2theta / d2**2

                    # calculate coefficients
                    r = (row * d1) + (d1 / 2)
                    c1 = 1 / r
                    c2 = 1 / r**2
                    theta = ((col - 1) * d2) + (d2 / 2)  # col-1 because this is wrapped in theta direction

                    # calculate normal vector x,y components
                    normalization_factor = np.sqrt(1 + c2 * deltheta**2 + delr**2)
                    norm_vec_x = (c1 * np.sin(theta) * deltheta) - (np.cos(theta) * delr) / normalization_factor
                    norm_vec_y = (-1 * c1 * np.cos(theta) * deltheta) - (np.sin(theta) * delr) / normalization_factor
                    norm_vec_z = 1 / normalization_factor

                    curvature_outputs[frm, row, col] = (del2r + c1 * delr + c2 * del2theta) / 2.0
                    kgauss_outputs[frm, row, col] = np.nan
                    normal_vector_outputs[frm, row, col * 3] = norm_vec_x
                    normal_vector_outputs[frm, row, col * 3 + 1] = norm_vec_y
                    normal_vector_outputs[frm, row, col * 3 + 2] = norm_vec_z

                else:
                    curvature_outputs[frm, row, col] = np.nan
                    kgauss_outputs[frm, row, col] = np.nan
                    normal_vector_outputs[frm, row, col * 3] = np.nan
                    normal_vector_outputs[frm, row, col * 3 + 1] = np.nan
                    normal_vector_outputs[frm, row, col * 3 + 2] = np.nan

    return curvature_outputs, kgauss_outputs, normal_vector_outputs


def empty_neighbor_test(curvature_inputs):
    data = np.isnan(curvature_inputs)
    nan_test = np.array(data, copy=True)
    nan_test2 = np.array(data, copy=True)
    knan_test = np.array(data, copy=True)

    shape = np.shape(data)
    dim1 = shape[1]
    dim2 = shape[2]
    dim3 = shape[0]

    for frm in range(dim3):
        for row in range(1, dim1 - 1):
            for col in range(1, dim2 - 1):
                if nan_test2[frm, row - 1, col] == True:
                    nan_test[frm, row, col] = True
                    knan_test[frm, row, col] = True
                elif nan_test2[frm, row + 1, col] == True:
                    nan_test[frm, row, col] = True
                    knan_test[frm, row, col] = True
                elif nan_test2[frm, row, col - 1] == True:
                    nan_test[frm, row, col] = True
                    knan_test[frm, row, col] = True
                elif nan_test2[frm, row, col + 1] == True:
                    nan_test[frm, row, col] = True
                    knan_test[frm, row, col] = True
                elif nan_test2[frm, row + 1, col + 1] == True:
                    knan_test[frm, row, col] = True
                elif nan_test2[frm, row - 1, col + 1] == True:
                    knan_test[frm, row, col] = True
                elif nan_test2[frm, row + 1, col - 1] == True:
                    knan_test[frm, row, col] = True
                elif nan_test2[frm, row - 1, col - 1] == True:
                    knan_test[frm, row, col] = True

    nan_test[:, 0, :] = True
    nan_test[:, dim1 - 1, :] = True
    nan_test[:, :, 0] = True
    nan_test[:, :, dim2 - 1] = True

    knan_test[:, 0, :] = True
    knan_test[:, dim1 - 1, :] = True
    knan_test[:, :, 0] = True
    knan_test[:, :, dim2 - 1] = True

    return nan_test, knan_test
