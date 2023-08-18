import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os
from utils import *
# from code_review2 import *


def init_curvature_data(height, polar, dims):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    # create arrays for storing curvature data
    if polar is True:
        curvature_inputs = np.zeros((N1_bins, N2_bins + 2, Nframes))
        curvature_outputs = np.zeros((N1_bins, N2_bins + 2, Nframes))
        kgauss_outputs = np.zeros((N1_bins, N2_bins + 2, Nframes))
        normal_vector_outputs = np.zeros((N1_bins, 3 * (N2_bins + 2), Nframes))
    elif polar is False:
        curvature_inputs = np.zeros((N1_bins + 2, N2_bins + 2, Nframes))
        curvature_outputs = np.zeros((N1_bins + 2, N2_bins + 2, Nframes))
        kgauss_outputs = np.zeros((N1_bins + 2, N2_bins + 2, Nframes))
        normal_vector_outputs = np.zeros((N1_bins + 2, 3 * (N2_bins + 2), Nframes))

    if polar is True:
        # wrap the inputs in the theta direction for calculating curvature
        curvature_inputs[:, 1:(N2_bins + 1), :] = height
        curvature_inputs[:, 0, :] = curvature_inputs[:, N2_bins, :]
        curvature_inputs[:, (N2_bins + 1), :] = curvature_inputs[:, 1, :]
    elif polar is False:
        # if cartesian, wrap in both directions
        curvature_inputs[1:(N1_bins + 1), 1:(N2_bins + 1), :] = height
        curvature_inputs[:, 0, :] = curvature_inputs[:, N2_bins, :]
        curvature_inputs[:, (N2_bins + 1), :] = curvature_inputs[:, 1, :]
        curvature_inputs[0, :, :] = curvature_inputs[N1_bins, :, :]
        curvature_inputs[(N1_bins + 1), :, :] = curvature_inputs[1, :, :]
        # and fill in the corners
        curvature_inputs[0, 0, :] = curvature_inputs[N1_bins, N2_bins, :]
        curvature_inputs[N1_bins + 1, N2_bins + 1, :] = curvature_inputs[1, 1, :]
        curvature_inputs[0, N2_bins + 1, :] = curvature_inputs[N1_bins, 1, :]
        curvature_inputs[N1_bins + 1, 0, :] = curvature_inputs[1, N2_bins, :]

    return curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs


def calculate_curvature(sys_name, bead, coordsys, inclusion, polar, dims, field_list, scale_dict):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    leaflist = field_list.copy()
    leaflist.append("zplus")

    for field in leaflist:
        field_height = np.load('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.height.npy')

        curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs = init_curvature_data(field_height, polar, dims)

        # if a bin is empty, you can't (nicely) measure the curvature of its neighbors
        nan_test, knan_test = empty_neighbor_test(curvature_inputs)

        # measure the laplacian and gaussian curvatures
        if polar is True:
            curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_polar(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims)
            # curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_polar(dims, curvature_inputs)
        elif polar is False:
            curvature_outputs, kgauss_outputs, normal_vector_outputs = measure_curvature_cart(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims)

        # unwrap along dim2 direction
        meancurvature = curvature_outputs[:, 1:N2_bins + 1, :]
        kcurvature = kgauss_outputs[:, 1:N2_bins + 1, :]
        normal_vectors = normal_vector_outputs[:, 3:3 * (N2_bins + 1), :]

        # if cartesian, unwrap along dim1 direction too
        if polar is False:
            meancurvature = meancurvature[1:N1_bins + 1, :, :]
            kcurvature = kcurvature[1:N1_bins + 1, :, :]
            normal_vectors = normal_vectors[1:N1_bins + 1, :, :]

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


def measure_curvature_cart(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    # mean curvature: Hxx + Hyy
    # gaussian curvature: HxxHyy - Hxy^2

    for frm in range(Nframes):
        for row in range(N1_bins + 2):
            for col in range(N2_bins + 2):
                if knan_test[row, col, frm] == False:

                    del2x = curvature_inputs[row - 1, col, frm] + curvature_inputs[row + 1, col, frm] - 2 * curvature_inputs[row, col, frm]
                    del2x = del2x / (d1**2)

                    del2y = curvature_inputs[row, col - 1, frm] + curvature_inputs[row, col + 1, frm] - 2 * curvature_inputs[row, col, frm]
                    del2y = del2y / (d2**2)

                    delxy = (curvature_inputs[row + 1, col + 1, frm] - curvature_inputs[row + 1, col, frm] - curvature_inputs[row, col + 1, frm] + 2 * curvature_inputs[row, col, frm] - curvature_inputs[row - 1, col, frm] - curvature_inputs[row, col - 1, frm] + curvature_inputs[row - 1, col - 1, frm])
                    delxy = delxy / (2 * d1 * d2)

                    # delxy = curvature_inputs[row+1,col+1,frm] - curvature_inputs[row+1,col-1,frm] - curvature_inputs[row-1,col+1,frm] + curvature_inputs[row-1,col-1,frm]
                    # delxy = delxy / (4*d1*d2)

                    delx = (curvature_inputs[row + 1, col, frm] - curvature_inputs[row - 1, col, frm]) / (2 * d1)

                    dely = (curvature_inputs[row, col + 1, frm] - curvature_inputs[row, col - 1, frm]) / (2 * d2)

                    normalization_factor = np.sqrt(1 + delx**2 + dely**2)
                    norm_vec_x = -1 * delx / normalization_factor
                    norm_vec_y = -1 * dely / normalization_factor
                    norm_vec_z = 1 / normalization_factor

                    curvature_outputs[row, col, frm] = (del2x + del2y) / 2.0
                    kgauss_outputs[row, col, frm] = del2x * del2y - delxy**2
                    normal_vector_outputs[row, col * 3, frm] = norm_vec_x
                    normal_vector_outputs[row, col * 3 + 1, frm] = norm_vec_y
                    normal_vector_outputs[row, col * 3 + 2, frm] = norm_vec_z

                elif nan_test[row, col, frm] == False:
                    del2x = curvature_inputs[row - 1, col, frm] + curvature_inputs[row + 1, col, frm] - 2 * curvature_inputs[row, col, frm]
                    del2x = del2x / d1**2

                    del2y = curvature_inputs[row, col - 1, frm] + curvature_inputs[row, col + 1, frm] - 2 * curvature_inputs[row, col, frm]
                    del2y = del2y / d2**2

                    delx = (curvature_inputs[row + 1, col, frm] - curvature_inputs[row - 1, col, frm]) / (2 * d1)

                    dely = (curvature_inputs[row, col + 1, frm] - curvature_inputs[row, col - 1, frm]) / (2 * d2)

                    normalization_factor = np.sqrt(1 + delx**2 + dely**2)
                    norm_vec_x = -1 * delx / normalization_factor
                    norm_vec_y = -1 * dely / normalization_factor
                    norm_vec_z = 1 / normalization_factor

                    curvature_outputs[row, col, frm] = (del2x + del2y) / 2.0
                    kgauss_outputs[row, col, frm] = np.nan
                    normal_vector_outputs[row, col * 3, frm] = norm_vec_x
                    normal_vector_outputs[row, col * 3 + 1, frm] = norm_vec_y
                    normal_vector_outputs[row, col * 3 + 2, frm] = norm_vec_z

                else:

                    curvature_outputs[row, col, frm] = np.nan
                    kgauss_outputs[row, col, frm] = np.nan
                    normal_vector_outputs[row, col * 3, frm] = np.nan
                    normal_vector_outputs[row, col * 3 + 1, frm] = np.nan
                    normal_vector_outputs[row, col * 3 + 2, frm] = np.nan

    return curvature_outputs, kgauss_outputs, normal_vector_outputs


def measure_curvature_polar(curvature_inputs, curvature_outputs, kgauss_outputs, normal_vector_outputs, nan_test, knan_test, dims):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    # mean curvature: 1/2 * [h_rr + 1/r(h_r) + 1/r**2(h_thetatheta)]
    # gaussian curvature: 1/r(h_r*h_rr) + 2/r**3(h_rtheta*h_theta) - 1/r**4(h_theta**2) - 1/r**2(h_rtheta**2 - h_rr*h_thetatheta)

    for frm in range(Nframes):
        for row in range(N1_bins):
            for col in range(N2_bins + 2):
                if knan_test[row, col, frm] == False:

                    # calculate d2h/dr2
                    del2r = curvature_inputs[row - 1, col, frm] + curvature_inputs[row + 1, col, frm] - 2 * curvature_inputs[row, col, frm]
                    del2r = del2r / d1**2

                    # calculate dh/dr
                    delr = (curvature_inputs[row + 1, col, frm] - curvature_inputs[row - 1, col, frm]) / (2 * d1)

                    # calculate d2h/drdtheta
                    delrdeltheta = (curvature_inputs[row + 1, col + 1, frm] - curvature_inputs[row + 1, col, frm] - curvature_inputs[row, col + 1, frm] + 2 * curvature_inputs[row, col, frm] - curvature_inputs[row - 1, col, frm] - curvature_inputs[row, col - 1, frm] + curvature_inputs[row - 1, col - 1, frm])
                    delrdeltheta = delrdeltheta / (2 * d1 * d2)

                    # calculate dh/dtheta
                    deltheta = (curvature_inputs[row, col + 1, frm] - curvature_inputs[row, col - 1, frm]) / (2 * d2)

                    # calculate d2h/dtheta2
                    del2theta = curvature_inputs[row, col - 1, frm] + curvature_inputs[row, col + 1, frm] - 2 * curvature_inputs[row, col, frm]
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
                    curvature_outputs[row, col, frm] = (del2r + c1 * delr + c2 * del2theta) / 2.0
                    kgauss_outputs[row, col, frm] = c1 * delr * del2r + 2 * c3 * delrdeltheta * deltheta - c4 * deltheta**2 - c2 * (delrdeltheta**2 - del2r * del2theta)
                    normal_vector_outputs[row, col * 3, frm] = norm_vec_x
                    normal_vector_outputs[row, col * 3 + 1, frm] = norm_vec_y
                    normal_vector_outputs[row, col * 3 + 2, frm] = norm_vec_z

                elif nan_test[row, col, frm] == False:

                    # calculate d2h/dr2
                    del2r = curvature_inputs[row - 1, col, frm] + curvature_inputs[row + 1, col, frm] - 2 * curvature_inputs[row, col, frm]
                    del2r = del2r / d1**2

                    # calculate dh/dr
                    delr = (curvature_inputs[row + 1, col, frm] - curvature_inputs[row - 1, col, frm]) / (2 * d1)

                    # calculate dh/dtheta
                    deltheta = (curvature_inputs[row, col + 1, frm] - curvature_inputs[row, col - 1, frm]) / (2 * d2)

                    # calculate d2h/dtheta2
                    del2theta = curvature_inputs[row, col - 1, frm] + curvature_inputs[row, col + 1, frm] - 2 * curvature_inputs[row, col, frm]
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

                    curvature_outputs[row, col, frm] = (del2r + c1 * delr + c2 * del2theta) / 2.0
                    kgauss_outputs[row, col, frm] = np.nan
                    normal_vector_outputs[row, col * 3, frm] = norm_vec_x
                    normal_vector_outputs[row, col * 3 + 1, frm] = norm_vec_y
                    normal_vector_outputs[row, col * 3 + 2, frm] = norm_vec_z

                else:
                    curvature_outputs[row, col, frm] = np.nan
                    kgauss_outputs[row, col, frm] = np.nan
                    normal_vector_outputs[row, col * 3, frm] = np.nan
                    normal_vector_outputs[row, col * 3 + 1, frm] = np.nan
                    normal_vector_outputs[row, col * 3 + 2, frm] = np.nan

    return curvature_outputs, kgauss_outputs, normal_vector_outputs


def empty_neighbor_test(curvature_inputs):
    data = np.isnan(curvature_inputs)
    nan_test = np.array(data, copy=True)
    nan_test2 = np.array(data, copy=True)
    knan_test = np.array(data, copy=True)

    shape = np.shape(data)
    dim1 = shape[0]
    dim2 = shape[1]
    dim3 = shape[2]

    for frm in range(dim3):
        for row in range(1, dim1 - 1):
            for col in range(1, dim2 - 1):
                if nan_test2[row - 1, col, frm] == True:
                    nan_test[row, col, frm] = True
                    knan_test[row, col, frm] = True
                elif nan_test2[row + 1, col, frm] == True:
                    nan_test[row, col, frm] = True
                    knan_test[row, col, frm] = True
                elif nan_test2[row, col - 1, frm] == True:
                    nan_test[row, col, frm] = True
                    knan_test[row, col, frm] = True
                elif nan_test2[row, col + 1, frm] == True:
                    nan_test[row, col, frm] = True
                    knan_test[row, col, frm] = True
                elif nan_test2[row + 1, col + 1, frm] == True:
                    knan_test[row, col, frm] = True
                elif nan_test2[row - 1, col + 1, frm] == True:
                    knan_test[row, col, frm] = True
                elif nan_test2[row + 1, col - 1, frm] == True:
                    knan_test[row, col, frm] = True
                elif nan_test2[row - 1, col - 1, frm] == True:
                    knan_test[row, col, frm] = True

    nan_test[0, :, :] = True
    nan_test[dim1 - 1, :, :] = True
    nan_test[:, 0, :] = True
    nan_test[:, dim2 - 1, :] = True

    knan_test[0, :, :] = True
    knan_test[dim1 - 1, :, :] = True
    knan_test[:, 0, :] = True
    knan_test[:, dim2 - 1, :] = True

    return nan_test, knan_test
