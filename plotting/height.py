import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import warnings
import glob
import os
from utils import *
# from code_review2 import *


def Make_surface_PDB(data, name, field, d1, d2, f, serial, bead, polar):
    resseqnum = 1
    atom_name = 'SURF'
    chain = 'X'
    row, col = data.shape
    beadnum = str(bead[1])
    beadname = "C" + beadnum + "  "

    for d1bin in range(row):
        for d2bin in range(col):
            if str(data[d1bin][d2bin]) != "nan":
                seriallen = 5 - (len(str(serial)))
                resseqlen = 4 - (len(str(resseqnum)))
                if polar == 1:
                    x = (d1 * d1bin + .5 * d1) * (np.cos(d2bin * d2 + 0.5 * d2))
                    y = (d1 * d1bin + .5 * d1) * (np.sin(d2bin * d2 + 0.5 * d2))
                else:
                    L1 = d1 * row
                    L2 = d2 * col
                    x = (d1 * d1bin + .5 * d1) - L1 / 2
                    y = (d2 * d2bin + .5 * d2) - L2 / 2
                x = coord_format(x)
                y = coord_format(y)
                z = coord_format(data[d1bin][d2bin])
                d1num = bin_format(d1bin)
                d2num = bin_format(d2bin)
                print('HETATM' + (' ' * seriallen) + str(serial) + ' ' + atom_name + ' ' + beadname + chain + (' ' * resseqlen) + str(resseqnum) + '    ' + x + y + z + d1num + d2num + '      ' + field[:4] + '  ', file=f)
                serial += 1
                resseqnum += 1
    return serial


def calculate_zplus(sys_name, bead, coordsys, inclusion, polar, dims, serial, pdb, scale_dict):
    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims
    zone = np.load('npy/' + sys_name + '.zone.' + bead + '.' + coordsys + '.height.npy')
    ztwo = np.load('npy/' + sys_name + '.ztwo.' + bead + '.' + coordsys + '.height.npy')

    zplus = (zone + ztwo) / 2

    avgzplus = calc_avg_over_time(zplus)

    # make plots!
    plot_maker(dim1vals, dim2vals, avgzplus, sys_name, 'zplus', scale_dict["height_max"], scale_dict["height_min"], inclusion, "avgHeight", bead, coordsys, scale_dict["colorbar"])

    # save as file for debugging / analysis AND make PDB!
    np.save('npy/' + sys_name + '.zplus.' + bead + '.' + coordsys + '.height.npy', zplus)
    np.savetxt('dat/' + sys_name + '.zplus.' + bead + '.' + coordsys + '.avgheight.dat', avgzplus, delimiter=',', fmt='%10.5f')
    serial = Make_surface_PDB(avgzplus, sys_name, 'zplus', d1, d2, pdb, serial, bead, polar)
    print(sys_name + ' ' + bead + " zplus height done!")


def analyze_height(sys_name, names_dict, coordsys, inclusion, polar, dims, field_list, scale_dict):
    serial = 1

    N1_bins, d1, N2_bins, d2, Nframes, dim1vals, dim2vals = dims

    pdbname = sys_name + "." + coordsys + ".avgheight.pdb"

    with open(pdbname, "w") as pdb:

        # print first line of pdb file
        print('CRYST1  150.000  150.000  110.000  90.00  90.00  90.00 P 1           1', file=pdb)

        # do height analysis
        for bead in names_dict['beads_list']:
            for field in field_list:

                # import traj values
                height_data = np.genfromtxt('tcl_output/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.height.dat', missing_values='nan', filling_values=np.nan)

                # create a new array that has each frame in a different array level
                height = np.zeros((N1_bins, N2_bins, Nframes))
                for frm in range(Nframes):
                    height[:, :, frm] = height_data[frm * N1_bins:(frm + 1) * N1_bins, 2:]

                # if a bin is occupied <10% of the time, it shouldn't be treated as part of the membrane
                pruned_height = mostly_empty(height)

                # take the average height over all frames
                avgHeight = calc_avg_over_time(pruned_height)

                # make plots!
                plot_maker(dim1vals, dim2vals, avgHeight, sys_name, field, scale_dict["height_max"], scale_dict["height_min"], inclusion, "avgHeight", bead, coordsys, scale_dict["colorbar"])

                # save as file for debugging / analysis AND make PDB!
                np.save('npy/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.height.npy', pruned_height)
                np.savetxt('dat/' + sys_name + '.' + field + '.' + bead + '.' + coordsys + '.avgheight.dat', avgHeight, delimiter=',', fmt='%10.5f')
                serial = Make_surface_PDB(avgHeight, sys_name, field, d1, d2, pdb, serial, bead, polar)
                print(sys_name + ' ' + bead + ' ' + field + " height done!")

            calculate_zplus(sys_name, bead, coordsys, inclusion, polar, dims, serial, pdb, scale_dict)

        # print last line of pdb file
        print('END', file=pdb)

    return
