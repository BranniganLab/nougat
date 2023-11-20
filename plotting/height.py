"""Functions related to height analysis and pdb creation."""

import numpy as np
from pathlib import Path
from utils import *


def Make_surface_PDB(data, field, d1, d2, f, serial, polar):
    resseqnum = 1
    atom_name = 'SURF'
    chain = 'X'
    row, col = data.shape
    beadname = "C1" + "  "

    for d1bin in range(row):
        for d2bin in range(col):
            if str(data[d1bin][d2bin]) != "nan":
                seriallen = 5 - (len(str(serial)))
                resseqlen = 4 - (len(str(resseqnum)))
                if polar:
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


def calculate_zplus(polar, serial, pdb, d1, d2, cwd):
    zone = np.load(cwd.joinpath("trajectory", "height", "zone.npy"))
    ztwo = np.load(cwd.joinpath("trajectory", "height", "ztwo.npy"))

    zplus = (zone + ztwo) / 2

    avgzplus = calc_avg_over_time(zplus)

    # save as file for debugging / analysis AND make PDB!
    np.save(cwd.joinpath("trajectory", "height", "zplus.npy"), zplus)
    np.save(cwd.joinpath("average", "height", "zplus.npy"), avgzplus)
    np.savetxt(cwd.joinpath("average", "height", "zplus.dat"), avgzplus, delimiter=',', fmt='%10.5f')
    if polar:
        avg_over_theta(cwd.joinpath("average", "height", "zplus"))
    serial = Make_surface_PDB(avgzplus, 'zplus', d1, d2, pdb, serial, polar)
    print(sys_name + " zplus height done!")


def parse_height(system_dict, polar, cwd):
    serial = 1

    N1_bins = system_dict['bin_info']['N1']
    d1 = system_dict['bin_info']['d1']
    N2_bins = system_dict['bin_info']['N2']
    d2 = system_dict['bin_info']['d2']

    with open("avgheight.pdb", "w") as pdb:

        # print first line of pdb file
        print('CRYST1  150.000  150.000  110.000  90.00  90.00  90.00 P 1           1', file=pdb)

        # do height analysis

        for field in ["zone", "ztwo", "zzero"]:

            # import traj values
            height_input_path = cwd.joinpath("tcl_output", "height", field + ".dat")
            height_data = np.genfromtxt(height_input_path, missing_values='nan', filling_values=np.nan)
            Nframes = int(np.shape(height_data)[0] / N1_bins)

            if 'nframes' in system_dict['bin_info']:
                # error check: do we get the same nframes from each field?
                if Nframes != system_dict['bin_info']['nframes']:
                    raise Exception
            else:
                system_dict['bin_info']['nframes'] = Nframes

            # create a new array that has each frame in a different array level
            height = np.zeros((Nframes, N1_bins, N2_bins))
            for frm in range(Nframes):
                height[frm, :, :] = height_data[frm * N1_bins: (frm + 1) * N1_bins, 2:]

            # if a bin is occupied <10% of the time, it shouldn't be treated as part of the membrane
            pruned_height = mostly_empty(height)

            # take the average height over all frames
            avgHeight = calc_avg_over_time(pruned_height)

            # save as file for debugging / analysis AND make PDB!
            np.save(cwd.joinpath("trajectory", "height", field + ".npy"), pruned_height)
            np.save(cwd.joinpath("average", "height", field + ".npy"), avgHeight)
            if polar:
                avg_over_theta(cwd.joinpath("average", "height", field))
            np.savetxt(cwd.joinpath("average", "height", field + ".dat"), avgHeight, delimiter=',', fmt='%10.5f')
            serial = Make_surface_PDB(avgHeight, field, d1, d2, pdb, serial, polar)
            print(field + " height done!")

        calculate_zplus(polar, serial, pdb, d1, d2, cwd)

        # print last line of pdb file
        print('END', file=pdb)

    return system_dict
