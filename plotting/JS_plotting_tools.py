"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
import os
from pathlib import Path
from utils import *

sys_dict = {
    "DT": "green",
    "lgDT": "green",
    "DL": "blue",
    "lgDL": "blue",
    "DX": "purple",
    "lgDX": "purple",
    "DB": "red",
    "lgDB": "red",
    "DY": "limegreen",
    "lgDY": "limegreen",
    "DO": "mistyrose",
    "lgDO": "mistyrose",
    "PO": "darkorange",
    "lgPO": "darkorange",
    "DP": "deepskyblue",
    "lgDP": "deepskyblue",
    "DG": "orchid",
    "lgDG": "orchid"
}

box_dict = {
    "large": "solid",
    "small": "dotted"
}

field_dict = {
    "zone": "solid",
    "ztwo": "solid",
    "zzero": "dotted",
    "zplus": "dashed"
}

bead_dict = {
    "DT": ['C2A.C2B'],
    "DL": ['C2A.C2B', 'C3A.C3B'],
    "DY": ['D2A.D2B', 'C3A.C3B'],
    "DO": ['D2A.D2B', 'C3A.C3B', 'C4A.C4B'],
    "PO": ['D2A.C2B', 'C3A.C3B', 'C4A.C4B'],
    "DP": ['C2A.C2B', 'C3A.C3B', 'C4A.C4B'],
    "DB": ['C2A.C2B', 'C3A.C3B', 'C4A.C4B', 'C5A.C5B'],
    "DG": ['C2A.C2B', 'D3A.D3B', 'C4A.C4B', 'C5A.C5B'],
    "DX": ['C2A.C2B', 'C3A.C3B', 'C4A.C4B', 'C5A.C5B', 'C6A.C6B']
}

small_sys_list = ["DT", "DY", "DL", "DO", "DP", "PO", "DG", "DB", "DX"]
large_sys_list = ["lgDT", "lgDY", "lgDL", "lgDO", "lgDP", "lgPO", "lgDG", "lgDB", "lgDX"]


def zoom_in(systems):
    max_scale_dict = {
        "height": 60,
        "thickness": 2,
        "curvature": 0.01,
        "Kcurvature": 0.001,
        "tail1": 0.6,
        "tail0": 0.6,
        "density": 2
    }
    min_scale_dict = {
        "height": -60,
        "thickness": 0,
        "curvature": -0.01,
        "Kcurvature": -0.001,
        "tail1": 0.0,
        "tail0": 0.0,
        "density": 0
    }
    values = ["height", "thickness", "curvature", "Kcurvature", "tail1", "tail0"]
    # values = ["height"]
    for system in systems:
        # os.chdir("dm1/"+system)
        os.chdir(system)
        os.chdir(system + '_polar_5_10_100_-1_1/')
        for field in ["zone", "ztwo"]:
            for value in values:
                print(value)
                if value == "tail0" or value == "tail1":
                    data = np.genfromtxt('dat/' + system + '.' + system[2:] + 'PC.' + value + '.' + field + '.polar.avgOrder.dat', delimiter=",", missing_values="nan", filling_values=np.nan)
                else:
                    data = np.genfromtxt('dat/' + system + '.' + field + '.C1A.C1B.polar.avg' + value + '.dat', delimiter=",", missing_values="nan", filling_values=np.nan)
                dim1 = np.linspace(0, 50, 16)
                dim2 = np.linspace(0, 2 * np.pi, 11)
                dim1vals, dim2vals = np.meshgrid(dim1, dim2, indexing='ij')
                print(data[:15, :].shape)
                plot_maker(dim1vals, dim2vals, data[:15, :], system, field, max_scale_dict[value], min_scale_dict[value], False, value, False, "polar", True)
        os.chdir('../..')


def sum_terms(systems, mol):
    for system in systems:
        filename_start = '/home/js2746/Bending/PC/whole_mols/' + mol + '/' + system + '/npy/' + system + '.'
        filename_end = '.C1A.C1B.cart.height.npy'
        try:
            zzero = np.load(filename_start + 'zzero' + filename_end)
            H1 = np.load(filename_start + 'zone' + '.C1A.C1B.cart.meancurvature.npy')
            H2 = np.load(filename_start + 'ztwo' + '.C1A.C1B.cart.meancurvature.npy')
        except:
            filename_end = '.C1A.C1B.D1A.D1B.cart.height.npy'
            zzero = np.load(filename_start + 'zzero' + filename_end)
            H1 = np.load(filename_start + 'zone' + '.C1A.C1B.D1A.D1B.cart.meancurvature.npy')
            H2 = np.load(filename_start + 'ztwo' + '.C1A.C1B.D1A.D1B.cart.meancurvature.npy')
        zplus = np.load(filename_start + 'zplus' + filename_end)
        epsilon = zplus - zzero
        epsilon2 = epsilon * epsilon
        H = H1 + H2
        H2 = H * H
        avgH2 = np.nanmean(H2, axis=2)
        sumavgH2 = np.nansum(avgH2)
        t0 = measure_t0(system, mol)
        print(t0)
        term6 = 2 * H * epsilon / t0
        term6 = np.nanmean(term6, axis=2)
        term7 = epsilon2 / (2 * t0**2)
        term7 = np.nanmean(term7, axis=2)
        sumterm6 = np.nansum(term6)
        sumterm7 = np.nansum(term7)
        print(sumterm6)
        print(sumterm7)
        print(sumavgH2)


def sum_over_K(systems):

    for system in systems:
        fig, ax = plt.subplots()
        for field in ["zone", "ztwo"]:
            data = np.load(system + '/npy/' + system + '.' + field + '.C1A.C1B.cart.gausscurvature.npy')
            k_list = []
            for frame in range(np.shape(data)[2]):
                frame_data = data[:, :, frame]
                k_sum = np.nansum(frame_data)
                k_list.append(k_sum)
            avg_list = rollingavg(k_list, 50)
        # if "lg" in system:
        # avg_list = avg_list/(40**2)
            # else:
            # avg_list = avg_list/(12**2)
            if field == "zone":
                style = "solid"
            else:
                style = "dashed"
            plt.plot(avg_list, linestyle=style)
            plt.ylim([0.0, .01])
        plt.savefig(system + ".sum_over_K.pdf")


def rollingavg(list_in, window_size):
    rollingavg = []
    for i in range(window_size):
        rollingavg.append(np.nan)
    for i in range(window_size, len(list_in) - window_size):
        windowavg = np.sum(list_in[i - window_size:i + window_size]) / (2 * window_size)
        rollingavg.append(windowavg)
    for i in range(window_size):
        rollingavg.append(np.nan)
    rollingavg = np.array(rollingavg)
    return rollingavg


def sum_over_H2(systems, system_names, groupname, nougvals, mol):
    sum_list = []
    for system, name in zip(systems, system_names):
        data = np.load(name + "/" + name + "_polar_" + nougvals + "/npy/" + name + '.avg_H_plus2.npy')
        areas = np.load(name + "/" + name + "_polar_" + nougvals + "/npy/" + name + '.' + coordsys + '.areas.npy')
        normalized = data / areas
        sumHplus2 = np.nansum(normalized)
        sum_list.append(sumHplus2)
    fig = plt.figure()
    plt.ylim(0, .01)
    plt.bar(systems, sum_list)
    plt.show()


def make_2d_series_over_time(path, quantity, polar, sys_name):
    """
    Make movie file of heatmap over trajectory.

    Parameters
    ----------
    path : string
        the path the to nougat outputs directory you want.
    quantity : string
        the name of the measurement.
    polar : bool
        Whether or not to use polar coordinates.
    sys_name : string
        the name you gave nougat.py when it made your files.

    Returns
    -------
    None.

    """
    cwd = os.getcwd()
    if cwd != path:
        os.chdir(path)

    if polar:
        coordsys = "polar"
    else:
        coordsys = "cart"

    # load the correct 2d data
    traj_data = np.load(path + "/npy/" + sys_name + "." + quantity + ".npy")
    nframes = np.shape(traj_data)[2]
    dims = bin_prep(sys_name, "C1A.C1B", polar, 'OFF')
    dim1vals = dims[5]
    dim2vals = dims[6]

    # for each frame, make plot
    image_dir = os.path.join(path, quantity + "_over_time")
    try:
        os.mkdir(image_dir)
    except IOError:
        pass
    os.chdir(image_dir)
    try:
        os.mkdir("pdf")
    except IOError:
        pass
    for frame in range(nframes):
        frame_data = traj_data[:, :, frame]
        plot_maker(dim1vals, dim2vals, frame_data, sys_name, str(frame), "auto", "auto", False, quantity, False, coordsys, True)

    os.chdir(cwd)


def run_sum_over_H2(mol):
    if mol == "5x29":
        sys_names = []
        for system in satsys:
            sys_names.append("lg" + system)
        sum_over_H2(satsys, sys_names, "satsys", "5_10_100_-1_1", "5x29")
        sys_names = []
        for system in monounsatsys:
            sys_names.append("lg" + system)
        sum_over_H2(monounsatsys, sys_names, "monounsatsys", "5_10_100_-1_1", "5x29")
    elif mol == "7k3g":
        sys_names = []
        for system in satsys:
            sys_names.append(system + "PC")
        sum_over_H2(satsys, sys_names, "satsys", "5_10_0_-1_1", "7k3g")
        sys_names = []
        for system in monounsatsys:
            sys_names.append(system + "PC")
        sum_over_H2(monounsatsys, sys_names, "monounsatsys", "5_10_0_-1_1", "7k3g")


def run_eps_corr_scatter(mol):
    if mol == "5x29":
        sys_names = []
        for system in satsys:
            sys_names.append("lg" + system)
        plot_eps_corr_scatter(satsys, sys_names, "satsys", "5_10_100_-1_1", "5x29")
        sys_names = []
        for system in monounsatsys:
            sys_names.append("lg" + system)
        plot_eps_corr_scatter(monounsatsys, sys_names, "monounsatsys", "5_10_100_-1_1", "5x29")
    elif mol == "7k3g":
        sys_names = []
        for system in satsys:
            sys_names.append(system + "PC")
        plot_eps_corr_scatter(satsys, sys_names, "satsys", "5_10_0_-1_1", "7k3g")
        sys_names = []
        for system in monounsatsys:
            sys_names.append(system + "PC")
        plot_eps_corr_scatter(monounsatsys, sys_names, "monounsatsys", "5_10_0_-1_1", "7k3g")


def plot_eps_corr_scatter(systems, system_names, groupname, nougvals, mol):
    for system, name in zip(systems, system_names):
        curvlist = []
        epslist = []
        curv = np.load(name + "/" + name + "_polar_" + nougvals + "/npy/" + name + '.H_plus.npy')
        eps = np.load(name + "/" + name + "_polar_" + nougvals + "/npy/" + name + '.epsilon.npy')
        for i in range(len(curv[:, 0, 0])):
            for j in range(len(curv[0, :, 0])):
                for k in range(len(curv[0, 0, :])):
                    if (np.isnan(curv[i, j, k])) or (np.isnan(eps[i, j, k])):
                        continue
                    else:
                        curvlist.append(curv[i, j, k])
                        epslist.append(eps[i, j, k])
        print(system, np.corrcoef(curvlist, epslist))
        # fig = plt.figure()
        # plt.scatter(curvlist,epslist)
        # m, b = np.polyfit(curvlist, epslist, 1)
        # axes = plt.gca()
        # x_vals = np.array(axes.get_xlim())
        # y_vals = b + m * x_vals
        # plt.plot(x_vals, y_vals, '-',color="green")
        # plt.show()


def plot_together(mols, paths, nougvals, xlim):
    mainpath = "/home/js2746/Bending/PC/whole_mols"
    dirname = os.path.join(mainpath, "avg_over_theta")
    try:
        os.mkdir(dirname)
    except:
        pass

    for quantity in ["avg_K_plus", "avg_K_minus", "corr_epst0_Kplus", "corr_mag_epst0_Hplus", "corr_epst0_Hplus", "avg_epsilon", "avg_rms_epsilon_over_t0", "avg_abs_epsilon", "avg_abs_epsilon_over_t0", "avg_epsilon2", "avg_H_plus", "avg_H_plus2", "avg_H_minus", "avg_H_minus2", "avg_epsilon_H", "avg_total_t", "avg_z_minus", "avg_z_minus2", "avg_z_minus_H_minus", "avg_epsilon_over_t0", "avg_epsilon_H_over_t0", "avg_epsilon2_over_t02", "avg_tilde_t", "avg_z_minus2_over_t02", "avg_z_minus_H_minus_over_t0"]:
        for series, seriesname in zip([satsys, monounsatsys], ["satsys", "monounsatsys"]):
            fig, axs = plt.subplots(2, sharex=True, sharey=True, layout="constrained")
            axs[0].set_xlim(0, xlim)
            limname = "zoomed_" + str(xlim)
            if quantity in min_scale_dict and xlim == 6:
                axs[0].set_ylim(min_scale_dict[quantity], max_scale_dict[quantity])
            for name in series:
                for mol, path, nougval in zip(mols, paths, nougvals):
                    if mol == "5x29":
                        sysname = "lg" + name
                        figax = 1
                    elif mol == "7k3g":
                        sysname = name + "PC"
                        figax = 0
                    npzfile = np.load(path + "/" + sysname + "/" + sysname + nougval + "/npy/" + sysname + "." + quantity + ".avg_over_theta.npz")
                    if mol == "5x29":
                        start = np.where(npzfile['x'] == 2.75)[0][0]
                    elif mol == "7k3g":
                        start = np.where(npzfile['x'] == 1.25)[0][0]
                    axs[figax].plot(npzfile['x'][start:], npzfile['y'][start:], color=colordict[name])
            # axs[0].axvline(x=1.25, color = 'black', linestyle='--')
            # axs[1].axvline(x=1.25, color = 'black', linestyle='--')
            axs[0].axvline(x=2.75, color='black', linestyle=':')
            axs[0].set_title("7k3g")
            axs[1].set_title("5x29")
            fig.supxlabel(r'$r \;(\mathrm{nm})$')
            fig.supylabel(legend_dict[quantity])
            plt.savefig(mainpath + "/avg_over_theta/" + seriesname + "_" + quantity + "_" + limname + ".pdf", dpi=700)
            plt.clf()
            plt.close()


def make_paper_writing_group_plot(comparison):
    """
    Make plot for paper writing group.

    Parameters
    ----------
    saturation : STRING
        "sat" or "unsat.

    Returns
    -------
    None.

    """
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex=True, figsize=(7, 7))
    if comparison == "sat":
        sys_list = ["lgDT", "lgDL", "lgDP", "lgDB", "lgDX"]
        nougat_values = "_polar_5_10_100_-1_1/npy/"
        path = "/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/"
    elif comparison == "unsat":
        sys_list = ["lgDY", "lgDO", "lgDG"]
        nougat_values = "_polar_10_10_100_150_1/npy/"
        path = "/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/"
    elif comparison == "stiffness":
        sys_list = ["100kjmol", "1000kjmol", "5000", "lgPO_15"]
        nougat_values = "_polar_10_10_0_-1_1/npy/"
        path = "/home/js2746/5x29_stiffness/"
    elif comparison == "elastic":
        sys_list = ["COMtiltspin", "new_gmx_pos"]
        nougat_values = "_polar_10_10_0_-1_2/npy/"
        path = "/home/js2746/Bending/PC/whole_mols/5x29/"
    elif comparison == "sat_elas":
        sys_list = ["DTPC", "DLPC", "DPPC", "DBPC", "DXPC"]
        nougat_values = "_polar_10_30_0_-1_1"
        path = Path("/home/js2746/local_deformations/")
    elif comparison == "unsat_elas":
        sys_list = ["DYPC", "DOPC", "DGPC"]
        nougat_values = "_polar_10_30_0_-1_1"
        path = Path("/home/js2746/local_deformations/")
    for system in sys_list:
        if system == "lgDY":
            dirname = "lgDY_15us"
            sysname = "lgDY_15us"
        elif system == "lgDO":
            dirname = "lgDO_30us"
            sysname = "lgDO_30us"
        elif system == "lgDG":
            dirname = "lgDG_15us"
            sysname = "lgDG_15"
        elif system == "lgPO_15" and comparison == "stiffness":
            dirname = "lgPO_50us"
            sysname = "lgPO_15"
            nougat_values = "_polar_10_10_100_150_1/npy/"
            path = "/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/"
        elif system in ["lgDT", "lgDL", "lgDP", "lgDB", "lgDX"]:
            dirname = system + "_2us"
            sysname = system
        elif system == "100kjmol":
            dirname = "100"
            sysname = system
        elif system == "1000kjmol":
            dirname = "1000"
            sysname = system
        else:
            dirname = system
            sysname = system
        fname = path.joinpath(sysname + nougat_values, "average")
        zone = np.load(fname.joinpath("height", "zone.avg_over_theta.npy"))
        zonestd = np.load(fname.joinpath("height", "zone.avg_over_theta.std.npy")) / np.sqrt(30)
        ztwo = np.load(fname.joinpath("height", "ztwo.avg_over_theta.npy"))
        ztwostd = np.load(fname.joinpath("height", "ztwo.avg_over_theta.std.npy")) / np.sqrt(30)
        tilde_t = np.load(fname.joinpath("misc", "avg_tilde_total_t.avg_over_theta.npy"))
        tilde_tstd = np.load(fname.joinpath("misc", "avg_tilde_total_t.avg_over_theta.std.npy")) / np.sqrt(30)
        epst0 = np.load(fname.joinpath("misc", "avg_epsilon_over_t0.avg_over_theta.npy"))
        epst0std = np.load(fname.joinpath("misc", "avg_epsilon_over_t0.avg_over_theta.std.npy")) / np.sqrt(30)
        tilde_eps2 = np.load(fname.joinpath("misc", "avg_rms_tilde_epsilon2.avg_over_theta.npy"))
        tilde_eps2std = np.load(fname.joinpath("misc", "avg_rms_tilde_epsilon2.avg_over_theta.std.npy")) / np.sqrt(30)
        tilde_Hplus2 = np.load(fname.joinpath("misc", "avg_rms_tilde_H_plus2.avg_over_theta.npy"))
        tilde_Hplus2std = np.load(fname.joinpath("misc", "avg_rms_tilde_H_plus2.avg_over_theta.std.npy")) / np.sqrt(30)
        corr_Hplus_eps = np.load(fname.joinpath("misc", "corr_eps_Hplus.avg_over_theta.npy"))
        corr_Hplus_epsstd = np.load(fname.joinpath("misc", "corr_eps_Hplus.avg_over_theta.std.npy")) / np.sqrt(30)
        epsHplus = np.load(fname.joinpath("misc", "avg_epsilon_H.avg_over_theta.npy"))
        epsHplusstd = np.load(fname.joinpath("misc", "avg_epsilon_H.avg_over_theta.npy")) / np.sqrt(30)

        # figure out what the x axis values should be
        tcl_output = np.genfromtxt(path.joinpath(sysname + nougat_values, "tcl_output", "height", "zone.dat"),
                                   missing_values='nan', filling_values=np.nan)
        Nr, dr, _, _, _, _ = dimensions_analyzer(tcl_output, True)
        rmin = 1.25
        xmin = dr / 2
        xmax = Nr * dr - xmin
        x_vals = np.linspace(xmin, xmax, Nr) / 10
        flag = True
        i = 0
        a = .4
        if rmin > x_vals[i]:
            while flag is True and i < len(x_vals):
                i += 1
                if rmin <= x_vals[i]:
                    flag = False

        x_vals = x_vals[i:]

        # height plot
        zone = zone[i:] / 10
        ztwo = ztwo[i:] / 10
        zonestd = zonestd[i:] / 10
        ztwostd = ztwostd[i:] / 10
        ax1.plot(x_vals, zone, color=colordict[system], label=mismatch_dict[system])
        ax1.plot(x_vals, ztwo, color=colordict[system], linestyle="dashed")
        ax1.fill_between(x_vals, (zone - zonestd), (zone + zonestd), alpha=a, color=colordict[system])
        ax1.fill_between(x_vals, (ztwo - ztwostd), (ztwo + ztwostd), alpha=a, color=colordict[system], linestyle="dashed")
        ax1.set_ylabel(r'$\langle \Delta z \rangle \; (\mathrm{nm})$', fontsize=10)
        ax1.tick_params(axis='both', which='major', labelsize=7)
        ax1.tick_params(axis='both', which='minor', labelsize=7)
        ax1.yaxis.set_major_formatter('{x:5.3f}')
        ax1.text(0.02, 0.95, "A", transform=ax1.transAxes, fontsize=10, va='top')
        # ax1.set_ylim(-5, 1.5)

        # tilde t plot
        tilde_t = tilde_t[i:]
        tilde_tstd = tilde_tstd[i:]
        ax2.plot(x_vals, tilde_t, color=colordict[system])
        ax2.fill_between(x_vals, (tilde_t - tilde_tstd), (tilde_t + tilde_tstd), alpha=a, color=colordict[system])
        ax2.set_ylabel(legend_dict['avg_tilde_total_t'], fontsize=10)
        ax2.axhline(1, color="gray", linestyle="--")
        ax2.tick_params(axis='both', which='major', labelsize=7)
        ax2.tick_params(axis='both', which='minor', labelsize=7)
        ax2.yaxis.set_major_formatter('{x:5.3f}')
        ax2.text(0.02, 0.95, "B", transform=ax2.transAxes, fontsize=10, va='top')
        # ax2.set_ylim(.8, 1.25)

        # eps/t0 plot
        epst0 = epst0[i:]
        epst0std = epst0std[i:]
        ax3.plot(x_vals, epst0, color=colordict[system])
        ax3.fill_between(x_vals, (epst0 - epst0std), (epst0 + epst0std), alpha=a, color=colordict[system])
        ax3.set_ylabel(legend_dict['avg_epsilon_over_t0'], fontsize=10)
        ax3.tick_params(axis='both', which='major', labelsize=7)
        ax3.tick_params(axis='both', which='minor', labelsize=7)
        ax3.axhline(0, color="gray", linestyle="--")
        ax3.yaxis.set_major_formatter('{x:5.3f}')
        ax3.text(0.02, 0.95, "C", transform=ax3.transAxes, fontsize=10, va='top')
        # ax3.set_ylim(-.1, .015)

        # tilde epsilon squared plot
        tilde_eps2 = tilde_eps2[i:]
        tilde_eps2std = tilde_eps2std[i:]
        ax4.plot(x_vals, tilde_eps2, color=colordict[system])
        ax4.fill_between(x_vals, (tilde_eps2 - tilde_eps2std), (tilde_eps2 + tilde_eps2std), alpha=a, color=colordict[system])
        ax4.set_ylabel(legend_dict['avg_rms_tilde_epsilon2'], fontsize=10)
        ax4.axhline(1, color="gray", linestyle="--")
        ax4.tick_params(axis='both', which='major', labelsize=7)
        ax4.tick_params(axis='both', which='minor', labelsize=7)
        ax4.yaxis.set_major_formatter('{x:5.3f}')
        ax4.text(0.02, 0.95, "D", transform=ax4.transAxes, fontsize=10, va='top')
        # ax4.set_ylim(0, 18)

        # tilde Hplus squared plot
        tilde_Hplus2 = tilde_Hplus2[i:]
        tilde_Hplus2std = tilde_Hplus2std[i:]
        ax5.plot(x_vals, tilde_Hplus2, color=colordict[system])
        ax5.fill_between(x_vals, (tilde_Hplus2 - tilde_Hplus2std), (tilde_Hplus2 + tilde_Hplus2std), alpha=a, color=colordict[system])
        ax5.set_ylabel(legend_dict['avg_rms_tilde_H_plus2'], fontsize=10)
        ax5.axhline(1, color="gray", linestyle="--")
        ax5.tick_params(axis='both', which='major', labelsize=7)
        ax5.tick_params(axis='both', which='minor', labelsize=7)
        ax5.yaxis.set_major_formatter('{x:5.3f}')
        ax5.text(0.02, 0.95, "E", transform=ax5.transAxes, fontsize=10, va='top')
        # ax5.set_ylim(0, 22)

        """
        # correlation between epsilon and Hplus plot
        corr_Hplus_eps = corr_Hplus_eps[i:]
        corr_Hplus_epsstd = corr_Hplus_epsstd[i:]
        ax6.plot(x_vals, corr_Hplus_eps, color=colordict[system])
        ax6.fill_between(x_vals, (corr_Hplus_eps - corr_Hplus_epsstd), (corr_Hplus_eps + corr_Hplus_epsstd), alpha=a, color=colordict[system])
        # ax6.set_ylim(-.25, 0)
        ax6.axhline(0, color="gray", linestyle="--")
        ax6.set_ylabel(legend_dict['corr_eps_Hplus'], fontsize=10)
        ax6.tick_params(axis='both', which='major', labelsize=7)
        ax6.tick_params(axis='both', which='minor', labelsize=7)
        ax6.yaxis.set_major_formatter('{x:<5.3f}')
        ax6.text(0.02, 0.95, "F", transform=ax6.transAxes, fontsize=10, va='top')
        """
        # epsilon times Hplus plot
        epsHplus = -1 * epsHplus[i:]
        epsHplusstd = -1 * epsHplusstd[i:]
        ax6.plot(x_vals, epsHplus, color=colordict[system])
        ax6.fill_between(x_vals, (epsHplus - epsHplusstd), (epsHplus + epsHplusstd), alpha=a, color=colordict[system])
        # ax6.set_ylim(-.25, 0)
        ax6.axhline(0, color="gray", linestyle="--")
        ax6.set_ylabel(legend_dict['neg_avg_epsilon_H'], fontsize=10)
        ax6.tick_params(axis='both', which='major', labelsize=7)
        ax6.tick_params(axis='both', which='minor', labelsize=7)
        ax6.yaxis.set_major_formatter('{x:<5.3f}')
        ax6.text(0.02, 0.95, "F", transform=ax6.transAxes, fontsize=10, va='top')

    ax1.set_xlim(0, xmax / 10 + 1)
    ax5.set_xlabel(r'$r\;(\mathrm{nm})$', fontsize=10)
    ax6.set_xlabel(r'$r\;(\mathrm{nm})$', fontsize=10)

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    # lgd = fig.legend(lines, labels, loc="upper center", bbox_to_anchor=(0.5, 1.075), ncol=len(labels), title=r'$t_R/t_0$')
    lgd = fig.legend(lines, labels, loc="upper center", bbox_to_anchor=(0.5, 1.075), ncol=len(labels), title='Percent Mismatch')

    fig.tight_layout()
    plt.savefig("/home/js2746/Desktop/comparisonfig_" + comparison + ".pdf", bbox_inches='tight')
    # plt.show()
    plt.clf()
    plt.close()


def plot_avg_H2_over_time(system, path, coordsys):
    H1 = np.load(path + system + ".zone.C1A.C1B." + coordsys + ".meancurvature.npy")
    H2 = np.load(path + system + ".ztwo.C1A.C1B." + coordsys + ".meancurvature.npy")
    Hplus = (H1 + H2) / 2
    Hplus2 = Hplus * Hplus
    nframes = np.shape(Hplus)[2]
    fig, axs = plt.subplots(2, sharex=True, sharey=False, layout="constrained")
    for item in [Hplus, Hplus2]:
        x = np.linspace(0, nframes - 1, nframes)
        x = x / 10.
        y = np.nanmean(item, axis=0)
        y = np.nanmean(y, axis=0)
        # remove zero point
        x = x[1:]
        y = y[1:]
        rollavg = rollingavg(y, 10)
        if item is Hplus:
            figax = 0
        else:
            figax = 1
        axs[figax].plot(x, y)
        axs[figax].plot(x, rollavg, color="red")
    axs[0].set_ylabel(r'$\langle H^+ \rangle \;(\mathrm{\dot A^{-1}})$')
    axs[1].set_ylabel(r'$\langle ( H^+ )^2 \rangle\; (\mathrm{\dot A^{-2}})$')
    fig.supxlabel(r'$t \;(\mu s)$')
    plt.savefig("H2_over_time.pdf", dpi=700)
    plt.clf()
    plt.close()


def normalize_by_bulk(mols, paths, nougvals, bulkpath, bulknougvals, xlim):
    mainpath = "/home/js2746/Bending/PC/whole_mols"
    dirname = os.path.join(mainpath, "avg_over_theta")
    try:
        os.mkdir(dirname)
    except:
        pass

    for quantity in ["avg_epsilon2", "avg_H_plus2"]:
        for series, seriesname in zip([satsys, monounsatsys], ["satsys", "monounsatsys"]):
            fig, axs = plt.subplots(2, sharex=True, sharey=True, layout="constrained")
            axs[0].set_xlim(0, xlim)
            limname = "zoomed_" + str(xlim)
            if xlim == 6:
                axs[0].set_ylim(min_scale_dict[quantity + "_normed6"], max_scale_dict[quantity + "_normed6"])
            else:
                axs[0].set_ylim(min_scale_dict[quantity + "_normed20"], max_scale_dict[quantity + "_normed20"])
            for name in series:
                for mol, path, nougval in zip(mols, paths, nougvals):
                    if mol == "5x29":
                        sysname = "lg" + name
                        figax = 1
                    elif mol == "7k3g":
                        sysname = name + "PC"
                        figax = 0
                    exp = np.load(path + "/" + sysname + "/" + sysname + nougval + "/npy/" + sysname + "." + quantity + ".npy")
                    bulk = np.load(bulkpath + "/lg" + name + "/lg" + name + bulknougvals + "/npy/lg" + name + "." + quantity + ".npy")
                    bulkval = np.mean(bulk[20:-2, :])
                    # bulkval = np.mean(bulk)
                    print(bulkval, name, quantity, mol)
                    normed = exp / bulkval
                    '''
					except ValueError as error:
						explen = np.shape(exp)[0]
						bulklen = np.shape(bulk)[0]
						print(name+": experimental is "+str(explen)+" long but bulk is "+str(bulklen)+" long.")
						minval = min(explen, bulklen)
						exp = exp[0:minval,:]
						bulk = bulk[0:minval,:]
						normed = exp/bulk
					'''
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        avg_vals = np.nanmean(normed, axis=1)
                    maxval = len(avg_vals)
                    x = np.arange(2.5, (maxval * 5 + 2.5), 5) / 10
                    if mol == "5x29":
                        start = np.where(x == 2.75)[0][0]
                    elif mol == "7k3g":
                        start = np.where(x == 1.25)[0][0]
                    axs[figax].plot(x[start:], avg_vals[start:], color=colordict[name])
            axs[0].axvline(x=2.75, color='black', linestyle='--')
            # axs[1].axvline(x=1.25, color = 'black', linestyle='--')
            # axs[1].axvline(x=2.75, color = 'black', linestyle=':')
            axs[0].set_title("7k3g")
            axs[1].set_title("5x29")
            fig.supxlabel(r'$r \;(\mathrm{nm})$')
            if quantity == "avg_epsilon2":
                fig.supylabel(r'$\langle \tilde \epsilon^2 \rangle$')
            elif quantity == "avg_H_plus2":
                fig.supylabel(r'$\langle (\tilde H^+)^2 \rangle$')
            plt.savefig(mainpath + "/avg_over_theta/" + seriesname + "_" + quantity + "_normed_" + limname + ".pdf", dpi=700)
            plt.clf()
            plt.close()


def plot_APL(path, sysname):
    """
    Plot Area Per Lipid from a nougat output

    Parameters
    ----------
    path : string
        Path to directory containing ___.APL.traj file
    sysname : string
        Name of system that you gave nougat

    Returns
    -------
    None.

    """
    APL_traj = np.loadtxt(path + sysname + ".APL.traj")
    fig, ax = plt.subplots(layout='constrained')
    X = APL_traj[:, 0]
    area = APL_traj[:, 1]
    ax.plot(X, area, color="blue")
    ax.plot(X, rollingavg(area, 20), color="red")
    ax.legend(["Raw", "Rolling Average (20)"])
    fig.supxlabel(r'$t \;(\mathrm{frames})$')
    fig.supylabel(r'$\mathrm{Area / lipid} \;(\mathrm{\dot A^2})$')
    plt.savefig(path + sysname + ".APL_traj.pdf", dpi=700)
    plt.clf()
    plt.close()


def plot_asymm_over_traj(path, sysname):
    """
    Plot the number of lipids per leaflet, plus lipids that have been discounted \
        or discarded from the leaflet sorter.

    Parameters
    ----------
    path : string
        Path to directory containing ___.asymm.traj file
    sysname : string
        Name of system that you gave nougat

    Returns
    -------
    None.

    """
    asymm_traj = np.loadtxt(path + sysname + ".asymm.traj")
    fig, axs = plt.subplots(4, sharex=True, layout='constrained')
    X = asymm_traj[:, 0] / 10  # magic number to convert from frames to us for my systems
    outer = asymm_traj[:, 1]
    inner = asymm_traj[:, 2]
    sideways = asymm_traj[:, 3]
    inside_inclusion = asymm_traj[:, 4]
    total_asymm = outer - inner
    total = outer + inner
    f_outer = outer / total
    f_inner = inner / total
    axs[0].plot(X, outer, color="blue", label='Outer Leaflet')
    axs[0].plot(X, inner, color="red", label="Inner leaflet")
    axs[1].plot(X, sideways, color="green", label="Unclassifiable")
    axs[1].plot(X, inside_inclusion, color="orange", label="Inside inclusion")
    axs[2].plot(X, total_asymm, color="black", label="Lipid number asymmetry")
    axs[3].plot(X, f_outer, color="blue", label="Outer Leaflet Fraction")
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[3].legend()
    fig.supxlabel(r'$t \;(\mathrm{\mu s})$')
    fig.supylabel(r'$nL$')
    plt.savefig(path + sysname + ".asymm_traj.pdf", dpi=700)


def compare_APLs(names, paths):
    """
    Create plot of average APL for several different systems

    Parameters
    ----------
    names : list
        list of system names you want to compare.
    path : string
        path to folder containing subdirectories of $names with APL.traj inside.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(layout="constrained",figsize=(4,2))
    for name, path in zip(names,paths):
        APL_traj = np.loadtxt(path + "/" + name + ".area.traj")
        X = APL_traj[:, 0] / 100
        Y = APL_traj[:, 2] / 100
        ax.plot(X, Y, color=colordict[name])
    # ax.legend(["position restraints", "position restraints + higher APL", "elastic network"])
    fig.supxlabel(r'$t \;(\mathrm{\mu s})$')
    fig.supylabel(r'$\Sigma \;(\mathrm{nm}^2)$')
    plt.savefig("/home/js2746/Desktop/comparison.APL_traj.pdf")
    plt.clf()
    plt.close()


def compare_P(names, path):
    """
    Create plot of pressure tensor components for several different systems

    Parameters
    ----------
    names : list
        list of system names you want to compare.
    path : string
        path to folder containing subdirectories of $names with $name.xvg inside.

    Returns
    -------
    None.

    """

    for name in names:
        fig, ax = plt.subplots()
        P_traj = np.loadtxt(path + "/" + name + ".xvg")
        X = P_traj[:, 0] / 1000000
        X = X[1:]
        PXX = P_traj[:, 1]
        PXX = PXX[1:]
        PYY = P_traj[:, 2]
        PYY = PYY[1:]
        PZZ = P_traj[:, 3]
        PZZ = PZZ[1:]
        if name == "elastic_cp":
            ax.set_title("Elastic Network")
        elif name == "APL0.595_cp":
            ax.set_title("Position Restraints - Original")
        else:
            ax.set_title("Position Restraints - Higher APL")
        ax.plot(X, rollingavg(PXX, 2000), color="red")
        ax.plot(X, rollingavg(PYY, 2000), color="blue")
        ax.plot(X, rollingavg(PZZ, 2000), color="green")
        ax.plot(X, rollingavg((PZZ - 0.5 * (PXX + PYY)), 2000), color="purple")
        ax.legend(["PXX", "PYY", "PZZ", "Gamma"])
        fig.supxlabel(r'$t \;(\mathrm{\mu s})$')
        fig.supylabel(r'$\mathrm{Pressure} \;(\mathrm{bar})$')
        plt.savefig(path + "/" + name + "_comparison.P_traj.pdf", dpi=700)
        plt.clf()
        plt.close()


def plot_P_v_A(names, path):
    """
    Create scatter plot of pressure (Y) versus box area (X) for several different systems

    Parameters
    ----------
    names : list
        list of system names you want to compare.
    path : string
        path to folder containing subdirectories of $names with $name.xvg inside.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots()
    for name in names:
        P_traj = np.loadtxt(path + "/" + name + ".xvg")
        A_traj = np.loadtxt(path + "/" + name[:-3] + ".area.traj")
        P_reduced = np.zeros((np.shape(A_traj)[0], 4))
        for row in range(np.shape(P_traj)[0]):
            if P_traj[row, 0] % 10000 == 0:
                P_reduced[int(P_traj[row, 0] / 10000), :] = P_traj[row, :]

        A_traj[:, 1] = A_traj[:, 1] / 100
        if name == "APL0.67_cp":
            A_traj[:, 1] = A_traj[:, 1] / 2397
        else:
            A_traj[:, 1] = A_traj[:, 1] / 2690
        PXX = P_reduced[:, 1]
        PYY = P_reduced[:, 2]
        PZZ = P_reduced[:, 3]
        GAMMA = PZZ - 0.5 * (PXX + PYY)
        GAMMA = GAMMA[20:]
        A_traj = A_traj[20:, 1]
        GAMMA = np.mean(GAMMA)
        A_traj = np.mean(A_traj)
        ax.scatter(A_traj, GAMMA, label=name[:-3])
    fig.legend()
    fig.supxlabel(r'$\mathrm{Area~per~lipid} \;(\mathrm{nm}^2)$')
    fig.supylabel(r'$\mathrm{Surface~Tension} \;(\mathrm{?})$')
    plt.savefig(path + "/scatter.A_P_traj.pdf")
    plt.clf()
    plt.close()


def plot_APL_v_nL(names, path):
    """
    Create plot of average APL versus number of lipids in system.

    Parameters
    ----------
    names : list
        list of system names you want to compare.
    path : string
        path to folder containing subdirectories of $names with APL.traj inside.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(layout="constrained")
    Y = []
    errors = []
    for name in names:
        APL_traj = np.loadtxt(path + name + "/" + name + ".APL.traj")
        Y.append(np.nanmean(APL_traj[:, 1]))
        errors.append(np.nanstd(APL_traj[:, 1]))
    X = np.int64(names)
    ax.errorbar(X, Y, errors)
    fig.supxlabel(r'$\mathrm{nL}$')
    fig.supylabel(r'$\mathrm{Avg. Area / lipid} \;(\mathrm{\dot A^2})$')
    plt.savefig(path + "comparison2.APL_traj.pdf", dpi=700)
    plt.clf()
    plt.close()


APL_color_dict = {
    "512": "red",
    "1024": "orange",
    "2048": "yellow",
    "4096": "green",
    "8192": "blue",
    "32768": "purple",
    "COMtiltspin": "black",
    "new_gmx": "grey",
    "100": "red",
    "250": "orange",
    "500": "green",
    "1000": "blue",
    "5000": "purple"
}

colordict = {
    "DT": "red",
    "DL": "orange",
    "DX": "purple",
    "DB": "blue",
    "DY": "orange",
    "DO": "green",
    "PO": "green",
    "DP": "green",
    "DG": "blue",
    "DTPC": "red",
    "DLPC": "orange",
    "DXPC": "purple",
    "DBPC": "blue",
    "DYPC": "orange",
    "DOPC": "green",
    "POPC": "green",
    "DPPC": "green",
    "DGPC": "blue",
    "lgDT": "red",
    "lgDL": "orange",
    "lgDX": "purple",
    "lgDB": "blue",
    "lgDY": "orange",
    "lgDO": "green",
    "lgPO": "green",
    "lgDP": "green",
    "lgDG": "blue",
    "100kjmol": "red",
    "1000kjmol": "blue",
    "5000": "purple",
    "lgPO_15": "green",
    "COMtiltspin": "red",
    "new_gmx_pos": "blue",
    "posres": "tab:orange",
    "APL_0.6": "tab:blue",
    "elasticNVT": "tab:green",
    "refcoord_scaling_no": "tab:red"
}
max_scale_dict = {
    "avg_epsilon_over_t0": .1,
    "avg_abs_epsilon": 4,
    "avg_abs_epsilon_over_t0": 0.4,
    "avg_epsilon2_over_t02": .15,
    "avg_epsilon_H_over_t0": 0,
    "avg_epsilon2": 5,
    "avg_epsilon2_normed6": 35,
    "avg_H_plus2_normed6": 15,
    "avg_epsilon2_normed20": 35,
    "avg_H_plus2_normed20": 30,
    "avg_H_plus2": .006,
    "avg_tilde_t": 1.05,
    "avg_epsilon": 2,
    "avg_total_t": 40,
    "corr_mag_epst0_Hplus": 0,
    "corr_epst0_Hplus": 0,
    "avg_rms_epsilon_over_t0": 6
}
min_scale_dict = {
    "avg_epsilon_over_t0": -.15,
    "avg_abs_epsilon": 0,
    "avg_abs_epsilon_over_t0": 0,
    "avg_epsilon2_over_t02": 0,
    "avg_epsilon_H_over_t0": -.015,
    "avg_epsilon2": 0,
    "avg_epsilon2_normed6": 0,
    "avg_H_plus2_normed6": 0,
    "avg_epsilon2_normed20": 0,
    "avg_H_plus2_normed20": 0,
    "avg_H_plus2": 0,
    "avg_tilde_t": 0,
    "avg_epsilon": -3,
    "avg_total_t": 0,
    "corr_mag_epst0_Hplus": -0.02,
    "corr_epst0_Hplus": -0.02,
    "avg_rms_epsilon_over_t0": .95
}

mismatch_dict = {
    "lgDT": "-6%",
    "lgDL": "-44%",
    "lgDP": "-59%",
    "lgDB": "-68%",
    "lgDX": "-74%",
    "lgDY": "-38%",
    "lgDO": "-54%",
    "lgDG": "-64%",
    "DTPC": "-6%",
    "DLPC": "-44%",
    "DPPC": "-59%",
    "DBPC": "-68%",
    "DXPC": "-74%",
    "DYPC": "-38%",
    "DOPC": "-54%",
    "DGPC": "-64%"
}

stiffness_dict = {
    "100kjmol": "100 kJ/mol",
    "1000kjmol": "1000 kJ/mol",
    "5000": "5000 kJ/mol",
    "lgPO_15": "position restraints",
    "COMtiltspin": "elastic network",
    "new_gmx_pos": "position restraints"
}

legend_dict = {
    "avg_epsilon_over_t0": r'$\langle \epsilon / t_0 \rangle$',
    "avg_abs_epsilon": r'$\langle | \epsilon | \rangle\; (\mathrm{\dot A})$',
    "avg_abs_epsilon_over_t0": r'$\langle | \epsilon / t_0 | \rangle$',
    "avg_epsilon2_over_t02": r'$\langle ( \epsilon / t_0 )^2 \rangle$',
    "avg_epsilon_H_over_t0": r'$\langle \epsilon H^+ / t_0 \rangle\; (\mathrm{\dot A^{-1}})$',
    "avg_epsilon2": r'$\langle \epsilon^2 \rangle\; (\mathrm{\dot A^2})$',
    "avg_tilde_epsilon2": r'$\langle \tilde \epsilon ^ 2 \rangle$',
    "avg_rms_tilde_epsilon2": r'$ \sqrt{\langle \tilde \epsilon^2 \rangle}$',
    "avg_H_plus2": r'$\langle ( H^+ )^2 \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_tilde_H_plus2": r'$\langle (\tilde H^+)^ 2 \rangle$',
    "avg_rms_tilde_H_plus2": r'$ \sqrt{\langle(\tilde H^+)^2\rangle}$',
    "avg_tilde_total_t": r'$\langle \tilde t \rangle$',
    "avg_epsilon": r'$\langle \epsilon \rangle\; (\mathrm{\dot A})$',
    "avg_total_t": r'$\langle t \rangle\; (\mathrm{\dot A})$',
    "corr_mag_epst0_Hplus": r'$\langle \delta_{| \epsilon| | H^+ | / t_0} \rangle; (\mathrm{\dot A^{-1}})$',
    "corr_epst0_Hplus": r'$ \langle \delta_{\epsilon H^+ / t_0} \rangle; ( \mathrm{\dot A^{-1}} )$',
    "avg_rms_epsilon_over_t0": r'$\langle \mathrm{rms}\;\epsilon / t_0 \rangle$',
    "avg_K_plus": r'$\langle K^+ \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_K_minus": r'$\langle K^- \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_H_plus": r'$\langle H^+ \rangle\; (\mathrm{\dot A^{-1}})$',
    "avg_H_minus": r'$\langle H^- \rangle\; (\mathrm{\dot A^{-1}})$',
    "avg_H_minus2": r'$\langle \left ( H^- \right )^2 \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_epsilon_H": r'$\langle  \epsilon H^+  \rangle$',
    "neg_avg_epsilon_H": r'$-\langle  \epsilon H^+  \rangle$',
    "avg_z_minus": r'$\langle z^- \rangle\; (\mathrm{\dot A})$',
    "avg_z_minus2": r'$\langle \left ( z^- \right )^2 \rangle\; (\mathrm{\dot A^2})$',
    "avg_z_minus_H_minus": r'$\langle z^- H^- \rangle$',
    "avg_z_minus2_over_t02": r'$\langle \left ( z^- / t_0 \right )^2 \rangle$',
    "avg_z_minus_H_minus_over_t0": r'$\langle z^- H^- / t_0 \rangle\; (\mathrm{\dot A^{-1}})$',
    "corr_epst0_Kplus": r'$\langle \delta_{\epsilon / t_0, K^+} \rangle\; (\mathrm{\dot A^{-2}})$',
    "corr_eps_Kplus": r'$\langle \delta_{\epsilon, K^+} \rangle; (\mathrm{\dot A^{-1}})$',
    "corr_mag_eps_Hplus": r'$\langle \delta_{| \epsilon |, | H^+ |} \rangle$',
    "corr_eps_Hplus": r'$\langle  \delta_{\epsilon, H^+} \rangle$'
}

allsys = ["DT", "DL", "DP", "DB", "DX", "DY", "DO", "DG"]
satsys = ["DT", "DL", "DP", "DB", "DX"]
monounsatsys = ["DY", "DO", "DG"]
coordsys = "polar"
# coordsys = "cart"
mols = ["5x29", "7k3g"]
paths = ["/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1", "/home/js2746/Bending/PC/whole_mols/7k3g/lgSims"]
nougvals = ["_polar_5_10_100_-1_1", "_polar_5_10_0_-1_1"]
bulkpath = "/home/js2746/Bending/PC/whole_mols/empty"
bulknougvals = "_polar_5_10_0_-1_1"
area_paths = ["/home/jesse/research/COMtiltspin/", "/home/jesse/research/new_gmx_pos/", "/home/jesse/research/5x29_stiffness/100/", "/home/jesse/research/5x29_stiffness/250/", "/home/jesse/research/5x29_stiffness/500/", "/home/jesse/research/5x29_stiffness/1000/", "/home/jesse/research/5x29_stiffness/5000/"]
# area_paths = ["/home/jesse/research/COMtiltspin/", "/home/jesse/research/new_gmx_pos/"]
# area_paths = ["/home/jesse/research/COMtiltspin/", "/home/jesse/research/new_gmx_pos/",]

if __name__ == "__main__":
    # for mol, path in zip(mols, paths):
    # run_avg_over_theta(mol,path)
    # for xlim in [6,20]:
    # plot_together(mols, paths, nougvals, xlim)
    # normalize_by_bulk(mols, paths, nougvals, bulkpath, bulknougvals, xlim)
    # run_sum_over_H2(mol)
    # run_eps_corr_scatter(mol)
    # plot_avg_H2_over_time("lgPO", "/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_42us/lgPO_cart_10_10_0_-1_1/npy/", "cart")
    # make_2d_series_over_time("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO/lgPO_polar_5_10_0_-1_1", "zone.C1A.C1B.polar.thickness", "polar", "lgPO")
    # plot_APL("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_42us/", 'lgPO')
    # compare_APLs(["512", "1024", "2048", "4096", "8192", "32768"], "/home/js2746/KC_project/")
    # compare_APLs(["APL0.595", "APL0.67", "elastic"], "/home/js2746/Bending/PC/whole_mols/5x29/APL")
    # compare_P(["APL0.595_cp", "APL0.67_cp", "elastic_cp"], "/home/js2746/Bending/PC/whole_mols/5x29/APL")
    # plot_P_v_A(["APL0.595_cp", "APL0.67_cp", "elastic_cp"], "/home/js2746/Bending/PC/whole_mols/5x29/APL")
    # plot_APL_v_nL(["512", "1024", "2048", "4096", "8192", "32768"], "/home/js2746/KC_project/")
    # plot_asymm_over_traj("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO_50us/", 'lgPO_50us')
    #make_paper_writing_group_plot("unsat_elas")
    #make_paper_writing_group_plot("sat_elas")
    # make_paper_writing_group_plot("elastic")
    compare_APLs(["posres", "APL_0.6", "refcoord_scaling_no", "elasticNVT"],["/home/js2746/gromacs_cell_artifact/posres/APL_0.6", "/home/js2746/gromacs_cell_artifact/elastic/NPT/APL_0.6", "/home/js2746/gromacs_cell_artifact/posres/refcoord_scaling/no", "/home/js2746/gromacs_cell_artifact/elastic/NVT/APL_0.6"])