"""
Created on Mon Jul 17 10:54:23 2023.

@author: js2746
"""

import matplotlib.pyplot as plt
import numpy as np
import warnings
import os
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
                plot_maker(dim1vals, dim2vals, data[:15, :], system, field, max_scale_dict[value], min_scale_dict[value], False, value, False, "polar")
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


def make_2d_series_over_time(path, quantity, coordsys, sys_name):
    """
    Make movie file of heatmap over trajectory.

    Parameters
    ----------
    path : string
        the path the to nougat outputs directory you want.
    quantity : string
        the name of the measurement.
    coordsys : string
        "polar" or "cart".
    sys_name : string
        the name you gave nougat.py when it made your files.

    Returns
    -------
    None.

    """
    cwd = os.getcwd()
    if cwd != path:
        os.chdir(path)

    # load the correct 2d data
    traj_data = np.load(path + "/npy/" + sys_name + "." + quantity + ".npy")
    nframes = np.shape(traj_data)[2]
    dims = bin_prep(sys_name, "C1A.C1B", coordsys, 'OFF')
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
        plot_maker(dim1vals, dim2vals, frame_data, sys_name, str(frame), "auto", "auto", False, quantity, False, coordsys)

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
                    npzfile = np.load(path + "/" + sysname + "/" + sysname + nougval + "/npy/" + sysname + "." + quantity + "avg_over_theta.npz")
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


def plot_avg_H2_over_time(system, path):
    H1 = np.load(path + system + ".zone.C1A.C1B.polar.meancurvature.npy")
    H2 = np.load(path + system + ".ztwo.C1A.C1B.polar.meancurvature.npy")
    Hplus = (H1 + H2) / 2
    Hplus2 = Hplus * Hplus
    nframes = np.shape(Hplus)[2]
    fig, axs = plt.subplots(2, sharex=True, sharey=False, layout="constrained")
    for item in [Hplus, Hplus2]:
        x = np.linspace(0, nframes - 1, nframes)
        x = x / 10.
        y = np.nanmean(item, axis=0)
        y = np.nanmean(y, axis=0)
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


colordict = {
    "DT": "red",
    "DL": "orange",
    "DX": "purple",
    "DB": "blue",
    "DY": "orange",
    "DO": "green",
    "PO": "green",
    "DP": "green",
    "DG": "blue"
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

legend_dict = {
    "avg_epsilon_over_t0": r'$\langle \epsilon / t_0 \rangle$',
    "avg_abs_epsilon": r'$\langle | \epsilon | \rangle\; (\mathrm{\dot A})$',
    "avg_abs_epsilon_over_t0": r'$\langle | \epsilon / t_0 | \rangle$',
    "avg_epsilon2_over_t02": r'$\langle ( \epsilon / t_0 )^2 \rangle$',
    "avg_epsilon_H_over_t0": r'$\langle \epsilon H^+ / t_0 \rangle\; (\mathrm{\dot A^{-1}})$',
    "avg_epsilon2": r'$\langle \epsilon^2 \rangle\; (\mathrm{\dot A^2})$',
    "avg_H_plus2": r'$\langle ( H^+ )^2 \rangle\; (\mathrm{\dot A^{-2}})$',
    "avg_tilde_t": r'$\langle \tilde t \rangle$',
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
    "corr_epst0_Kplus": r'$\langle \epsilon K^+ / t_0 \rangle - \langle \epsilon / t_0 \rangle \langle K^+ \rangle\; (\mathrm{\dot A^{-2}})$'
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

if __name__ == "__main__":
    # for mol, path in zip(mols, paths):
    # run_avg_over_theta(mol,path)
    # for xlim in [6,20]:
    # plot_together(mols, paths, nougvals, xlim)
    # normalize_by_bulk(mols, paths, nougvals, bulkpath, bulknougvals, xlim)
    # run_sum_over_H2(mol)
    # run_eps_corr_scatter(mol)
    # plot_avg_H2_over_time("lgPO", "/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/mixed_dm_flags/lgPO/lgPO_polar_5_10_0_-1_1/npy/")
    make_2d_series_over_time("/home/js2746/Bending/PC/whole_mols/5x29/40nmSystems/dm1/lgPO/lgPO_polar_5_10_0_-1_1", "zone.C1A.C1B.polar.thickness", "polar", "lgPO")
