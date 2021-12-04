# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

"""
                       This is the transform_data module

This module process numeric trajectories and either prints the result or
display the output into a text area embedded in a canvas or plots it

The following is the list of function for this module:

1. normalize
2. calc_cross_corr
3. calc_covariance2
4. calc_covariance
5. fano_factor
6. prob_density_calc
7. prob_density_calc2
8. prob_density_calc3
9. ave_traj_calc

"""

import matplotlib.pyplot as plt
import numpy as np
from BioSANS2020.analysis.numeric.sample_points import sample_points
from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.gui_functs.scrollable_text import INSERT
from BioSANS2020.gui_functs.draw_figure import draw_figure
# from mpl_toolkits.mplot3d import Axes3D
# from scipy.stats import norm
from BioSANS2020.myglobal import mglobals as globals2


def normalize(vect):
    """ returns the normalized form of the input vector v """
    norm = np.linalg.norm(vect)
    if norm == 0:
        return vect
    return vect / norm


def calc_cross_corr(edata, items):
    """This function calculates the cross correlation of edata[0] and
    returns a plot of the correlation as a function of lags.
    Args:
        edata : two dimensional array of data & labels (data, label).
                data[0] is time and data[1:] are trajectories
        items : three dimensional array of [canvas, scroll_x, scroll_y]
    """
    if len(edata[0]) != 3:
        ndata, slabels = sample_points(edata)
        # time = ndata[0]
        data_length = len(ndata)
        dlen = len(slabels)

        ddm = [0] * len(slabels)
        for ind in range(1, data_length):
            ith_data = ndata[ind]
            for jnd in range(dlen):
                ddm[jnd] = ith_data[jnd] + ddm[jnd]

        ddm = [ddm[jnd] / (data_length - 1) for jnd in range(dlen)]
        vvv = [[[0] for xvar in range(dlen)] for xvar in range(dlen)]
        for ind in range(1, data_length):
            ith_data = ndata[ind]
            for jnd in range(dlen):
                for knd in range(dlen):
                    vvv[jnd][knd] = vvv[jnd][knd] + np.correlate(
                        normalize(ith_data[jnd][-1000:] - ddm[jnd][-1000:]),
                        normalize(ith_data[knd][-1000:] - ddm[knd][-1000:]),
                        "full",
                    )
        vvv = np.array(vvv) / (data_length - 1)
        for jnd in range(dlen):
            for knd in range(dlen):
                plt.figure(figsize=(9.68, 3.95))
                plt.xlabel("lag")
                plt.ylabel(slabels[jnd] + "-" + slabels[knd] + " correlation")
                xvar = int(len(vvv[jnd, knd]) / 2)
                xvar = [-y for y in list(range(xvar))][::-1] \
                    + [0] + list(range(xvar))
                line = plt.plot(xvar, vvv[jnd, knd])
                fig = plt.gcf()
                plt.tight_layout()
                globals2.plotted.append([plt.gca(), fig, [line]])
                draw_figure(items, fig)
                # plt.close()


def calc_covariance2(edata):
    """This function calculates the covariance of edata[0] and prints
    the result in a terminal window.
    Args:
        edata : two dimensional array of data & labels [data, label].
                data[0] is time and data[1:] are trajectories
    """
    ndata, slabels = sample_points(edata)
    # time = ndata[0]
    data_length = len(ndata)
    dlen = len(slabels)

    ddm = [0] * len(slabels)
    for ind in range(1, data_length):
        data_slice = np.array(ndata[ind])[:, -100:]
        for jnd in range(dlen):
            ddm[jnd] = data_slice[jnd] + ddm[jnd]

    ddm = [ddm[jnd] / (data_length - 1) for jnd in range(dlen)]
    ddm = np.mean(np.array(ddm), 1).reshape(dlen, 1)

    cumulative_cov = 0
    for ind in range(1, data_length):
        data_slice = np.array(ndata[ind])[:, -100:] - ddm
        ddd = np.zeros((dlen, dlen))
        for jnd in range(dlen):
            for knd in range(dlen):
                ddd[jnd, knd] = np.mean(data_slice[jnd] * data_slice[knd])
        cumulative_cov = cumulative_cov + ddd
    cumulative_cov = cumulative_cov / (data_length - 1)
    print("covariance")
    print(slabels)
    print(cumulative_cov)
    print()
    vvv = np.zeros((dlen, dlen))
    for ind in range(dlen):
        for jnd in range(dlen):
            vvv[ind, jnd] = cumulative_cov[ind, jnd] \
                / np.sqrt(cumulative_cov[ind, ind] * cumulative_cov[jnd, jnd])
    print("\ncorrelation")
    print(slabels)
    print(vvv)


def calc_covariance(edata, items, points=100):
    """This function calculates the covariance of edata[0] and prints
    the result in a text area embedded in a canvas.
    Args:
        edata : two dimensional array of data & labels [data, label].
                data[0] is time and data[1:] are trajectories
        items : three dimensional array of [canvas, scroll_x, scroll_y]
        points : last number of points considered in covariance
         calculation from -points to the end of array ([-points:] slice)
    """
    data, slabels = edata
    data_length = len(data)
    ddm = 0
    for ind in range(1, data_length):
        data_slice = data[ind][-points:, 1:]
        ddm = ddm + data_slice

    ddm = np.mean(ddm / (data_length - 1), 0)
    print(ddm)

    cumulative_cov = 0
    dlen = len(slabels)
    for knd in range(1, data_length):
        data_slice = data[knd][-points:, 1:] - ddm
        ddd = np.zeros((dlen, dlen))
        for ind in range(dlen):
            for jnd in range(dlen):
                ddd[ind, jnd] = np.mean(data_slice[:, ind]
                                        * data_slice[:, jnd])
        cumulative_cov = cumulative_cov + ddd
    cumulative_cov = cumulative_cov / (data_length - 1)

    text = prepare_scroll_text(items)
    text.insert(INSERT, "covariance\n\n")
    for ind in range(dlen):
        for jnd in range(ind, dlen):
            if str(cumulative_cov[ind, jnd]) not in {"None", "nan", "0.0"}:
                label = slabels[ind] + " and " + slabels[jnd]
                text.insert(INSERT, label.ljust(50) + " = "
                            + str(cumulative_cov[ind, jnd]) + "\n")

    vvv = np.zeros((dlen, dlen))
    text.insert(INSERT, "\ncorrelation\n\n")
    for ind in range(dlen):
        for jnd in range(ind, dlen):
            if str(cumulative_cov[ind, jnd]) not in {"None", "nan", "0.0"}:
                vvv[ind, jnd] = cumulative_cov[ind, jnd] \
                                / np.sqrt(cumulative_cov[ind, ind]
                                          * cumulative_cov[jnd, jnd])
                label = slabels[ind] + " and " + slabels[jnd]
                text.insert(INSERT, label.ljust(50) + " = "
                            + str(vvv[ind, jnd]) + "\n")

    data_slice = 0
    for ind in range(data_length):
        data_slice = data_slice + (data[ind][-points:, 1:] - ddm) ** 2
    fano_f = np.mean(data_slice / (data_length - 1), 0) / ddm
    text.insert(INSERT, "\nFano Factor\n\n")
    for ind in range(dlen):
        if str(fano_f[ind]) not in {"None", "nan"}:
            text.insert(INSERT, slabels[ind].ljust(50) + " = "
                        + str(fano_f[ind]) + "\n")

    globals2.plot_i = globals2.plot_i + 1


def fano_factor(edata, items, points=100):
    """This function calculates the fano-factor of edata[0] and prints
    the result in a text area embedded in a canvas.
    Args:
        edata : two dimensional array of data & labels [data, label].
                data[0] is time and data[1:] are trajectories
        items : three dimensional array of [canvas, scroll_x, scroll_y]
        points : last number of points considered in fano-factor
         calculation from -points to the end of array ([-points:] slice)
    """
    data, slabels = edata
    data_length = len(data)
    ddm = 0
    for ind in range(data_length):
        data_slice = data[ind][-points:, 1:]
        ddm = ddm + data_slice

    ddm = np.mean(ddm / (data_length - 1), 0)

    data_slice = 0
    for ind in range(data_length):
        data_slice = data_slice + (data[ind][-points:, 1:] - ddm) ** 2
    fano_f = np.mean(data_slice / (data_length - 1), 0) / ddm

    text = prepare_scroll_text(items)
    dlen = len(slabels)
    text.insert(INSERT, "Fano Factor\n\n")
    for ind in range(dlen):
        if str(fano_f[ind]) not in {"None", "nan"}:
            text.insert(INSERT, slabels[ind].ljust(50) + " = "
                        + str(fano_f[ind]) + "\n")

    globals2.plot_i = globals2.plot_i + 1


def prob_density_calc(edata, items):
    """This function calculates the probability density of edata[0] and
    returns a plot of the probability density.
    Args:
        edata : two dimensional array of data & labels (data, label).
                data[0] is time and data[1:] are trajectories
        items : three dimensional array of [canvas, scroll_x, scroll_y]
    """
    data, slabels = edata
    data_length = len(data)
    dlen = len(slabels)

    data_slice = data[0][:, 1:]
    for ind in range(1, data_length):
        data_slice = np.concatenate((data_slice, data[ind][:, 1:]), axis=0)

    for jnd in range(dlen):
        plt.figure(figsize=(9.68, 3.95))
        plt.xlabel("conc(" + slabels[jnd] + ")")
        plt.ylabel("Prob")
        line = plt.hist(data_slice[:, jnd], bins=100, density=True,
                        orientation="vertical")
        plt.legend([slabels[jnd]])
        plt.tight_layout()
        fig = plt.gcf()
        globals2.plotted.append([plt.gca(), fig, [line]])
        draw_figure(items, fig)
        # plt.close()


def prob_density_calc2(edata, items):
    """This function calculates the probability density of edata[0] and
    returns a plot of the probability density.
    Args:
        edata : two dimensional array of data & labels (data, label).
                data[0] is time and data[1:] are trajectories
        items : three dimensional array of [canvas, scroll_x, scroll_y]
    """
    data, slabels = edata
    data_length = len(data)
    dlen = len(slabels)

    data_slice = data[0][:, 1:]
    data_list = data[0][:, 0].flatten()
    # tlen = len(data_list)
    for ind in range(1, data_length):
        data_slice = np.concatenate((data_slice, data[ind][:, 1:]), axis=0)
        data_list = np.concatenate((data_list, data[ind][:, 0]))

    for jnd in range(dlen):
        # r = data_slice
        zzs, xxs, yys = np.histogram2d(data_list, data_slice[:, jnd].flatten(),
                                       bins=(100, 100))
        plt.figure(figsize=(9.68, 3.95))
        plt.xlabel("time")
        plt.ylabel("conc(" + slabels[jnd] + ")")
        cntr = plt.contour(
            xxs[1:],
            yys[1:],
            zzs.T,
            list(range(int(np.min(zzs)), int(np.max(zzs)) + 1, 2)),
            cmap="Spectral_r",
        )
        fig = plt.gcf()
        fig.colorbar(cntr)
        plt.tight_layout()
        globals2.plotted.append([plt.gca(), fig, [cntr]])
        draw_figure(items, fig)
        # plt.close()


def prob_density_calc3(edata, items, bins=50):
    """This function calculates the probability density of edata[0] per
    bins and returns a plot of the probability density.
    Args:
        edata : two dimensional array of data & labels (data, label).
                data[0] is time and data[1:] are trajectories
        items : three dimensional array of [canvas, scroll_x, scroll_y]
        bins : number of bins an entire trajectory will be discretized
    """
    if len(edata[0]) != 3:
        ndata, slabels = sample_points(edata)
        time = ndata[0]
        data_length = len(ndata)
        dlen = len(slabels)

        ddd = [[] for xvar in range(dlen)]
        for knd in range(1, data_length):
            data_slice = ndata[knd]
            for jnd in range(dlen):
                ddd[jnd].append(data_slice[jnd])
        ddd = np.array(ddd)

        vvv = {}
        for ind in range(dlen):
            for jnd in range(len(time)):
                vvv[(ind, jnd)] = np.histogram(ddd[ind][:, jnd], bins=bins)

        tlen = len(time)
        strt = int(tlen / 10)
        for ind in range(dlen):
            lines = []
            fig = plt.figure(figsize=(9.68, 3.95))
            axf = fig.gca(projection="3d")
            for jnd in range(strt, tlen, int(tlen / 5)):
                xxs = vvv[(ind, jnd)][1][1:]
                yys = vvv[(ind, jnd)][0]
                width = xxs[1] - xxs[0]
                line = axf.bar(xxs, yys, zs=time[jnd], zdir="xvar",
                               width=max(width, 0.8), alpha=1)
                lines.append(line)
            axf.set_xlabel("time")
            axf.set_ylabel("conc(" + slabels[ind] + ")")
            axf.set_zlabel("freq")
            axf.view_init(elev=40, azim=-120)
            globals2.plotted.append([axf, fig, lines])
            draw_figure(items, fig)
            # plt.close()


def ave_traj_calc(edata, items):
    """This function calculates the average trajectory of edata[0] and
    returns a plot of the average trajectory as a function of time.
    Args:
        edata : two dimensional array of data & labels (data, label).
                data[0] is time and data[1:] are trajectories
        items : three dimensional array of [canvas, scroll_x, scroll_y]
    """
    if len(edata[0]) != 3:
        ndata, slabels = sample_points(edata)
        time = ndata[0]
        data_length = len(ndata)
        dlen = len(slabels)

        ddm = [0] * len(slabels)
        for ind in range(1, data_length):
            data_slice = ndata[ind]
            for jnd in range(dlen):
                ddm[jnd] = data_slice[jnd] + ddm[jnd]

        ddm = [ddm[jnd] / (data_length - 1) for jnd in range(dlen)]

        plt.figure(figsize=(9.68, 3.95))
        plt.xlabel("time (sec)")
        plt.ylabel("conc")
        lines = []
        for jnd in range(dlen):
            line = plt.plot(time, ddm[jnd], label=slabels[jnd])
            lines.append(line)
        plt.legend()
        fig = plt.gcf()
        plt.tight_layout()
        globals2.plotted.append([plt.gca(), fig, lines])
        draw_figure(items, fig)
        # plt.close()
