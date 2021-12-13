# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

"""
                       This is the plot_traj module

This module process numeric trajectories and plots them

The following is the list of function for this module:

1. plot_traj
2. plot_traj2

"""

# import multiprocessing
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits import mplot3d
from numpy import max as npmax

from BioSANS2020.gui_functs.draw_figure import draw_figure


def plot_traj(data, slabels, items, plotted, mix_plot=True, logx=False,
              logy=False, normalize=False, si_ticked=None):
    """This function plots time series trajectories.
    Args:
        data : three dimensional data. The outermost index represents
               the trajectories. Each row of trajectory contains a 2 x 2
               numpy array where the rows is the time and the columns
               are the labels. The first column is time and the next few
               columns are the trajectories.
        slabels : labels of the trajectories in data
        items : 3 item list of [canvas, scroll_x, scroll_y]
        plotted : an array where the figures plotted are pushed i.e.
                  plotted.append([plt.gca(), fig, lines]). This is used
                  for managing the plots in the GUI.
        mix_plot : if true, the all data will be plotted together in one
                   plot otherwise in seperate plots
        logx : if true, uses logscale in x-axis
        logy : if true, uses logscale in y-axis
        normalize : normalized the data based on max(data)
        si_ticked : list of indexes in slabels to be shown in plot
    """
    miter = len(data)
    col = ['C' + str(i) for i in range(100)]
    if si_ticked is None:
        si_ticked = range(len(slabels))
    if mix_plot:
        plt.figure(figsize=(9.68, 3.95))
        plt.xlabel("time (sec)")
        plt.ylabel("conc")
        lines = []
        if logx:
            plt.xscale('log')
        if logy:
            plt.yscale('log')
        for j in range(miter):

            for i in si_ticked:
                if normalize:
                    line = plt.plot(
                        data[j][:, 0], data[j][:, i + 1]
                        / (npmax(data[j][:, i + 1]) + 1.0e-30), col[i])
                else:
                    line = plt.plot(data[j][:, 0], data[j][:, i + 1], col[i])
        plt.legend([slabels[i] for i in si_ticked])
        plt.tight_layout()
        fig = plt.gcf()
        lines.append(line)
        plotted.append([plt.gca(), fig, lines])
        # fig_canvas_agg = draw_figure(items, fig)
        draw_figure(items, fig)
        # plt.close()
    else:
        lines = []
        for i in range(len(slabels)):
            plt.figure(figsize=(9.68, 3.95))
            plt.xlabel("time")
            plt.ylabel("conc")
            if logx:
                plt.xscale('log')
            if logy:
                plt.yscale('log')
            for j in range(miter):
                if normalize:
                    line = plt.plot([t / 3600 for t in data[j][0]],
                                    data[j][1][:, i]
                                    / (npmax(data[j][1][:, i])
                                       + 1.0e-30), col[i])
                else:
                    line = plt.plot([t / 3600 for t in data[j][0]],
                                    data[j][1][:, i], col[i])
                lines.append(line)

            plt.legend([slabels[i] for i in si_ticked])
            plt.tight_layout()
            fig = plt.gcf()
            plotted.append([plt.gca(), fig, lines])
            # fig_canvas_agg = draw_figure(items, fig)
            draw_figure(items, fig)
            # plt.close()


def plot_traj2(data, slabels, items, plotted, logx=False, logy=False,
               normalize=False, xlabel="conc", ylabel="conc",
               zlabel="conc", trange=None):
    """This function plots time series trajectories or phase portriat.
    Args:
        data : three dimensional data. The outermost index represents
               the trajectories. Each row of trajectory contains a 2 x 2
               numpy array where the rows is the time and the columns
               are the labels. The first column is time and the next few
               columns are the trajectories.
        slabels : labels of the trajectories in data
        items : 3 item list of [canvas, scroll_x, scroll_y]
        plotted : an array where the figures plotted are pushed i.e.
                  plotted.append([plt.gca(), fig, lines]). This is used
                  for managing the plots in the GUI.
        logx : if true, uses logscale in x-axis
        logy : if true, uses logscale in y-axis
        normalize : normalized the data based on max(data)
        xlabel : label of x-axis
        ylabel : label of y-axis
        zlabel : label of y-axis
        trange : slice indexes of the trajectory i.e. -1000:-1
    """
    miter = len(data)
    # col = ['C' + str(i) for i in range(100)]
    plt.figure(figsize=(9.68, 3.95))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if trange:
        tran = str(trange).split(":")
        tran = slice(int(tran[0]), int(tran[1]))
    else:
        tran = slice(0, -1)
    if zlabel == "None":
        lines = []
        if logx:
            plt.xscale('log')
        if logy:
            plt.yscale('log')
        for j in range(miter):
            ind1 = -1 if xlabel == "time" else slabels.index(xlabel)
            ind2 = -1 if ylabel == "time" else slabels.index(ylabel)
            if normalize:
                line = plt.plot(data[j][tran, ind1 + 1],
                                data[j][tran, ind2 + 1] /
                                (npmax(data[j][tran, ind2 + 1])
                                 + 1.0e-30))
            else:
                line = plt.plot(data[j][tran, ind1 + 1],
                                data[j][tran, ind2 + 1])
        plt.tight_layout()
        fig = plt.gcf()
        lines.append(line)
        plotted.append([plt.gca(), fig, lines])
        # fig_canvas_agg = draw_figure(items, fig)
        draw_figure(items, fig)
        # plt.close()
    else:
        axf = plt.axes(projection='3d')
        for j in range(miter):
            ind1 = -1 if xlabel == "time" else slabels.index(xlabel)
            ind2 = -1 if ylabel == "time" else slabels.index(ylabel)
            ind3 = -1 if zlabel == "time" else slabels.index(zlabel)
            axf.plot3D(data[j][tran, ind1 + 1], data[j]
                       [tran, ind2 + 1], data[j][tran, ind3 + 1])
        axf.set_xlabel(xlabel)
        axf.set_ylabel(ylabel)
        axf.set_zlabel(zlabel)
        plt.show()
