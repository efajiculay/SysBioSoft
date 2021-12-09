"""

                 This module is the param_slider module

The  sole purpose of this  module is  to visually  modify parameters and
compare the result to a plotted data as the estimate changes in the plot

The functions in this modules are;

1. load_data
2. param_ode_model
3. update_range
4. submit
5. update
5. param_ode_int

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from tkinter import filedialog
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox  # , RadioButtons

from BioSANS2020.propagation.propensity import propensity_vec
from BioSANS2020.propagation.propensity import propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals
from BioSANS2020.myglobal import mglobals as globals2


def load_data():
    """This function reads the trajectory file containing the data which
    should follow the following example format;

    time	A	B
    0.0	100.0	0.0
    0.25	88.24969025197632	11.750309748023732
    0.5	77.88007831231087	22.119921687689185
    0.75	68.72892784164061	31.27107215835944
    1.0	60.65306592491437	39.346934075085684
    1.25	53.526142785532	46.473857214468055
    1.5	47.236655135816875	52.76334486418318
    1.75	41.68620193454698	58.31379806545308
    2.0	36.78794415253036	63.21205584746969
    2.25	32.46524678349081	67.53475321650924
    ...

    The format  above have  a header where the first  column is time and
    all other columns are species or components. The rows are the values
    of measurements  corresponding to  the header. Each row is delimited
    by  tab character "\t". If  the data file is in excel, just copy the
    data  from  excel to a text editor, save  it with a  filename and it
    will already  be tab delimited. If there  are several  replicates of
    the  data, just append them to  the end without  the header and this
    function can still handle all replicates.

    time	A	B
    0.0	100.0	0.0                                 # first replicate
    ...                                             # continuation
    0.0	100.0	0.0                                 # second replicate
    ...                                             # continuation

    Args:
        file (string, optional): trajectory file. Defaults to None.

    Returns:
        tuple: tuple  of  data and  labels (data, labels). The data is a
            list of all  the trajectories in the file and the labels are
            the corresponding name of the columns in the trajectory
    """
    try:
        end_time = 0
        file = filedialog.askopenfilename(title="Select file")
        with open(file, "r") as f_file:
            data = []
            ddvar = []
            row1 = str(f_file.readline()).strip()
            slabels = row1.split("\t")[1:]
            for row in f_file:
                cols = [float(xvar) for xvar in row.split("\t")]
                if end_time > cols[0] and ddvar:
                    data.append(np.array(ddvar))
                    ddvar = []
                ddvar.append(cols)
                end_time = cols[0]
            data.append(np.array(ddvar))
            return (data, slabels)
    except BaseException:
        return None


def param_ode_model(z_var, t_var, sp_comp, ks_dict, r_dict, p_dict,
                    v_stoich, molar=False):
    """This fuction returns the differential equation of components with
    respect to time.

    Args:
        z_var (list): list of initial concentration
        t_var (list): time stamp
        sp_comp (dict): dictionary of appearance or position of species
            or component in the reaction tag of BioSANS topology file.

            For example;

                #REACTIONS
                A => B, -kf1
                B => C, kf2

            The value of sp_comp is

                sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}

                A appears in first reaction with index 0
                B appears in first and second reaction with index 0, 1
                C appears in second reaction with index 1
        ks_dict (dict): dictionary of rate constant that appears in each
            reactions. For the above reaction, the value of ks_dict is

            For example;

                #REACTIONS
                A => B , 0.3        # first reaction
                B <=> C, 0.1, 0.2   # second reaction

            The value of ks_dict is

                ks_dict = {
                    0 : [0.3],      # first reaction
                    1 : [0.1, 0.2]  # second reaction
                }

        r_dict (dict): dictionary of reactant stoichiometry. For example

            r_dict = {
                0: {'A': 1},  # first reaction, coefficient of A is 1
                1: {'B': 1}   # second reaction, coefficient of B is 1
            }

        p_dict (dict): dictionary of product stoichiometry. For example

            p_dict = {
                0: {'B': 1},  # first reaction, coefficient of B is 1
                1: {'C': 1}   # second reaction, coefficient of C is 1
            }
        v_stoich (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.

    Returns:
        dxdt (numpy.ndarray): value of time derivatives of components
    """
    spc = [s for s in sp_comp]
    conc = {spc[xvar]: z_var[xvar] for xvar in range(len(spc))}
    if not molar:
        d_prop = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        d_prop = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(v_stoich, d_prop).reshape(len(z_var))
    for xvar in globals2.conBoundary:
        ind = spc.index(xvar)
        dxdt[ind] = 0
    return dxdt


def update_range(dk_var, valc):
    """[summary]

    Args:
        dk_var ([type]): [description]
        valc ([type]): [description]
    """
    # print(dk_var, type(dk_var), valc, type(valc), sep="\n\n")
    dk_var.valmin = valc[0]
    dk_var.valmax = valc[1]
    dk_var.ax.set_xlim(dk_var.valmin, dk_var.valmax)


def submit(text, dk_var):
    valc = [float(xvar.strip()) for xvar in text.split(",")]
    update_range(dk_var, valc)


def update(ks_dict, kk_list, fig, slabels, lvar, sp_comp,
           r_dict, p_dict, v_stoich, molar, t_var, z_var):
    cc_var = 0
    for ih_var, _ in enumerate(ks_dict):
        if len(ks_dict[ih_var]) == 1:
            ks_dict[ih_var][0] = kk_list[cc_var].val
            cc_var = cc_var + 1
        else:
            ks_dict[ih_var][0] = kk_list[cc_var].val
            cc_var = cc_var + 1
            ks_dict[ih_var][1] = kk_list[cc_var].val
            cc_var = cc_var + 1
    res = odeint(param_ode_model, z_var, t_var,
                 args=(sp_comp, ks_dict, r_dict, p_dict, v_stoich, molar))
    for ih_var in range(len(z_var)):
        lvar[ih_var].set_ydata(res[:, ih_var])
        lvar[ih_var].set_label(slabels[ih_var])
    fig.canvas.draw_idle()


def param_ode_int(conc, t_var, sp_comp, ks_dict, r_dict, p_dict,
                  v_stoich, molar=False, rfile="", set_p=None):

    if set_p is None:
        set_p = []
    spc = [s for s in sp_comp]
    z_var = [conc[a] for a in sp_comp]
    c_miss = []
    for i, _ in enumerate(z_var):
        if z_var[i] < 0:
            c_miss.append(spc[i])
    k_miss = []
    for i, _ in enumerate(ks_dict):
        if len(ks_dict[i]) == 1:
            if ks_dict[i][0] < 0:
                k_miss.append((i, 0))
                ks_dict[i][0] = 2
        elif len(ks_dict[i]) == 2:
            if ks_dict[i][0] < 0:
                k_miss.append((i, 0))
                ks_dict[i][0] = 2
            if ks_dict[i][1] < 0:
                k_miss.append((i, 1))
                ks_dict[i][0] = 2
    # npar = len(c_miss) + len(k_miss)

    data2 = load_data()
    if data2:
        si_2 = data2[1]
        t_var = data2[0][0][:, 0]
        data = data2[0][0][:, 1:]

    get_globals(rfile)
    z_var = [conc[a] for a in sp_comp]
    slabels = [slabels for slabels in sp_comp]

    nz_var = []
    reserve_events_words = {"t", "time", "status",
                            "status2", "timer", "finish", "delay", "dtime"}
    for row in range(v_stoich.shape[0]):
        key = slabels[row].strip().split("_")[0]
        if key not in reserve_events_words:
            nz_var.append(row)
    sp_comp = {slabels[z_var]: sp_comp[slabels[z_var]] for z_var in nz_var}
    slabels = [slabels for slabels in sp_comp]

    scount = 0
    for ih_var, _ in enumerate(ks_dict):
        if len(ks_dict[ih_var]) == 1:
            scount = scount + 1
        else:
            scount = scount + 2

    res = odeint(param_ode_model, z_var, t_var,
                 args=(sp_comp, ks_dict, r_dict, p_dict, v_stoich, molar))
    fig, _ = plt.subplots(figsize=(12, 4))
    plt.subplots_adjust(left=0.1, right=0.70)
    if data2:
        plt.ylim(0, np.max(data) * 1.05)
    # if set_p[0] == True:
    if set_p[0]:
        plt.xscale("log")
    # if set_p[1] == True:
    if set_p[1]:
        plt.yscale("log")

    if data2:
        for ih_var, _ in enumerate(si_2):
            plt.plot(t_var, data[:, ih_var], lw=1, ls='--',
                     label=si_2[ih_var] + "_True")

    lvar = []
    for ih_var in nz_var:
        lvar.append(plt.plot(t_var, res[:, ih_var], lw=2,
                             label=slabels[ih_var])[0])
    plt.legend()

    axcolor = 'lightgoldenrodyellow'
    kk_list = []
    ck_list = []
    cc_var = 0
    for ih_var, _ in enumerate(ks_dict):
        strt1 = round(abs(ks_dict[ih_var][0] * 0.10), 2)
        ends1 = round(abs(ks_dict[ih_var][0] * 1.90), 2)
        if len(ks_dict[ih_var]) == 1:
            sk_plot = plt.axes([0.73, 0.93 - 0.05 * cc_var, 0.1, 0.04],
                               facecolor=axcolor)
            s_lider = Slider(sk_plot, 'kf' + str(ih_var + 1), strt1, ends1,
                             valinit=ks_dict[ih_var][0],
                             valstep=ks_dict[ih_var][0] / 100)
            kk_list.append(s_lider)
            ck_list.append(
                TextBox(plt.axes([0.93, 0.93 - 0.05 * cc_var, 0.08, 0.04]),
                        'range', initial=str(strt1) + ', ' + str(ends1)))
            cc_var = cc_var + 1
        else:
            sk_plot = plt.axes([0.73, 0.93 - 0.05 * cc_var, 0.1, 0.04],
                               facecolor=axcolor)
            s_lider = Slider(sk_plot, 'kf' + str(ih_var + 1), strt1, ends1,
                             valinit=ks_dict[ih_var][0],
                             valstep=ks_dict[ih_var][0] / 100)
            kk_list.append(s_lider)
            ck_list.append(
                TextBox(plt.axes([0.93, 0.93 - 0.05 * cc_var, 0.08, 0.04]),
                        'range', initial=str(strt1) + ', ' + str(ends1)))
            cc_var = cc_var + 1

            strt2 = round(abs(ks_dict[ih_var][1] * 0.10), 2)
            ends2 = round(abs(ks_dict[ih_var][1] * 1.90), 2)

            sk_plot = plt.axes([0.73, 0.93 - 0.05 * cc_var, 0.1, 0.04],
                               facecolor=axcolor)
            s_lider = Slider(sk_plot, 'kb' + str(ih_var + 1), strt2, ends2,
                             valinit=ks_dict[ih_var][1],
                             valstep=ks_dict[ih_var][1] / 100)
            kk_list.append(s_lider)
            ck_list.append(
                TextBox(plt.axes([0.93, 0.93 - 0.05 * cc_var, 0.08, 0.04]),
                        'range', initial=str(strt2) + ', ' + str(ends2)))
            cc_var = cc_var + 1

    for ih_var in range(len(ck_list)):
        kk_list[ih_var].on_changed(
            lambda val: update(ks_dict, kk_list, fig,
                               slabels, lvar, sp_comp, r_dict,
                               p_dict, v_stoich, molar, t_var, z_var)
        )

    for ih_var in range(len(ck_list)):
        exec('ck_list[' + str(ih_var)
             + '].on_submit(lambda v: submit(v, kk_list['
             + str(ih_var) + ']))',
             {'kk_list': kk_list, 'submit': submit, 'ck_list': ck_list})

    plt.show()
    return [0]
