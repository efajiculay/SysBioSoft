"""

                        This is the process module

This module reads BioSANS topology file, grab the components or species,
rate  constants, stoichiometric  matrix,  propensity   vector, algebraic
rules, conditional statements, and othe types of definitions into a dic-
tionary.  This module  calls the process_hub  module that distribute the
tasks to other modules.

The functions in this module are as follows; #

1. eval_dict
2. tofloat
3. is_number
4. process


"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

# import warnings
from tkinter import messagebox as message_upon_error
# import re
# warnings.filterwarnings('ignore')

# import time
import numpy as np

from BioSANS2020.prepcodes.processes_hub import process_hub
# from BioSANS2020.myglobal import mglobals as globals2
# from BioSANS2020.myglobal import proc_global
from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT


def eval_dict(to_eval, loc_dict):
    """This function takes a string expression and return the evaluated
    expression using SBML_FUNCT_DICT and the locals() dictionary where
    eval_dict is called.

    Args:
        to_eval (str): the expression to evaluate
        loc_dict (dict): local dictionary from the calling function

    Returns:
        multitype: result of eval command
    """
    return eval(to_eval, loc_dict, SBML_FUNCT_DICT)


def tofloat(val, loc_dict):
    """This function attempts to convert the input val into float

    Args:
        val (str): the expression to evaluate
        loc_dict (dict): local dictionary from the calling function

    Returns:
        float: float equivalent of val
    """
    try:
        return float(val)
    except ValueError:
        return float(eval_dict(val, loc_dict))


def is_number(xvar):
    """This function checks if a string xvar is float

    Args:
        xvar (str): input string expression or number

    Returns:
        bool: True if xvar cal be converted to float otherwise False
    """
    try:
        float(xvar)
        return True
    except ValueError:
        return False


def process(
        rfile="Reactions",
        miter=1,
        conc_unit="molecules",
        v_volms=1.0e-20,
        tend=1,
        del_coef=10,
        normalize=False,
        logx=False,
        logy=False,
        method="CLE",
        tlen=1000,
        mix_plot=True,
        save=True,
        out_fname="",
        plot_show=True,
        time_unit="time (sec)",
        vary="",
        vary2="",
        mult_proc=False,
        implicit=False,
        items=None,
        exp_data_file=None,
        c_input={}
):
    """[summary]

    Args:
        rfile (str): file name of BioSANS topology file.
        miter (int, optional): Number of iteration or trajectory samples
            for stochastic integration
        conc_unit (bool, optional): "mole","molar", or "molecules" - the
            unit used in any amount in topology file. Defaults to
            "molecules".
        v_volms (float, optional): the volume of compartment used in the
            simulation. Defaults to 1.0e-20.
        tend (float): trajectory simulation end time. Defaults to 1.
        del_coef (float, optional): factor for modifying time steps used
            in the integration/propagation of ODE. Defaults to 10.
        normalize (bool, optional): True  will  be normalized the y axis
            based on max value . Defaults to False.
        logx (bool, optional): If True, the x-axis will be in log scale.
            Defaults to False.
        logy (bool, optional): if True, the y-axis will be in log scale.
            Defaults to False.
        method (str, optional): Defaults to "CLE". Any of the option in
            the list of available method keywords is listed below;

            Stochastic (refer to section 10.2.4)

            1.	"CLE"            - Molecules(micro), tau-adaptive
            2.	"CLE2"           - Molecules(micro), cle-fixIntvl
            3.	"Gillespie_"     - Molecules(micro), Direct method
            4.	"Tau-leaping"    - Molecules(micro),
                                   Not swapping with Gillespie
            5.	"Tau-leaping2"   - Molecules(micro),
                                   Swapping with Gillespie
            6.	"Sim-TauLeap"    - Molecules(micro), Simplified,
                                   Swapping with Gillespie

            Deterministic (refer to section 10.2.1)

            7.	"Euler-1"        - Molecules(micro), tau-adaptive-1
            8.	"Euler-2"        - Molar (macro), tau-adaptive-1
            9.	"Euler-3"        - Mole (macro), tau-adaptive-1
            10.	"Euler2-1"	     - Molecules(micro), tau-adaptive-2
            11.	"Euler2-2"       - Molar (macro), tau-adaptive-2
            12.	"Euler2-3"       - Mole (macro), tau-adaptive-2
            13.	"ODE-1"          - Molecules(micro),
                                   using ode_int from scipy
            14.	"ODE-2"          - Molar(macro),
                                   using ode_int from scipy
            15.	"ODE-3"          - Mole(macro), using ode_int from scipy
            16.	"rk4-1"          - Molecules(micro), fix-interval
            17.	"rk4-2"          - Molar(macro), fix-interval
            18.	"rk4-3"          - Mole(macro), fix-interval
            19.	"rk4-1a"         - Molecules(micro), tau-adaptive
            20.	"rk4-2a"         - Molar(macro), tau-adaptive
            21.	"rk4-3a"         - Mole(macro), tau-adaptive

            Linear Noise Approximation (refer to 10.1.2 & 10.2.2)

            22.	"LNA"             - Numeric, values
            23.	"LNA-vs"          - Symbolic, values, Macroscopic
            24.	"LNA-ks"          - Symbolic, f(ks), Macroscopic
            25.	"LNA-xo"          - Symbolic, f(xo), Macroscopic
            26.	"LNA2"            - Symbolic, f(xo,ks), Microscopic
            27.	"LNA3"            - Symbolic, f(xo,ks), Macroscopic
            28.	"LNA(t)"          - COV-time-dependent, Macroscopic
            29.	"LNA2(t)"         - FF-time-dependent, Macroscopic

            Network Localization (refer to 10.1.3)

            30.	"NetLoc1"         - Symbolic, Macroscopic
            31.	"NetLoc2"         - Numeric, Macroscopic

            Parameter estimation (refer to 10.2.3)

            32.	"k_est1"          - MCEM, Macroscopic
            33.	"k_est2"          - MCEM, Microscopic
            34.	"k_est3"          - NM-Diff. Evol., Macroscopic
            35.	"k_est4"          - NM-Diff. Evol., Microscopic
            36.	"k_est5"          - Parameter slider/scanner
            37.	"k_est6"          - Nelder-Mead (NM), Macroscopic
            38.	"k_est7"          - Nelder-Mead (NM), Microscopic
            39.	"k_est8"          - Powell, Macroscopic
            40.	"k_est9"          - Powell, Microscopic
            41.	"k_est10"         - L-BFGS-B, Macroscopic
            42.	"k_est11"         - L-BFGS-B, Microscopic

            Symbolic/Analytical expression of species (refer to 10.1.1)

            43.	"Analyt"          - Pure Symbolic :f(t,xo,k)
            44.	"Analyt-ftx"      - Semi-Symbolic :f(t,xo)
            45.	"SAnalyt"         - Semi-Symbolic :f(t)
            46.	"SAnalyt-ftk"     - Semi-Symbolic :f(t,k)
            47.	"Analyt2"         - Creates commands for wxmaxima

        tlen (int, optional): number of integration steps reported in
            the final result. Defaults to 1000.
        mix_plot (bool, optional): If True, all species are plotted in
            one plot/figure. Defaults to True.
        save (bool, optional): If True, the resulting trajectory in the
            simulation will be saved as a file. Defaults to True.
        out_fname (str, optional): output filename. Defaults to "".
        plot_show (bool, optional): If True, an image of the plots will
            be created in the directory of the topology file.
        time_unit (str, optional): Defaults to "time (sec)".
        vary (str, optional): Varying initial concentration.
            Defaults to "".
        vary2 (str, optional): [description]. Defaults to "".
        mult_proc (bool, optional): If True, trajectories will be propa-
            gated on parallel. Defaults to False.
        implicit (bool, optional): True means report in time intervals
            similar to the input time intervals even if actual step is
            more or less. Defaults to False.
        items (tuple, optional): (canvas, scroll_x, scroll_y).
            Defaults to None.
        exp_data_file ([type], optional): Experimental data file contai-
            ning True or accepted trajectories. Defaults to None.
        c_input (dict, optional): [description]. Defaults to {}.

    Returns:
        list: list of simulated trajecotry.
            data[j][0] - time for trajectory j
            data[j][1][:, i] - trajectories of each component i
    """

    avogadros_num = 6.022e+23

    # globals2.MODIFIED = {}
    # globals2.PROP_MODIFIED = {}

    # for func in globals2.EXEC_FUNCTIONS:
    #     exec(func, globals(), SBML_FUNCT_DICT)

    if conc_unit == "molar":
        if method in ["ODE-1", "rk4-1", "rk4-1a", "Euler-1",
                      "Euler2-1", "CLE", "CLE2", "Tau-leaping",
                      "Tau-leaping2", "Sim-TauLeap", "Gillespie_"]:
            avgdrs_times_vol = avogadros_num * v_volms
            factor_stch = 2
        elif method in ["ODE-3", "rk4-3", "rk4-3a",
                        "Euler-3", "Euler2-3"]:
            avgdrs_times_vol = v_volms
            factor_stch = 1
        else:
            avgdrs_times_vol = 1
            factor_stch = 1
    elif conc_unit == "molecules":
        if method in ["ODE-2", "rk4-2", "rk4-2a",
                      "Euler-2", "Euler2-2"]:
            avgdrs_times_vol = 1 / avogadros_num * v_volms
            factor_stch = 1 / 2
        elif method in ["ODE-3", "rk4-3", "rk4-3a",
                        "Euler-3", "Euler2-3"]:
            avgdrs_times_vol = 1 / avogadros_num
            factor_stch = 1 / 2
        else:
            avgdrs_times_vol = 1
            factor_stch = 1
    elif conc_unit == "moles":
        if method in ["ODE-1", "rk4-1", "rk4-1a", "Euler-1", "Euler2-1",
                      "CLE", "CLE2", "Tau-leaping", "Tau-leaping2",
                      "Sim-TauLeap", "Gillespie_"]:
            avgdrs_times_vol = avogadros_num
            factor_stch = 2
        if method in ["ODE-2", "rk4-2", "rk4-2a",
                      "Euler-2", "Euler2-2"]:
            avgdrs_times_vol = 1 / v_volms
            factor_stch = 1
        else:
            avgdrs_times_vol = 1
            factor_stch = 1
    with open(rfile, "r") as file:
        rows = []
        conc = {}

        last = ""
        for row in file:
            if last == "Function_Definitions":
                if row.strip() != "" and row[0] != "#":
                    exec(row.strip(), locals(), SBML_FUNCT_DICT)
                elif row[0] == "#":
                    last = "#"
            elif last == "#":
                if row.strip() != "" and row[0] != "@":
                    rows.append(row)
                elif row.strip() != "" and row[0] != "#":
                    last = "@"
            elif last == "@":
                if row.strip() != "" and row[0] != "@":
                    cvar = row.split(",")
                    conc[cvar[0].strip()] = tofloat(cvar[1], locals())
                    # if len(cvar)>=3:
                    # cc = ",".join(cvar[2:])
                    # cc2 = cc.split(":")[0].replace("lambda","") \
                    #     .split(",")
                    # globals2.MODIFIED[cvar[0].strip()] = \
                    #     [cc2,eval_dict(cc)]
            elif row[0] == "#":
                last = "#"
                # gg = row.split(",")[1:]
                # try:
                # for xvar in gg:
                # xx = xvar.split("=")
                # globals2.SETTINGS[xx[0].strip()] = xx[1].strip()
                # except:
                # pass
            elif row[0] == "@":
                last = "@"
            elif row.strip() == "Function_Definitions:":
                last = "Function_Definitions"

        file.close()

    ks_dict = {}
    r_dict = {}
    p_dict = {}
    sp_comp = {}
    rxn_rows = len(rows)
    for ih_ind in range(rxn_rows):
        r_dict[ih_ind] = {}
        p_dict[ih_ind] = {}
        col_row = rows[ih_ind].split(":::::")
        row = col_row[0].strip().split(",")
        # if len(col_row)>1:
        # krow = col_row[1].strip().split(":::")
        # if len(krow)==2:
        # cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
        # cc3 = krow[1].split(":")[0].replace("lambda","").split(",")
        # globals2.PROP_MODIFIED["Prop_"+str(ih_ind)] = [
        #     (cc2,eval_dict(krow[0])),(cc3,eval_dict(krow[1]))]
        # else:
        # cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
        # globals2.PROP_MODIFIED["Prop_"+str(ih_ind)] = \
        #     [(cc2,eval_dict(krow[0]))]

        if len(row) == 3:
            ks_dict[ih_ind] = [
                tofloat(row[1], locals()),
                tofloat(row[2], locals())]
        else:
            ks_dict[ih_ind] = [tofloat(row[1], locals())]
        col_var = row[0].split("<=>")
        if len(col_var) == 1:
            col_var = row[0].split("=>")

        sp_c = col_var[0]
        svar = sp_c.strip().split()
        if len(svar) > 1:
            last = 1
            for xvar in svar:
                if not is_number(xvar) and xvar != "+":
                    r_dict[ih_ind][xvar] = last
                    last = 1
                    if xvar in sp_comp:
                        sp_comp[xvar].add(ih_ind)
                    else:
                        sp_comp[xvar] = {ih_ind}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = svar[0]
            r_dict[ih_ind][xvar] = 1
            if xvar in sp_comp:
                sp_comp[xvar].add(ih_ind)
            else:
                sp_comp[xvar] = {ih_ind}

        sp_c = col_var[1]
        svar = sp_c.strip().split()
        if len(svar) > 1:
            last = 1
            for xvar in svar:
                if (not is_number(xvar) or xvar.lower() == "e") \
                        and xvar != "+":
                    p_dict[ih_ind][xvar] = last
                    last = 1
                    if xvar in sp_comp:
                        sp_comp[xvar].add(ih_ind)
                    else:
                        sp_comp[xvar] = {ih_ind}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = svar[0]
            p_dict[ih_ind][xvar] = 1
            if xvar in sp_comp:
                sp_comp[xvar].add(ih_ind)
            else:
                sp_comp[xvar] = {ih_ind}
    # print(sp_comp)

    stoch_var = []

    for sp_c in sp_comp:
        row = []
        for r_ind in range(rxn_rows):
            prod = p_dict[r_ind][sp_c] if sp_c in p_dict[r_ind] else 0
            rect = r_dict[r_ind][sp_c] if sp_c in r_dict[r_ind] else 0
            row.append(prod - rect)
            if len(ks_dict[r_ind]) == 2:
                row.append(-row[-1])
        stoch_var.append(row)

    stoch_var = np.array(stoch_var)

    ksn_dict = {}
    concn = {}
    try:
        for r_ind in range(rxn_rows):
            ksn_dict[r_ind] = [0] * 2
            if len(r_dict[r_ind]) == 1:
                for xvar in r_dict[r_ind]:
                    if r_dict[r_ind][xvar] == 0:
                        ksn_dict[r_ind][0] = ks_dict[r_ind][0] \
                            * avgdrs_times_vol
                    elif r_dict[r_ind][xvar] == 1:
                        ksn_dict[r_ind][0] = ks_dict[r_ind][0]
                    elif r_dict[r_ind][xvar] == 2:
                        ksn_dict[r_ind][0] = factor_stch * ks_dict[r_ind][0] \
                            / avgdrs_times_vol
                    else:
                        ksn_dict[r_ind][0] = ks_dict[r_ind][0]
                    concn[xvar] = conc[xvar] * avgdrs_times_vol

            elif len(r_dict[r_ind]) == 2:
                ksn_dict[r_ind][0] = ks_dict[r_ind][0] / avgdrs_times_vol
                for xvar in r_dict[r_ind]:
                    concn[xvar] = conc[xvar] * avgdrs_times_vol

            if len(p_dict[r_ind]) == 1:
                for xvar in p_dict[r_ind]:
                    if p_dict[r_ind][xvar] == 0:
                        if len(ks_dict[r_ind]) == 2:
                            ksn_dict[r_ind][1] = ks_dict[r_ind][1] \
                                * avgdrs_times_vol
                        else:
                            ksn_dict[r_ind] = [ksn_dict[r_ind][0]]
                    elif p_dict[r_ind][xvar] == 1:
                        if len(ks_dict[r_ind]) == 2:
                            ksn_dict[r_ind][1] = ks_dict[r_ind][1]
                        else:
                            ksn_dict[r_ind] = [ksn_dict[r_ind][0]]
                    elif p_dict[r_ind][xvar] == 2:
                        if len(ks_dict[r_ind]) == 2:
                            ksn_dict[r_ind][1] = factor_stch \
                                * ks_dict[r_ind][1] / avgdrs_times_vol
                        else:
                            ksn_dict[r_ind] = [ksn_dict[r_ind][0]]
                    else:
                        if len(ks_dict[r_ind]) == 2:
                            ksn_dict[r_ind][1] = ks_dict[r_ind][1]
                        else:
                            ksn_dict[r_ind] = [ksn_dict[r_ind][0]]
                    concn[xvar] = conc[xvar] * avgdrs_times_vol
            elif len(p_dict[r_ind]) == 2:
                if len(ks_dict[r_ind]) == 2:
                    ksn_dict[r_ind][1] = ks_dict[r_ind][1] / avgdrs_times_vol
                else:
                    ksn_dict[r_ind] = [ksn_dict[r_ind][0]]
                for xvar in p_dict[r_ind]:
                    concn[xvar] = conc[xvar] * avgdrs_times_vol

        for xvar in conc:
            if xvar not in concn:
                concn[xvar] = conc[xvar]

        # for xvar in c_input:
            # concn[xvar] = c_input[xvar]

        # t_o = time.time()
        tvar = np.linspace(0, tend, int(tlen + 1))
        data = process_hub(
            tvar, sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
            v_volms, miter, logx, logy, del_coef, normalize, method,
            mix_plot, save, out_fname, plot_show, time_unit, vary,
            mult_proc, items, vary2, implicit, rfile, exp_data_file)
        # print(time.time()-t_o,"Process time")

        return data
    except Exception as error:
        message_upon_error.showinfo(
            "showinfo",
            "Check your topology files for missing species \
                in reaction and concentration tag : " + str(error))
