"""

                        This is the process module

This module reads BioSANS topology file, grab the components or species,
rate  constants, stoichiometric  matrix,  propensity   vector, algebraic
rules, conditional statements, and othe types of definitions into a dic-
tionary.  This module  calls the process_hub  module that distribute the
tasks to other modules.

The functions in this module are as follows;

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
    return eval(to_eval, loc_dict, SBML_FUNCT_DICT)


def tofloat(val, loc_dict):
    try:
        return float(eval_dict(val, loc_dict))
    except BaseException:
        return float(val)


def is_number(xvar):
    try:
        float(xvar)
        return True
    except ValueError:
        return False


def process(
        rfile="Reactions",
        miter=1,
        in_molar=True,
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

    avogadros_num = 6.022e+23

    # globals2.MODIFIED = {}
    # globals2.PROP_MODIFIED = {}

    # for func in globals2.EXEC_FUNCTIONS:
    #     exec(func, globals(), SBML_FUNCT_DICT)

    if in_molar == "molar":
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
    elif in_molar == "molecules":
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
    elif in_molar == "moles":
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
