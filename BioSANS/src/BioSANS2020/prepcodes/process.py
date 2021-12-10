# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

# import warnings
from tkinter import messagebox as message_upon_error
# import re
# warnings.filterwarnings('ignore')

import time
import numpy as np

from BioSANS2020.prepcodes.processes_hub import process_hub
# from BioSANS2020.myglobal import mglobals as globals2
# from BioSANS2020.myglobal import proc_global
from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT


INF = np.inf
NaN = np.nan
inf = INF
nan = NaN


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

    A = 6.022e+23

    # globals2.MODIFIED = {}
    # globals2.PROP_MODIFIED = {}

    # for func in globals2.EXEC_FUNCTIONS:
    # exec(func, globals(), SBML_FUNCT_DICT)

    if in_molar == "molar":
        if method in ["ODE-1", "rk4-1", "rk4-1a", "Euler-1",
                      "Euler2-1", "CLE", "CLE2", "Tau-leaping",
                      "Tau-leaping2", "Sim-TauLeap", "Gillespie_"]:
            AV = A * v_volms
            fa = 2
        elif method in ["ODE-3", "rk4-3", "rk4-3a", "Euler-3", "Euler2-3"]:
            AV = v_volms
            fa = 1
        else:
            AV = 1
            fa = 1
    elif in_molar == "molecules":
        if method in ["ODE-2", "rk4-2", "rk4-2a", "Euler-2", "Euler2-2"]:
            AV = 1 / A * v_volms
            fa = 1 / 2
        elif method in ["ODE-3", "rk4-3", "rk4-3a", "Euler-3", "Euler2-3"]:
            AV = 1 / A
            fa = 1 / 2
        else:
            AV = 1
            fa = 1
    elif in_molar == "moles":
        if method in ["ODE-1", "rk4-1", "rk4-1a", "Euler-1", "Euler2-1",
                      "CLE", "CLE2", "Tau-leaping", "Tau-leaping2",
                      "Sim-TauLeap", "Gillespie_"]:
            AV = A
            fa = 2
        if method in ["ODE-2", "rk4-2", "rk4-2a", "Euler-2", "Euler2-2"]:
            AV = 1 / v_volms
            fa = 1
        else:
            AV = 1
            fa = 1
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
                    # globals2.MODIFIED[cvar[0].strip()] = [cc2,eval_dict(cc)]
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
    Sp = {}
    Rxn = len(rows)
    for ih in range(Rxn):
        r_dict[ih] = {}
        p_dict[ih] = {}
        col_row = rows[ih].split(":::::")
        row = col_row[0].strip().split(",")
        # if len(col_row)>1:
        # krow = col_row[1].strip().split(":::")
        # if len(krow)==2:
        # cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
        # cc3 = krow[1].split(":")[0].replace("lambda","").split(",")
        # globals2.PROP_MODIFIED["Prop_"+str(ih)] = [
        #     (cc2,eval_dict(krow[0])),(cc3,eval_dict(krow[1]))]
        # else:
        # cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
        # globals2.PROP_MODIFIED["Prop_"+str(ih)] = \
        #     [(cc2,eval_dict(krow[0]))]

        if len(row) == 3:
            ks_dict[ih] = [
                tofloat(row[1], locals()),
                tofloat(row[2], locals())]
        else:
            ks_dict[ih] = [tofloat(row[1], locals())]
        var = row[0].split("<=>")
        if len(var) == 1:
            var = row[0].split("=>")

        sp = var[0]
        s = sp.strip().split()
        if len(s) > 1:
            last = 1
            for xvar in s:
                if not is_number(xvar) and xvar != "+":
                    r_dict[ih][xvar] = last
                    last = 1
                    if xvar in Sp:
                        Sp[xvar].add(ih)
                    else:
                        Sp[xvar] = {ih}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = s[0]
            r_dict[ih][xvar] = 1
            if xvar in Sp:
                Sp[xvar].add(ih)
            else:
                Sp[xvar] = {ih}

        sp = var[1]
        s = sp.strip().split()
        if len(s) > 1:
            last = 1
            for xvar in s:
                if (not is_number(xvar) or xvar.lower() == "e") \
                        and xvar != "+":
                    p_dict[ih][xvar] = last
                    last = 1
                    if xvar in Sp:
                        Sp[xvar].add(ih)
                    else:
                        Sp[xvar] = {ih}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = s[0]
            p_dict[ih][xvar] = 1
            if xvar in Sp:
                Sp[xvar].add(ih)
            else:
                Sp[xvar] = {ih}
    # print(Sp)

    V = []

    for sp in Sp:
        row = []
        for r in range(Rxn):
            prod = p_dict[r][sp] if sp in p_dict[r] else 0
            rect = r_dict[r][sp] if sp in r_dict[r] else 0
            row.append(prod - rect)
            if len(ks_dict[r]) == 2:
                row.append(-row[-1])
        V.append(row)

    V = np.array(V)

    Ksn = {}
    concn = {}
    try:
        for r in range(Rxn):
            Ksn[r] = [0] * 2
            if len(r_dict[r]) == 1:
                for xvar in r_dict[r]:
                    if r_dict[r][xvar] == 0:
                        Ksn[r][0] = ks_dict[r][0] * AV
                    elif r_dict[r][xvar] == 1:
                        Ksn[r][0] = ks_dict[r][0]
                    elif r_dict[r][xvar] == 2:
                        Ksn[r][0] = fa * ks_dict[r][0] / AV
                    else:
                        Ksn[r][0] = ks_dict[r][0]
                    concn[xvar] = conc[xvar] * AV

            elif len(r_dict[r]) == 2:
                Ksn[r][0] = ks_dict[r][0] / AV
                for xvar in r_dict[r]:
                    concn[xvar] = conc[xvar] * AV

            if len(p_dict[r]) == 1:
                for xvar in p_dict[r]:
                    if p_dict[r][xvar] == 0:
                        if len(ks_dict[r]) == 2:
                            Ksn[r][1] = ks_dict[r][1] * AV
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    elif p_dict[r][xvar] == 1:
                        if len(ks_dict[r]) == 2:
                            Ksn[r][1] = ks_dict[r][1]
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    elif p_dict[r][xvar] == 2:
                        if len(ks_dict[r]) == 2:
                            Ksn[r][1] = fa * ks_dict[r][1] / AV
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    else:
                        if len(ks_dict[r]) == 2:
                            Ksn[r][1] = ks_dict[r][1]
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    concn[xvar] = conc[xvar] * AV
            elif len(p_dict[r]) == 2:
                if len(ks_dict[r]) == 2:
                    Ksn[r][1] = ks_dict[r][1] / AV
                else:
                    Ksn[r] = [Ksn[r][0]]
                for xvar in p_dict[r]:
                    concn[xvar] = conc[xvar] * AV

        for xvar in conc:
            if xvar not in concn:
                concn[xvar] = conc[xvar]

        # for xvar in c_input:
            # concn[xvar] = c_input[xvar]

        t_o = time.time()
        t = np.linspace(0, tend, int(tlen + 1))
        data = process_hub(
            t, Sp, Ksn, concn, r_dict, p_dict, V, v_volms, miter,
            logx, logy, del_coef, normalize, method, mix_plot,
            save, out_fname, plot_show, time_unit, vary,
            mult_proc, items, vary2, implicit, rfile, exp_data_file)
        # print(time.time()-t_o,"Process time")

        return data
    except Exception as e:
        message_upon_error.showinfo(
            "showinfo",
            "Check your topology files for missing species \
                in reaction and concentration tag : " + str(e))
