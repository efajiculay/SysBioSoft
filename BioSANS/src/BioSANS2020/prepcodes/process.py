#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

#!/usr/bin/env python
# coding: utf-8

#import warnings
import re
# warnings.filterwarnings('ignore')

import numpy as np
import time
from BioSANS2020.prepcodes.processes_hub import *
#from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.math_functs.sbmlMath import *

from tkinter import messagebox as message_upon_error

INF = np.inf
NaN = np.nan
inf = INF
nan = NaN


def tofloat(val):
    try:
        return float(eval(val))
    except:
        return float(val)


def isNumber(x):
    try:
        tofloat(x)
        return True
    except:
        return False


def process(
        rfile="Reactions",
        miter=1,
        inMolar=True,
        Vm=1.0e-20,
        tn=1,
        delX=10,
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
        expDataFile=None,
        Cinput={}
):

    A = 6.022e+23

    #globals2.modified = {}
    #globals2.PropModified = {}

    # for func in globals2.execFunctions:
    # exec(func,globals())

    if inMolar == "molar":
        if method in ["ODE-1", "rk4-1", "rk4-1a", "Euler-1", "Euler2-1", "CLE", "CLE2", "Tau-leaping", "Tau-leaping2", "Sim-TauLeap", "Gillespie_"]:
            AV = A*Vm
            fa = 2
        elif method in ["ODE-3", "rk4-3", "rk4-3a", "Euler-3", "Euler2-3"]:
            AV = Vm
            fa = 1
        else:
            AV = 1
            fa = 1
    elif inMolar == "molecules":
        if method in ["ODE-2", "rk4-2", "rk4-2a", "Euler-2", "Euler2-2"]:
            AV = 1/A*Vm
            fa = 1/2
        elif method in ["ODE-3", "rk4-3", "rk4-3a", "Euler-3", "Euler2-3"]:
            AV = 1/A
            fa = 1/2
        else:
            AV = 1
            fa = 1
    elif inMolar == "moles":
        if method in ["ODE-1", "rk4-1", "rk4-1a", "Euler-1", "Euler2-1", "CLE", "CLE2", "Tau-leaping", "Tau-leaping2", "Sim-TauLeap", "Gillespie_"]:
            AV = A
            fa = 2
        if method in ["ODE-2", "rk4-2", "rk4-2a", "Euler-2", "Euler2-2"]:
            AV = 1/Vm
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
                    exec(row.strip(), globals())
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
                    conc[cvar[0].strip()] = tofloat(cvar[1])
                    # if len(cvar)>=3:
                    #cc = ",".join(cvar[2:])
                    #cc2 = cc.split(":")[0].replace("lambda","").split(",")
                    #globals2.modified[cvar[0].strip()] = [cc2,eval(cc)]
            elif row[0] == "#":
                last = "#"
                #gg = row.split(",")[1:]
                # try:
                # for x in gg:
                #xx = x.split("=")
                #globals2.settings[xx[0].strip()] = xx[1].strip()
                # except:
                # pass
            elif row[0] == "@":
                last = "@"
            elif row.strip() == "Function_Definitions:":
                last = "Function_Definitions"

        file.close()

    Ks = {}
    Rr = {}
    Rp = {}
    Sp = {}
    Rxn = len(rows)
    for ih in range(Rxn):
        Rr[ih] = {}
        Rp[ih] = {}
        col_row = rows[ih].split(":::::")
        row = col_row[0].strip().split(",")
        # if len(col_row)>1:
        #krow = col_row[1].strip().split(":::")
        # if len(krow)==2:
        #cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
        #cc3 = krow[1].split(":")[0].replace("lambda","").split(",")
        #globals2.PropModified["Prop_"+str(ih)] = [(cc2,eval(krow[0])),(cc3,eval(krow[1]))]
        # else:
        #cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
        #globals2.PropModified["Prop_"+str(ih)] = [(cc2,eval(krow[0]))]

        if len(row) == 3:
            Ks[ih] = [tofloat(row[1]), tofloat(row[2])]
        else:
            Ks[ih] = [tofloat(row[1])]
        var = row[0].split("<=>")
        if len(var) == 1:
            var = row[0].split("=>")

        sp = var[0]
        s = sp.strip().split()
        if len(s) > 1:
            last = 1
            for x in s:
                if not isNumber(x) and x != "+":
                    Rr[ih][x] = last
                    last = 1
                    if x in Sp:
                        Sp[x].add(ih)
                    else:
                        Sp[x] = {ih}
                elif isNumber(x):
                    last = tofloat(x)
        else:
            x = s[0]
            Rr[ih][x] = 1
            if x in Sp:
                Sp[x].add(ih)
            else:
                Sp[x] = {ih}

        sp = var[1]
        s = sp.strip().split()
        if len(s) > 1:
            last = 1
            for x in s:
                if (not isNumber(x) or x.lower() == "e") and x != "+":
                    Rp[ih][x] = last
                    last = 1
                    if x in Sp:
                        Sp[x].add(ih)
                    else:
                        Sp[x] = {ih}
                elif isNumber(x):
                    last = tofloat(x)
        else:
            x = s[0]
            Rp[ih][x] = 1
            if x in Sp:
                Sp[x].add(ih)
            else:
                Sp[x] = {ih}
    # print(Sp)

    V = []

    for sp in Sp:
        row = []
        for r in range(Rxn):
            prod = Rp[r][sp] if sp in Rp[r] else 0
            rect = Rr[r][sp] if sp in Rr[r] else 0
            row.append(prod-rect)
            if len(Ks[r]) == 2:
                row.append(-row[-1])
        V.append(row)

    V = np.array(V)

    Ksn = {}
    concn = {}
    try:
        for r in range(Rxn):
            Ksn[r] = [0]*2
            if len(Rr[r]) == 1:
                for x in Rr[r]:
                    if Rr[r][x] == 0:
                        Ksn[r][0] = Ks[r][0]*AV
                    elif Rr[r][x] == 1:
                        Ksn[r][0] = Ks[r][0]
                    elif Rr[r][x] == 2:
                        Ksn[r][0] = fa*Ks[r][0]/AV
                    else:
                        Ksn[r][0] = Ks[r][0]
                    concn[x] = conc[x]*AV

            elif len(Rr[r]) == 2:
                Ksn[r][0] = Ks[r][0]/AV
                for x in Rr[r]:
                    concn[x] = conc[x]*AV

            if len(Rp[r]) == 1:
                for x in Rp[r]:
                    if Rp[r][x] == 0:
                        if len(Ks[r]) == 2:
                            Ksn[r][1] = Ks[r][1]*AV
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    elif Rp[r][x] == 1:
                        if len(Ks[r]) == 2:
                            Ksn[r][1] = Ks[r][1]
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    elif Rp[r][x] == 2:
                        if len(Ks[r]) == 2:
                            Ksn[r][1] = fa*Ks[r][1]/AV
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    else:
                        if len(Ks[r]) == 2:
                            Ksn[r][1] = Ks[r][1]
                        else:
                            Ksn[r] = [Ksn[r][0]]
                    concn[x] = conc[x]*AV
            elif len(Rp[r]) == 2:
                if len(Ks[r]) == 2:
                    Ksn[r][1] = Ks[r][1]/AV
                else:
                    Ksn[r] = [Ksn[r][0]]
                for x in Rp[r]:
                    concn[x] = conc[x]*AV

        for x in conc:
            if x not in concn:
                concn[x] = conc[x]

        # for x in Cinput:
            #concn[x] = Cinput[x]

        t_o = time.time()
        t = np.linspace(0, tn, int(tlen+1))
        data = process_hub(t, Sp, Ksn, concn, Rr, Rp, V, Vm, miter, logx, logy, delX, normalize, method, mix_plot,
                           save, out_fname, plot_show, time_unit, vary, mult_proc, items, vary2, implicit, rfile, expDataFile)
        #print(time.time()-t_o,"Process time")

        return data
    except Exception as e:
        message_upon_error.showinfo(
            "showinfo", "Check your topology files for missing species in reaction and concentration tag : "+str(e))
