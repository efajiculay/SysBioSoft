#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import *
from BioSANS2020.propagation.recalculate_globals import *
from BioSANS2020.propagation.deterministic.LNAapprox import *
from BioSANS2020.myglobal import mglobals as globals2


def Euler_model(Sp, ks_dict, conc, r_dict, p_dict, V, dt, del_coef, molar=False):

    if not molar:
        D = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        D = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    fx = np.matmul(V, D)
    h1 = np.abs(fx) + 1.0e-30
    dt = min(del_coef * np.min(1 / h1), dt)
    dxdt = fx * dt
    if globals2.CON_BOUNDARY:
        Spc = [x for x in Sp]
        for x in globals2.CON_BOUNDARY:
            ind = Spc.index(x)
            dxdt[ind] = 0
    return [dxdt.reshape(len(Sp)), dt]


def Euler_int(t, Sp, ks_dict, sconc, r_dict, p_dict, V, del_coef=10, LNAsolve=False,
              items=None, implicit=False, molar=False, rfile=""):
    get_globals(rfile)
    tnew = []
    dt = t[-1] - t[-2]
    yconc = {x: sconc[x] for x in sconc}
    conc = {x: sconc[x] for x in sconc}
    apply_rules(conc, yconc)
    stch_var = [[conc[z] for z in Sp]]

    if not implicit:
        tnow = t[0]
        tnew.append(tnow)
        while tnow < t[-1]:
            m = np.nan_to_num(Euler_model(
                Sp, ks_dict, conc, r_dict, p_dict, V, dt, del_coef, molar))
            mm = m[0].reshape(1, len(m[0]))[0]
            tnow = tnow + m[1]
            ind = 0
            for sp in Sp:
                conc[sp] = stch_var[-1][ind] + mm[ind]
                ind = ind + 1
            apply_rules(conc, yconc)
            stch_var.append([conc[z] for z in Sp])
            tnew.append(tnow)
        if LNAsolve:
            return LNA_steady_state(t, Sp, ks_dict, conc, r_dict, p_dict, V, items=items)
    else:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1
        C = [conc[z] for z in Sp]
        while tnew[-1] < t[-1]:
            m = np.nan_to_num(Euler_model(
                Sp, ks_dict, conc, r_dict, p_dict, V, dt, del_coef, molar))
            tnow = tnow + m[1]
            if tnow > t[tindex]:
                tnow = tnow - m[1]
                dt2 = t[tindex] - tnow
                m = np.nan_to_num(Euler_model(
                    Sp, ks_dict, conc, r_dict, p_dict, V, dt2, del_coef, molar))
                mm = m[0].reshape(1, len(m[0]))[0]
                tnow = tnow + m[1]

                ind = 0
                for sp in Sp:
                    conc[sp] = C[ind] + mm[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                C = [conc[z] for z in Sp]
                for sp in Sp:
                    conc[sp] = max(0, conc[sp])
                stch_var.append(C)

                tnew.append(tnow)
                tindex = tindex + 1
            else:
                mm = m[0].reshape(1, len(m[0]))[0]
                ind = 0
                for sp in Sp:
                    conc[sp] = C[ind] + mm[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                C = [conc[z] for z in Sp]
                for sp in Sp:
                    conc[sp] = max(0, conc[sp])
    return (tnew, np.array(stch_var))


def Euler2_model(Sp, ks_dict, conc, r_dict, p_dict, V, molar=False):
    if not molar:
        D = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        D = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(V, D)
    if globals2.CON_BOUNDARY:
        Spc = [x for x in Sp]
        for x in globals2.CON_BOUNDARY:
            ind = Spc.index(x)
            dxdt[ind] = 0
    return dxdt.reshape(len(Sp))


def EulerWerEst(h, Sp, ks_dict, conc, r_dict, p_dict, V, molar):
    k1 = h * Euler2_model(Sp, ks_dict, conc, r_dict, p_dict, V, molar)
    yconc = {}
    ind = 0
    for x in Sp:
        yconc[x] = conc[x] + k1[ind]
        ind = ind + 1
    err = 0.5 * h * Euler2_model(Sp, ks_dict, yconc, r_dict, p_dict, V, molar) - 0.5 * k1
    return [yconc, err]


def EulerHelp(htry, eps, yscal, Sp, ks_dict, conc, r_dict, p_dict, V, molar):
    SAFETY = 0.9
    PGROW = -0.2
    PSHRNK = -0.25
    ERRCON = 1.89e-4
    errmax = 1000
    h = htry
    while errmax > 1:
        ytemp, yerr = EulerWerEst(h, Sp, ks_dict, conc, r_dict, p_dict, V, molar)
        errmax = max(0, np.max(yerr / yscal)) / eps
        if errmax > 1:
            h = SAFETY * h * (errmax**PSHRNK)
            if h < 0.1 * h:
                h = 0.1 * h
        else:
            if errmax > ERRCON:
                hnext = SAFETY * h * (errmax**PGROW)
            else:
                hnext = 5.0 * h
            hdid = h
            y = ytemp

    return [y, hdid, hnext, yerr]


def Euler2_int(t, Sp, ks_dict, conc, r_dict, p_dict, V, yscal=10, LNAsolve=False,
               items=None, implicit=False, molar=False, rfile=""):
    get_globals(rfile)
    tnew = []
    delt = t[-1] - t[-2]
    stch_var = [[conc[z] for z in Sp]]
    eps = 1.0e-4

    if not implicit:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1
        while tnew[-1] < t[-1]:
            conc_old = {x: conc[x] for x in conc}
            conc, dt, delt, e = EulerHelp(
                delt, eps, yscal, Sp, ks_dict, conc, r_dict, p_dict, V, molar)
            tnow = tnow + dt
            if tnow > t[tindex]:
                tnow = tnow - dt
                dt = t[tindex] - tnow
                conc, e = EulerWerEst(dt, Sp, ks_dict, conc_old, r_dict, p_dict, V, molar)
                delt = 5 * dt
                tindex = tindex + 1
                tnow = tnow + dt
            stch_var.append([conc[z] for z in Sp])
            tnew.append(tnow)

        if LNAsolve:
            return LNA_steady_state(t, Sp, ks_dict, conc, r_dict, p_dict, V, items=items)
    else:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1
        while tnew[-1] < t[-1]:
            conc_old = {x: conc[x] for x in conc}
            conc, dt, delt, e = EulerHelp(
                delt, eps, yscal, Sp, ks_dict, conc, r_dict, p_dict, V, molar)
            tnow = tnow + dt
            if tnow > t[tindex]:
                tnow = tnow - dt
                dt = t[tindex] - tnow
                conc, e = EulerWerEst(dt, Sp, ks_dict, conc_old, r_dict, p_dict, V, molar)
                delt = 5 * dt
                stch_var.append([conc[z] for z in Sp])
                tnow = tnow + dt
                tnew.append(tnow)
                tindex = tindex + 1
    return (tnew, np.array(stch_var))
