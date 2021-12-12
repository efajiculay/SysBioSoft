#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals, \
    apply_rules
from BioSANS2020.propagation.deterministic.LNAapprox import \
    LNA_steady_state
from BioSANS2020.myglobal import mglobals as globals2


def euler_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                stch_var, d_time, del_coef, molar=False):

    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    fofx = np.matmul(stch_var, prop_flux)
    h1_div = np.abs(fofx) + 1.0e-30
    d_time = min(del_coef * np.min(1 / h1_div), d_time)
    dxdt = fofx * d_time
    if globals2.CON_BOUNDARY:
        spc = [x for x in sp_comp]
        for x in globals2.CON_BOUNDARY:
            ind = spc.index(x)
            dxdt[ind] = 0
    return [dxdt.reshape(len(sp_comp)), d_time]


def euler_int(
        t, sp_comp, ks_dict, sconc, r_dict, p_dict, stch_var,
        del_coef=10, LNAsolve=False, items=None, implicit=False,
        molar=False, rfile=""):
    get_globals(rfile)
    tnew = []
    d_time = t[-1] - t[-2]
    yconc = {x: sconc[x] for x in sconc}
    conc = {x: sconc[x] for x in sconc}
    apply_rules(conc, yconc)
    stch_var = [[conc[z] for z in sp_comp]]

    if not implicit:
        tnow = t[0]
        tnew.append(tnow)
        while tnow < t[-1]:
            m = np.nan_to_num(euler_model(
                sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, d_time, del_coef, molar))
            mm = m[0].reshape(1, len(m[0]))[0]
            tnow = tnow + m[1]
            ind = 0
            for sp in sp_comp:
                conc[sp] = stch_var[-1][ind] + mm[ind]
                ind = ind + 1
            apply_rules(conc, yconc)
            stch_var.append([conc[z] for z in sp_comp])
            tnew.append(tnow)
        if LNAsolve:
            return LNA_steady_state(
                t, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, items=items)
    else:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1
        C = [conc[z] for z in sp_comp]
        while tnew[-1] < t[-1]:
            m = np.nan_to_num(euler_model(
                sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, d_time, del_coef, molar))
            tnow = tnow + m[1]
            if tnow > t[tindex]:
                tnow = tnow - m[1]
                dt2 = t[tindex] - tnow
                m = np.nan_to_num(euler_model(
                    sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, dt2, del_coef, molar))
                mm = m[0].reshape(1, len(m[0]))[0]
                tnow = tnow + m[1]

                ind = 0
                for sp in sp_comp:
                    conc[sp] = C[ind] + mm[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                C = [conc[z] for z in sp_comp]
                for sp in sp_comp:
                    conc[sp] = max(0, conc[sp])
                stch_var.append(C)

                tnew.append(tnow)
                tindex = tindex + 1
            else:
                mm = m[0].reshape(1, len(m[0]))[0]
                ind = 0
                for sp in sp_comp:
                    conc[sp] = C[ind] + mm[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                C = [conc[z] for z in sp_comp]
                for sp in sp_comp:
                    conc[sp] = max(0, conc[sp])
    return (tnew, np.array(stch_var))


def euler2_model(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, molar=False):
    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(stch_var, prop_flux)
    if globals2.CON_BOUNDARY:
        spc = [x for x in sp_comp]
        for x in globals2.CON_BOUNDARY:
            ind = spc.index(x)
            dxdt[ind] = 0
    return dxdt.reshape(len(sp_comp))


def euler_wer_est(h, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, molar):
    k1 = h * euler2_model(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, molar)
    yconc = {}
    ind = 0
    for x in sp_comp:
        yconc[x] = conc[x] + k1[ind]
        ind = ind + 1
    err = 0.5 * h * euler2_model(sp_comp,
                                 ks_dict,
                                 yconc,
                                 r_dict,
                                 p_dict,
                                 stch_var,
                                 molar) - 0.5 * k1
    return [yconc, err]


def euler_help(htry, eps, yscal, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, molar):
    SAFETY = 0.9
    PGROW = -0.2
    PSHRNK = -0.25
    ERRCON = 1.89e-4
    errmax = 1000
    h = htry
    while errmax > 1:
        ytemp, yerr = euler_wer_est(
            h, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, molar)
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


def euler2_int(
        t, sp_comp, ks_dict, conc, r_dict, p_dict,
        stch_var, yscal=10, LNAsolve=False, items=None, implicit=False,
        molar=False, rfile=""):
    get_globals(rfile)
    tnew = []
    delt = t[-1] - t[-2]
    stch_var = [[conc[z] for z in sp_comp]]
    eps = 1.0e-4

    if not implicit:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1
        while tnew[-1] < t[-1]:
            conc_old = {x: conc[x] for x in conc}
            conc, d_time, delt, evar = euler_help(
                delt, eps, yscal, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, molar)
            tnow = tnow + d_time
            if tnow > t[tindex]:
                tnow = tnow - d_time
                d_time = t[tindex] - tnow
                conc, evar = euler_wer_est(
                    d_time, sp_comp, ks_dict, conc_old, r_dict, p_dict, stch_var, molar)
                delt = 5 * d_time
                tindex = tindex + 1
                tnow = tnow + d_time
            stch_var.append([conc[z] for z in sp_comp])
            tnew.append(tnow)

        if LNAsolve:
            return LNA_steady_state(
                t, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, items=items)
    else:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1
        while tnew[-1] < t[-1]:
            conc_old = {x: conc[x] for x in conc}
            conc, d_time, delt, evar = euler_help(
                delt, eps, yscal, sp_comp, ks_dict,
                conc, r_dict, p_dict, stch_var, molar)
            tnow = tnow + d_time
            if tnow > t[tindex]:
                tnow = tnow - d_time
                d_time = t[tindex] - tnow
                conc, evar = euler_wer_est(
                    d_time, sp_comp, ks_dict, conc_old, r_dict, p_dict,
                    stch_var, molar)
                delt = 5 * d_time
                stch_var.append([conc[z] for z in sp_comp])
                tnow = tnow + d_time
                tnew.append(tnow)
                tindex = tindex + 1
    return (tnew, np.array(stch_var))
