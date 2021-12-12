#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals
from BioSANS2020.myglobal import mglobals as globals2


def cle_model(Sp, ks_dict, conc, r_dict, p_dict, V, dt, del_coef, reg=False):
    D = propensity_vec(ks_dict, conc, r_dict, p_dict)
    G = np.sqrt(D)
    nlen = np.random.randn(len(D)).reshape(len(D), 1)
    fx = np.matmul(V, D)
    gx = np.matmul(V, G * nlen)
    if not reg:
        h1 = np.abs(D) + 1.0e-30  # np.abs(fx)+1.0e-30
        h2 = np.abs(G) + 1.0e-30  # np.abs(gx)+1.0e-30
        dt = min(del_coef * min(np.min(1 / h1), np.min(1 / h2)), dt)
    sqdt = np.sqrt(dt)

    ind = 0
    for sp in Sp:
        key = sp.strip().split("_")[0]
        if key in reserve_events_words:
            gx[ind] = 0
        ind = ind + 1

    fd = fx * dt + gx * sqdt
    return [fd.reshape(len(Sp)), dt]


def cle_calculate(t, Sp, ks_dict, sconc, r_dict, p_dict, V, del_coef=10,
                  rr=1, implicit=False, rfile=""):
    get_globals(rfile)
    np.random.seed(int(rr * 100))
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
            m = np.nan_to_num(cle_model(Sp, ks_dict, conc, r_dict, p_dict, V, dt, del_coef))
            mm = m[0].reshape(1, len(m[0]))[0]
            tnow = tnow + m[1]
            ind = 0
            for sp in Sp:
                conc[sp] = stch_var[-1][ind] + mm[ind]
                ind = ind + 1
            apply_rules(conc, yconc)
            stch_var.append([conc[z] for z in Sp])
            for sp in Sp:
                conc[sp] = max(0, conc[sp])
            tnew.append(tnow)
    else:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1
        C = [conc[z] for z in Sp]
        while tnew[-1] < t[-1]:
            m = np.nan_to_num(cle_model(Sp, ks_dict, conc, r_dict, p_dict, V, dt, del_coef))
            tnow = tnow + m[1]
            if tnow > t[tindex]:
                tnow = tnow - m[1]
                dt2 = t[tindex] - tnow
                m = np.nan_to_num(
                    cle_model(Sp, ks_dict, conc, r_dict, p_dict, V, dt2, del_coef, True))
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


def cle2_calculate(t, Sp, ks_dict, sconc, r_dict, p_dict, V, del_coef=1, rr=1, rfile=""):
    get_globals(rfile)
    np.random.seed(int(rr * 100))
    div = max(1, int(1 / del_coef))
    dt = (t[-1] - t[-2]) / div
    yconc = {x: sconc[x] for x in sconc}
    conc = {x: sconc[x] for x in sconc}
    apply_rules(conc, yconc)
    stch_var = [[conc[z] for z in Sp]]
    tnow = t[0]
    tnew = [tnow]
    dv = 0
    C = stch_var[-1]
    while abs(tnow - t[-1]) > 1.0e-10:
        m = np.nan_to_num(cle_model(Sp, ks_dict, conc, r_dict, p_dict, V, dt, 1, True))
        mm = m[0].reshape(1, len(m[0]))[0]
        tnow = tnow + m[1]
        ind = 0
        for sp in Sp:
            conc[sp] = C[ind] + mm[ind]
            ind = ind + 1
        apply_rules(conc, yconc)
        C = [conc[z] for z in Sp]
        if dv == div - 1:
            stch_var.append(C)
            tnew.append(tnow)
            dv = 0
        else:
            dv = dv + 1
        for sp in Sp:
            conc[sp] = max(0, conc[sp])

    return (tnew, np.array(stch_var))
