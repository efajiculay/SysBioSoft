#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import *
from BioSANS2020.propagation.recalculate_globals import *
from BioSANS2020.myglobal import mglobals as globals2


def rk4_model(Sp, Ks, conc, Rr, Rp, V, t, molar=False):
    if not molar:
        D = propensity_vec(Ks, conc, Rr, Rp)
    else:
        D = propensity_vec_molar(Ks, conc, Rr, Rp)
    dxdt = np.matmul(V, D).reshape(len(Sp))
    if globals2.conBoundary:
        Spc = [sp for sp in Sp]
        for x in globals2.conBoundary:
            ind = Spc.index(x)
            dxdt[ind] = 0
    return dxdt


def runge_kutta_forth(Sp, Ks, conc, Rr, Rp, V, t, delt, nSp, molar=False):

    yconc = {x: conc[x] for x in nSp}

    da_dt = rk4_model(Sp, Ks, conc, Rr, Rp, V, t, molar)
    newA1 = da_dt * delt
    ind = 0
    for sp in Sp:
        yconc[sp] = conc[sp] + 0.5 * newA1[ind]
        ind = ind + 1

    da_dt = rk4_model(Sp, Ks, yconc, Rr, Rp, V, t + 0.5 * delt, molar)
    newA2 = da_dt * delt
    ind = 0
    for sp in Sp:
        yconc[sp] = conc[sp] + 0.5 * newA2[ind]
        ind = ind + 1

    da_dt = rk4_model(Sp, Ks, yconc, Rr, Rp, V, t + 0.5 * delt, molar)
    newA3 = da_dt * delt
    ind = 0
    for sp in Sp:
        yconc[sp] = conc[sp] + 0.5 * newA3[ind]
        ind = ind + 1

    da_dt = rk4_model(Sp, Ks, yconc, Rr, Rp, V, t + delt, molar)
    newA4 = da_dt * delt

    NEW = (newA1 + 2.0 * newA2 + 2.0 * newA3 + newA4) / 6.0

    newA = yconc
    ind = 0
    for sp in Sp:
        newA[sp] = conc[sp] + NEW[ind]
        ind = ind + 1

    return [newA, t + delt]


def rungek4_int(conc, time, Sp, Ks, Rr, Rp, V, molar=False, delx=1, rfile=""):
    get_globals(rfile)
    tn = time[-1]
    div = max(1, int(1 / delx))
    delt = (time[-1] - time[-2]) / div
    yconc = {x: conc[x] for x in conc}
    slabels = [a for a in Sp]
    apply_rules(conc, yconc, [0], [conc[a] for a in Sp], slabels)
    z = [conc[a] for a in Sp]
    Z = [z]
    t = 0.0
    Tc = [t]
    nSp = [x for x in conc if x not in Sp]
    dv = 0
    Tc2 = [t]
    Z2 = [z]

    while abs(t - tn) > 1.0e-10:
        conc, t = runge_kutta_forth(
            Sp, Ks, conc, Rr, Rp, V, t, delt, nSp, molar)
        # conc, yerr = rkck(delt,Sp,Ks,conc,Rr,Rp,V,t,nSp,molar) #higher order runge-kutta
        # t = t + delt                                           #higher order
        # runge-kutta
        apply_rules(conc, yconc, Tc2, Z2, slabels)
        Z2.append([conc[a] for a in Sp])
        Tc2.append(t)
        if dv == div - 1:
            Z.append([conc[a] for a in Sp])
            Tc.append(t)
            dv = 0
        else:
            dv = dv + 1
    Z = np.array(Z)
    return [Tc, Z]


def rkck(h, Sp, Ks, conc, Rr, Rp, V, t, nSp, molar):

    A2 = 0.2
    A3 = 0.3
    A4 = 0.6
    A5 = 1.0
    A6 = 0.875

    B21 = 0.2
    B31 = 3.0 / 40.0
    B32 = 9.0 / 40.0

    B41 = 0.3
    B42 = -0.9
    B43 = 1.2

    B51 = -11.9 / 54.0
    B52 = 2.5
    B53 = -70.0 / 27.0
    B54 = 35.0 / 27.0

    B61 = 1631.0 / 55296.0
    B62 = 175.0 / 512.0
    B63 = 575.0 / 13824.0
    B64 = 44275.0 / 110592.0
    B65 = 253.0 / 4096.0

    C1 = 37.0 / 378.0
    C3 = 250.0 / 621.0
    C4 = 125.0 / 594.0
    C6 = 512.0 / 1771.0

    DC1 = C1 - 2825.0 / 27648.0
    DC3 = C3 - 18575.0 / 48384.0
    DC4 = C4 - 13525.0 / 55296.0
    DC5 = -277.0 / 14336.0
    DC6 = C6 - 0.25

    dydx = rk4_model(Sp, Ks, conc, Rr, Rp, V, t, molar)

    ytemp = {x: conc[x] for x in nSp}
    ind = 0
    for sp in Sp:
        ytemp[sp] = conc[sp] + B21 * h * dydx[ind]
        ind = ind + 1
    ak2 = rk4_model(Sp, Ks, ytemp, Rr, Rp, V, t + A2 * h, molar)

    ind = 0
    for sp in Sp:
        ytemp[sp] = conc[sp] + h * (B31 * dydx[ind] + B32 * ak2[ind])
        ind = ind + 1
    ak3 = rk4_model(Sp, Ks, ytemp, Rr, Rp, V, t + A3 * h, molar)

    ind = 0
    for sp in Sp:
        ytemp[sp] = conc[sp] + h * \
            (B41 * dydx[ind] + B42 * ak2[ind] + B43 * ak3[ind])
        ind = ind + 1
    ak4 = rk4_model(Sp, Ks, ytemp, Rr, Rp, V, t + A4 * h, molar)

    ind = 0
    for sp in Sp:
        ytemp[sp] = conc[sp] + h * \
            (B51 * dydx[ind] + B52 * ak2[ind] +
             B53 * ak3[ind] + B54 * ak4[ind])
        ind = ind + 1
    ak5 = rk4_model(Sp, Ks, ytemp, Rr, Rp, V, t + A5 * h, molar)

    ind = 0
    for sp in Sp:
        ytemp[sp] = conc[sp] + h * (B61 * dydx[ind] + B62 *
                                    ak2[ind] + B63 * ak3[ind] + B64 * ak4[ind] + B65 * ak5[ind])
        ind = ind + 1
    ak6 = rk4_model(Sp, Ks, ytemp, Rr, Rp, V, t + A6 * h, molar)

    yout = ytemp
    yerr = []
    ind = 0
    for sp in Sp:
        yout[sp] = conc[sp] + h * \
            (C1 * dydx[ind] + C3 * ak3[ind] + C4 * ak4[ind] + C6 * ak6[ind])
        yerr.append(h * (DC1 * dydx[ind] + DC3 * ak3[ind] + DC4 *
                         ak4[ind] + DC5 * ak5[ind] + DC6 * ak6[ind]))
        ind = ind + 1

    return [yout, yerr]


def rkqs(htry, eps, yscal, Sp, Ks, conc, Rr, Rp, V, t, nSp, molar):
    SAFETY = 0.9
    PGROW = -0.2
    PSHRNK = -0.25
    ERRCON = 1.89e-4
    errmax = 1000
    h = htry
    while errmax > 1:
        ytemp, yerr = rkck(h, Sp, Ks, conc, Rr, Rp, V, t, nSp, molar)
        errmax = max(0, np.max(np.array(yerr) / np.array(yscal))) / eps
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
            conc = ytemp

    return [conc, hdid, hnext, yerr]


def rungek4a_int(t, Sp, Ks, conc, Rr, Rp, V, yscal=10,
                 molar=False, implicit=False, rfile=""):
    get_globals(rfile)
    tnew = []
    delt = t[-1] - t[-2]
    eps = 1.0e-8
    nSp = [x for x in conc if x not in Sp]
    yconc = {x: conc[x] for x in Sp}
    slabels = [a for a in Sp]
    apply_rules(conc, yconc, [0], [conc[a] for a in Sp], slabels)
    y = [conc[a] for a in Sp]
    stch_var = [y]
    if not implicit:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1

        while abs((tnew[-1] - t[-1]) / t[-1]) > 1.0e-5:

            y_old = conc
            conc, dt, delt, e = rkqs(
                delt, eps, yscal, Sp, Ks, conc, Rr, Rp, V, tnow, nSp, molar)
            tnow = tnow + dt
            apply_rules(conc, yconc)
            if tnow > t[tindex]:
                tnow = tnow - dt
                dt = t[tindex] - tnow
                conc, tnow = runge_kutta_forth(
                    Sp, Ks, y_old, Rr, Rp, V, tnow, dt, nSp, molar)
                delt = 5 * dt
                apply_rules(conc, yconc)
                tindex = tindex + 1
            stch_var.append([conc[a] for a in Sp])
            tnew.append(tnow)
    else:
        tnow = t[0]
        tnew.append(tnow)
        tindex = 1

        while abs((tnew[-1] - t[-1]) / t[-1]) > 1.0e-5:

            y_old = conc
            conc, dt, delt, e = rkqs(
                delt, eps, yscal, Sp, Ks, conc, Rr, Rp, V, tnow, nSp, molar)
            tnow = tnow + dt
            apply_rules(conc, yconc)
            if tnow > t[tindex]:
                tnow = tnow - dt
                dt = t[tindex] - tnow
                conc, tnow = runge_kutta_forth(
                    Sp, Ks, y_old, Rr, Rp, V, tnow, dt, nSp, molar)
                delt = 5 * dt
                apply_rules(conc, yconc)
                stch_var.append([conc[a] for a in Sp])
                tnew.append(tnow)
                tindex = tindex + 1
    return [tnew, np.array(stch_var)]
