#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from scipy.integrate import odeint
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals
from BioSANS2020.propagation.deterministic.LNAapprox import lna_ss_jacobian, LNA_model_ss
from BioSANS2020.myglobal import mglobals as globals2
import numpy as np


def LNA_ode_model(z, t, Sp, ks_dict, r_dict, p_dict, V, molar=False):
    Spc = [s for s in Sp]
    conc = {Spc[x]: z[x] for x in range(len(Spc))}
    if not molar:
        D = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        D = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(V, D).reshape(len(z))
    for x in globals2.CON_BOUNDARY:
        ind = Spc.index(x)
        dxdt[ind] = 0
    return dxdt


def LNA_cov_model(A, B, C):
    return np.matmul(A, C) + np.matmul(C, np.transpose(A)) + B


def LNA_non_steady_state_old(
        conc, t, Sp, ks_dict, r_dict, p_dict, V, molar=True, rfile="", del_coef=10):
    get_globals(rfile)
    z = [conc[a] for a in Sp]
    zoftime = odeint(LNA_ode_model, z, t, args=(Sp, ks_dict, r_dict, p_dict, V, molar))

    half = []
    SiNew = []
    Sps = [s for s in Sp]
    lenSps = len(Sps)
    for i in range(lenSps):
        s = Sps[i]
        for j in range(i, lenSps):
            p = Sps[j]
            SiNew.append("cov(" + s + "," + p + ")")
            half.append(lenSps * i + j)

    C = np.zeros((len(z), len(z)))
    Cov = [[x for x in C.flatten()[half]]]
    dt = t[-1] - t[-2]
    tnew = [0]

    for stch_var in zoftime[1:]:
        ind = 0
        for sp in Sp:
            conc[sp] = stch_var[ind]
            ind = ind + 1
        AA = lna_ss_jacobian(LNA_model_ss, stch_var, Sp, V, ks_dict, r_dict, p_dict)
        f = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)
        BB = np.matmul(np.matmul(V, np.diag(f.flatten())), V.T)

        A = np.nan_to_num(AA)
        B = np.nan_to_num(BB)

        fx = LNA_cov_model(A, B, C)
        C = C + fx * dt
        Cov.append([x for x in C.flatten()[half]])

    return [np.array(Cov), SiNew, t]


def LNA_non_steady_state(conc, t, Sp, ks_dict, r_dict, p_dict, V,
                         molar=True, rfile="", del_coef=10):
    get_globals(rfile)
    z = [conc[a] for a in Sp]
    #zoftime = odeint(LNA_ode_model,z,t, args=(Sp,ks_dict,r_dict,p_dict,V,molar))

    half = []
    SiNew = []
    Sps = [s for s in Sp]
    lenSps = len(Sps)
    for i in range(lenSps):
        s = Sps[i]
        for j in range(i, lenSps):
            p = Sps[j]
            SiNew.append("cov(" + s + "," + p + ")")
            half.append(lenSps * i + j)

    C = np.zeros((len(z), len(z)))
    Cov = [[x for x in C.flatten()[half]]]
    dt = t[-1] - t[-2]
    tnew = [0]

    # for stch_var in zoftime[1:]:
    stch_var = np.array(z)
    while tnew[-1] < t[-1]:
        ind = 0
        for sp in Sp:
            conc[sp] = stch_var[ind]
            ind = ind + 1
        AA = lna_ss_jacobian(LNA_model_ss, stch_var, Sp, V, ks_dict, r_dict, p_dict)
        f = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)
        BB = np.matmul(np.matmul(V, np.diag(f.flatten())), V.T)

        A = np.nan_to_num(AA)
        B = np.nan_to_num(BB)

        fxCovr = LNA_cov_model(A, B, C)
        fxConc = LNA_ode_model(stch_var, t, Sp, ks_dict, r_dict, p_dict, V, molar)

        h1 = np.abs(fxCovr) + 1.0e-30
        dt = max(min(del_coef * np.min(1 / h1), dt), 1.0e-4)
        tnew.append(tnew[-1] + dt)
        stch_var = stch_var + fxConc * dt
        C = C + fxCovr * dt
        Cov.append([x for x in C.flatten()[half]])

    return [np.array(Cov), SiNew, tnew]


def LNA_non_steady_state2(conc, t, Sp, ks_dict, r_dict, p_dict, V,
                          molar=True, rfile="", del_coef=10):
    z = [conc[a] for a in Sp]
    Cov, slabels, tnew = LNA_non_steady_state(
        conc, t, Sp, ks_dict, r_dict, p_dict, V, molar, rfile, del_coef)
    zoftime = odeint(LNA_ode_model, z, tnew, args=(Sp, ks_dict, r_dict, p_dict, V, molar))

    SiNew = []
    Sps = [s for s in Sp]
    lenSps = len(Sps)
    FFdiv = []

    for x in slabels:
        SiNew.append(x.replace("cov", "FF"))

    for stch_var in zoftime:
        row = []
        for i in range(lenSps):
            s = Sps[i]
            for j in range(i, lenSps):
                p = Sps[j]
                Sij = stch_var[i] * stch_var[j]
                row.append(np.sqrt(Sij if Sij != 0 else 1))
        FFdiv.append(row)

    return [Cov / np.array(FFdiv), SiNew, tnew]
