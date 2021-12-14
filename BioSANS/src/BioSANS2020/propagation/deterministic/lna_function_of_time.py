# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from scipy.integrate import odeint
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals
from BioSANS2020.propagation.deterministic.lna_approx \
    import lna_ss_jacobian, lna_model_ss
from BioSANS2020.myglobal import mglobals as globals2



def lna_ode_model(zlist, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var, molar=False):
    spc = list(sp_comp.keys())
    conc = {spc[xvar]: zlist[xvar] for xvar in range(len(spc))}
    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(stch_var, prop_flux).reshape(len(zlist))
    for xvar in globals2.CON_BOUNDARY:
        ind = spc.index(xvar)
        dxdt[ind] = 0
    return dxdt


def lna_cov_model(a_jac, b_diff, cov):
    return np.matmul(a_jac, cov) + np.matmul(cov, np.transpose(a_jac)) + b_diff


def lna_non_steady_state_old(
        conc, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var,
        molar=True, rfile="", del_coef=10):
    get_globals(rfile)
    zlist = [conc[a] for a in sp_comp]
    zoftime = odeint(
        lna_ode_model, zlist, t_var,
        args=(sp_comp, ks_dict, r_dict, p_dict, stch_var, molar))

    half = []
    si_new = []
    sps = [svar for svar in sp_comp]
    len_sps = len(sps)
    for i in range(len_sps):
        svar = sps[i]
        for j in range(i, len_sps):
            pvar = sps[j]
            si_new.append("cov(" + svar + "," + pvar + ")")
            half.append(len_sps * i + j)

    cvar = np.zeros((len(zlist), len(zlist)))
    cov_var = [[xvar for xvar in cvar.flatten()[half]]]
    dtime = t_var[-1] - t_var[-2]
    tnew = [0]

    for s_traj in zoftime[1:]:
        ind = 0
        for spi in sp_comp:
            conc[spi] = s_traj[ind]
            ind = ind + 1
        aa_jac = lna_ss_jacobian(
            lna_model_ss, s_traj, sp_comp, stch_var, ks_dict, r_dict, p_dict)
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)
        bb_diff = np.matmul(np.matmul(stch_var, np.diag(prop_flux.flatten())), stch_var.T)

        a_jac = np.nan_to_num(aa_jac)
        b_diff = np.nan_to_num(bb_diff)

        fofx = lna_cov_model(a_jac, b_diff, cvar)
        cvar = cvar + fofx * dtime
        cov_var.append([xvar for xvar in cvar.flatten()[half]])

    return [np.array(cov_var), si_new, t_var]


def lna_non_steady_state(conc, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var,
                         molar=True, rfile="", del_coef=10):
    get_globals(rfile)
    zlist = [conc[a] for a in sp_comp]
    # zoftime = odeint(
    #     lna_ode_model,zlist,t_var, args=(sp_comp,ks_dict,r_dict,p_dict,stch_var,molar))

    half = []
    si_new = []
    sps = list(sp_comp.keys())  # [svar for svar in sp_comp]
    len_sps = len(sps)
    for i in range(len_sps):
        svar = sps[i]
        for j in range(i, len_sps):
            pvar = sps[j]
            si_new.append("cov(" + svar + "," + pvar + ")")
            half.append(len_sps * i + j)

    cvar = np.zeros((len(zlist), len(zlist)))
    cov_var = [[xvar for xvar in cvar.flatten()[half]]]
    dtime = t_var[-1] - t_var[-2]
    tnew = [0]

    # for s_traj in zoftime[1:]:
    s_traj = np.array(zlist)
    while tnew[-1] < t_var[-1]:
        ind = 0
        for spi in sp_comp:
            conc[spi] = s_traj[ind]
            ind = ind + 1
        aa_jac = lna_ss_jacobian(
            lna_model_ss, s_traj, sp_comp, stch_var, ks_dict, r_dict, p_dict)
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)
        bb_diff = np.matmul(np.matmul(stch_var, np.diag(prop_flux.flatten())), stch_var.T)

        a_jac = np.nan_to_num(aa_jac)
        b_diff = np.nan_to_num(bb_diff)

        fx_covr = lna_cov_model(a_jac, b_diff, cvar)
        fx_conc = lna_ode_model(
            s_traj, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var, molar)

        h1_var = np.abs(fx_covr) + 1.0e-30
        dtime = max(min(del_coef * np.min(1 / h1_var), dtime), 1.0e-4)
        tnew.append(tnew[-1] + dtime)
        s_traj = s_traj + fx_conc * dtime
        cvar = cvar + fx_covr * dtime
        cov_var.append([xvar for xvar in cvar.flatten()[half]])

    return [np.array(cov_var), si_new, tnew]


def lna_non_steady_state2(conc, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var,
                          molar=True, rfile="", del_coef=10):
    zlist = [conc[a] for a in sp_comp]
    cov_var, slabels, tnew = lna_non_steady_state(
        conc, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var, molar, rfile, del_coef)
    zoftime = odeint(
        lna_ode_model, zlist, tnew, args=(sp_comp, ks_dict, r_dict, p_dict, stch_var, molar))

    si_new = []
    sps = [svar for svar in sp_comp]
    len_sps = len(sps)
    ff_div = []

    for xvar in slabels:
        si_new.append(xvar.replace("cov", "FF"))

    for s_traj in zoftime:
        row = []
        for i in range(len_sps):
            # svar = sps[i]
            for j in range(i, len_sps):
                # pvar = sps[j]
                s_ij = s_traj[i] * s_traj[j]
                row.append(np.sqrt(s_ij if s_ij != 0 else 1))
        ff_div.append(row)

    return [cov_var / np.array(ff_div), si_new, tnew]
