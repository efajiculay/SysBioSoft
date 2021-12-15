"""

             This module is the lna_function_of_time module

     This module performs numerical propagation of linear noise approxi-
mation or LNA by exploiting the following relationship.

                   dC/dt = AC + CA.T + BB

where A is defined as d(S*f)/dx where S is the stoichiometric  matrix, f
are the propensities or fluxes, and x are the components or species. The
flux is a function of species x and rate constant k. BB is the diffusion
matrix equivalent to S*diag(f)*S.T where diag(f) is a square matrix with
zero non-diagonal elements and f[i] in each diagonal elements. C is the
covariance matrix and t is time.

The following are the list of function for this module.


"""

from scipy.integrate import odeint
import numpy as np
# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.propagation.deterministic.lna_approx \
    import lna_ss_jacobian, lna_model_ss
from BioSANS2020.propagation.recalculate_globals import get_globals
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar


def lna_ode_model(zlist, t_var, sp_comp, ks_dict, r_dict, p_dict,
                  stch_var, molar=False):
    """This  function returns the derivative of components with respect
    to time at a particular state of the system based on the inputs.

    Args:
        zlist (list): list of components or species amounts
        sp_comp (dict): dictionary of appearance or position of species
            or component in the reaction tag of BioSANS topology file.

            For example;

                #REACTIONS
                A => B, -kf1    # negative means to be estimated
                B => C, kf2

            The value of sp_comp is

                sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}

                A appears in first reaction with index 0
                B appears in first and second reaction with index 0, 1
                C appears in second reaction with index 1
        ks_dict (dict): dictionary of rate constant that appear in each
            reactions.

            For example;

                #REACTIONS
                A => B , 0.3        # first reaction
                B <=> C, 0.1, 0.2   # second reaction

            The value of ks_dict is

                ks_dict = {
                    0 : [0.3],      # first reaction
                    1 : [0.1, 0.2]  # second reaction
                }
        r_dict (dict): dictionary of reactant stoichiometry. For example

            r_dict = {
                0: {'A': 1},  # first reaction, coefficient of A is 1
                1: {'B': 1}   # second reaction, coefficient of B is 1
            }
        p_dict (dict): dictionary of product stoichiometry. For example

            p_dict = {
                0: {'B': 1},  # first reaction, coefficient of B is 1
                1: {'C': 1}   # second reaction, coefficient of C is 1
            }
        stch_var (np.ndarray): stoichiometric matrix
        molar (boolean) : if yes, uses macroscopic propensities else
            uses microscopic propensities
    Returns:
        np.ndarray: derivative of species with respect to time.d(S*f)/dt
    """
    spc = list(sp_comp.keys())
    conc = {spc[xvar]: zlist[xvar] for xvar in range(len(spc))}
    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        prop_flux = propensity_vec_molar(
            ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(stch_var, prop_flux).reshape(len(zlist))
    for xvar in globals2.CON_BOUNDARY:
        ind = spc.index(xvar)
        dxdt[ind] = 0
    return dxdt


def lna_cov_model(a_jac, b_diff, cov):
    """Returns the evaluated dC/dt at a particular instant based on the
    current state of a_jac, b_diff, cov

    Args:
        a_jac (np.ndarray):  d(S*f)/dx or jacobian of ODE with respect
            to components or species
        b_diff (np.ndarray): diffusuin matrix
        cov (np.ndarray): covariance matrix

    Returns:
        [type]: [description]
    """
    return np.matmul(a_jac, cov) + np.matmul(cov, np.transpose(a_jac)) + b_diff


def lna_non_steady_state(conc, t_var, sp_comp, ks_dict, r_dict, p_dict,
                         stch_var, molar=True, rfile="", del_coef=10):
    """This function returns covariance trajectoy based on LNA

    Args:
        conc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
        t_var (list): time stamp of trajectories
        sp_comp (dict): dictionary of appearance or position of species
            or component in the reaction tag of BioSANS topology file.

            For example;

                #REACTIONS
                A => B, -kf1    # negative means to be estimated
                B => C, kf2

            The value of sp_comp is

                sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}

                A appears in first reaction with index 0
                B appears in first and second reaction with index 0, 1
                C appears in second reaction with index 1
        ks_dict (dict): dictionary of rate constant that appear in each
            reactions.

            For example;

                #REACTIONS
                A => B , 0.3        # first reaction
                B <=> C, 0.1, 0.2   # second reaction

            The value of ks_dict is

                ks_dict = {
                    0 : [0.3],      # first reaction
                    1 : [0.1, 0.2]  # second reaction
                }
        r_dict (dict): dictionary of reactant stoichiometry. For example

            r_dict = {
                0: {'A': 1},  # first reaction, coefficient of A is 1
                1: {'B': 1}   # second reaction, coefficient of B is 1
            }
        p_dict (dict): dictionary of product stoichiometry. For example

            p_dict = {
                0: {'B': 1},  # first reaction, coefficient of B is 1
                1: {'C': 1}   # second reaction, coefficient of C is 1
            }
        stch_var (np.ndarray): stoichiometric matrix
        molar : True if conc. is in molar otherwise False
        rfile (str): file name of BioSANS topology file.
        del_coef (float, optional): factor for modifying time steps used
            in the integration/propagation of ODE. Defaults to 10.

    Returns:
        list: [covariance, labels, time]
    """
    get_globals(rfile)
    zlist = [conc[a] for a in sp_comp]
    # zoftime = odeint(
    #     lna_ode_model,zlist,t_var,
    #     args=(sp_comp,ks_dict,r_dict,p_dict,stch_var,molar))

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
        bb_diff = np.matmul(np.matmul(stch_var, np.diag(prop_flux.flatten())),
                            stch_var.T)

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


def lna_non_steady_state2(conc, t_var, sp_comp, ks_dict, r_dict, p_dict,
                          stch_var, molar=True, rfile="", del_coef=10):
    """This function returns Fan-factor trajectoy based on LNA

    Args:
        conc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
        t_var (list): time stamp of trajectories
        sp_comp (dict): dictionary of appearance or position of species
            or component in the reaction tag of BioSANS topology file.

            For example;

                #REACTIONS
                A => B, -kf1    # negative means to be estimated
                B => C, kf2

            The value of sp_comp is

                sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}

                A appears in first reaction with index 0
                B appears in first and second reaction with index 0, 1
                C appears in second reaction with index 1
        ks_dict (dict): dictionary of rate constant that appear in each
            reactions.

            For example;

                #REACTIONS
                A => B , 0.3        # first reaction
                B <=> C, 0.1, 0.2   # second reaction

            The value of ks_dict is

                ks_dict = {
                    0 : [0.3],      # first reaction
                    1 : [0.1, 0.2]  # second reaction
                }
        r_dict (dict): dictionary of reactant stoichiometry. For example

            r_dict = {
                0: {'A': 1},  # first reaction, coefficient of A is 1
                1: {'B': 1}   # second reaction, coefficient of B is 1
            }
        p_dict (dict): dictionary of product stoichiometry. For example

            p_dict = {
                0: {'B': 1},  # first reaction, coefficient of B is 1
                1: {'C': 1}   # second reaction, coefficient of C is 1
            }
        stch_var (np.ndarray): stoichiometric matrix
        molar : True if conc. is in molar otherwise False
        rfile (str): file name of BioSANS topology file.
        del_coef (float, optional): factor for modifying time steps used
            in the integration/propagation of ODE. Defaults to 10.

    Returns:
        list: [fanp_factor, labels, time]
    """
    zlist = [conc[a] for a in sp_comp]
    cov_var, slabels, tnew = lna_non_steady_state(
        conc, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var, molar,
        rfile, del_coef)
    zoftime = odeint(
        lna_ode_model, zlist, tnew, args=(sp_comp, ks_dict, r_dict, p_dict,
                                          stch_var, molar))

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


# def lna_non_steady_state_old(
#         conc, t_var, sp_comp, ks_dict, r_dict, p_dict, stch_var,
#         molar=True, rfile="", del_coef=10):
#     get_globals(rfile)
#     zlist = [conc[a] for a in sp_comp]
#     zoftime = odeint(
#         lna_ode_model, zlist, t_var,
#         args=(sp_comp, ks_dict, r_dict, p_dict, stch_var, molar))
#
#     half = []
#     si_new = []
#     sps = [svar for svar in sp_comp]
#     len_sps = len(sps)
#     for i in range(len_sps):
#         svar = sps[i]
#         for j in range(i, len_sps):
#             pvar = sps[j]
#             si_new.append("cov(" + svar + "," + pvar + ")")
#             half.append(len_sps * i + j)
#
#     cvar = np.zeros((len(zlist), len(zlist)))
#     cov_var = [[xvar for xvar in cvar.flatten()[half]]]
#     dtime = t_var[-1] - t_var[-2]
#     tnew = [0]
#
#     for s_traj in zoftime[1:]:
#         ind = 0
#         for spi in sp_comp:
#             conc[spi] = s_traj[ind]
#             ind = ind + 1
#         aa_jac = lna_ss_jacobian(
#             lna_model_ss, s_traj, sp_comp, stch_var, ks_dict, r_dict, p_dict)
#         prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)
#         bb_diff = np.matmul(np.matmul(stch_var, np.diag(prop_flux.flatten())),
#                             stch_var.T)
#
#         a_jac = np.nan_to_num(aa_jac)
#         b_diff = np.nan_to_num(bb_diff)
#
#         fofx = lna_cov_model(a_jac, b_diff, cvar)
#         cvar = cvar + fofx * dtime
#         cov_var.append([xvar for xvar in cvar.flatten()[half]])
#
#     return [np.array(cov_var), si_new, t_var]
