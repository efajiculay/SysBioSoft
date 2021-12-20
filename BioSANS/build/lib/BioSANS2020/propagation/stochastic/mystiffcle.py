"""

                  This module is the mystiffcle module

This can propagate non-stiff to moderately stiff stochastic simulation
using the chemical langevine equation. Here two versions are provided

1) Tau-adaptive CLE
2) Fix-inreval CLE

The following are the list of function for this module.

1. cle_model
2. cle_calculate
3. cle2_calculate


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec
from BioSANS2020.propagation.recalculate_globals \
    import get_globals, apply_rules, reserve_events_words
# from BioSANS2020.myglobal import mglobals as globals2


def cle_model(sp_comp, ks_dict, conc, r_dict, p_dict,
              stch_var, dtime, del_coef, reg=False):
    """This functions prepare the CLE model for integration

    Args:
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
        conc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
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
        stch_var (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        dtime (float): step-size
        del_coef (float): step-size factor or modifier
        reg (bool, optional): If True, the model is for fix-interval CLE
            . Defaults to False.

    Returns:
        np.ndarray: f_d = fofx * dtime + gofx * sqdt
    """
    prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict)
    sqrt_prop = np.sqrt(prop_flux)
    nlen = np.random.randn(len(prop_flux)).reshape(len(prop_flux), 1)
    fofx = np.matmul(stch_var, prop_flux)
    gofx = np.matmul(stch_var, sqrt_prop * nlen)
    if not reg:
        h_1 = np.abs(prop_flux) + 1.0e-30  # np.abs(fofx)+1.0e-30
        h_2 = np.abs(sqrt_prop) + 1.0e-30  # np.abs(gofx)+1.0e-30
        dtime = min(del_coef * min(np.min(1 / h_1), np.min(1 / h_2)), dtime)
    sqdt = np.sqrt(dtime)

    ind = 0
    for s_p in sp_comp:
        key = s_p.strip().split("_")[0]
        if key in reserve_events_words:
            gofx[ind] = 0
        ind = ind + 1

    f_d = fofx * dtime + gofx * sqdt
    return [f_d.reshape(len(sp_comp)), dtime]


def cle_calculate(tvar, sp_comp, ks_dict, sconc, r_dict, p_dict, stch_var,
                  del_coef=10, rand_seed=1, implicit=False, rfile=""):
    """This functions performs the tau-adaptive CLE integration
    Args:
        tvar (list): time stamp of simulation
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
        sconc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
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
        stch_var (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        del_coef (float): step-size factor or modifier
        rand_seed (float): random  seed value picked  at random for each
            trajectory. They have been sampled from the calling program.
        implicit (bool, optional): True means report in time intervals
            similar to the input time intervals even if actual step is
            more or less. Defaults to False.
        rfile (string, optional): name of topology file where some
            parameters or components are negative indicating  they  have
            to be estimated. Defaults to None.

    Returns:
        tuple: (time, trajectories)
    """
    get_globals(rfile)
    np.random.seed(int(rand_seed * 100))
    tnew = []
    dtime = tvar[-1] - tvar[-2]
    yconc = {x: sconc[x] for x in sconc}
    conc = {x: sconc[x] for x in sconc}
    apply_rules(conc, yconc)
    res_traj = [[conc[z] for z in sp_comp]]

    if not implicit:
        tnow = tvar[0]
        tnew.append(tnow)
        while tnow < tvar[-1]:
            mvar = np.nan_to_num(
                cle_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                          stch_var, dtime, del_coef))
            mmvar = mvar[0].reshape(1, len(mvar[0]))[0]
            tnow = tnow + mvar[1]
            ind = 0
            for s_p in sp_comp:
                conc[s_p] = res_traj[-1][ind] + mmvar[ind]
                ind = ind + 1
            apply_rules(conc, yconc)
            res_traj.append([conc[z] for z in sp_comp])
            for s_p in sp_comp:
                conc[s_p] = max(0, conc[s_p])
            tnew.append(tnow)
    else:
        tnow = tvar[0]
        tnew.append(tnow)
        tindex = 1
        cvar = [conc[z] for z in sp_comp]
        while tnew[-1] < tvar[-1]:
            mvar = np.nan_to_num(
                cle_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                          stch_var, dtime, del_coef))
            tnow = tnow + mvar[1]
            if tnow > tvar[tindex]:
                tnow = tnow - mvar[1]
                dt2 = tvar[tindex] - tnow
                mvar = np.nan_to_num(
                    cle_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                              stch_var, dt2, del_coef, True))
                mmvar = mvar[0].reshape(1, len(mvar[0]))[0]
                tnow = tnow + mvar[1]

                ind = 0
                for s_p in sp_comp:
                    conc[s_p] = cvar[ind] + mmvar[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                cvar = [conc[z] for z in sp_comp]
                for s_p in sp_comp:
                    conc[s_p] = max(0, conc[s_p])
                res_traj.append(cvar)

                tnew.append(tnow)
                tindex = tindex + 1
            else:
                mmvar = mvar[0].reshape(1, len(mvar[0]))[0]
                ind = 0
                for s_p in sp_comp:
                    conc[s_p] = cvar[ind] + mmvar[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                cvar = [conc[z] for z in sp_comp]
                for s_p in sp_comp:
                    conc[s_p] = max(0, conc[s_p])
    return (tnew, np.array(res_traj))


def cle2_calculate(tvar, sp_comp, ks_dict, sconc, r_dict, p_dict,
                   stch_var, del_coef=1, rand_seed=1, rfile=""):
    """This functions performs the fix-interval CLE integration
    Args:
        tvar (list): time stamp of simulation
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
        sconc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
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
        stch_var (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        del_coef (float): step-size factor or modifier
        rand_seed (float): random  seed value picked  at random for each
            trajectory. They have been sampled from the calling program.
        rfile (string, optional): name of topology file where some
            parameters or components are negative indicating  they  have
            to be estimated. Defaults to None.

    Returns:
        tuple: (time, trajectories)
    """
    get_globals(rfile)
    np.random.seed(int(rand_seed * 100))
    div = max(1, int(1 / del_coef))
    dtime = (tvar[-1] - tvar[-2]) / div
    yconc = {x: sconc[x] for x in sconc}
    conc = {x: sconc[x] for x in sconc}
    apply_rules(conc, yconc)
    res_traj = [[conc[z] for z in sp_comp]]
    tnow = tvar[0]
    tnew = [tnow]
    d_v = 0
    cvar = res_traj[-1]
    while abs(tnow - tvar[-1]) > 1.0e-10:
        mvar = np.nan_to_num(
            cle_model(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                      dtime, 1, True))
        mmvar = mvar[0].reshape(1, len(mvar[0]))[0]
        tnow = tnow + mvar[1]
        ind = 0
        for s_p in sp_comp:
            conc[s_p] = cvar[ind] + mmvar[ind]
            ind = ind + 1
        apply_rules(conc, yconc)
        cvar = [conc[z] for z in sp_comp]
        if d_v == div - 1:
            res_traj.append(cvar)
            tnew.append(tnow)
            d_v = 0
        else:
            d_v = d_v + 1
        for s_p in sp_comp:
            conc[s_p] = max(0, conc[s_p])

    return (tnew, np.array(res_traj))
