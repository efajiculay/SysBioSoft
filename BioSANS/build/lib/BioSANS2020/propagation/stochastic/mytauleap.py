"""

                  This module is the mytauleap module

The purpose of this module is to propagate stochastic trajectories using
the tau-leaping algorithm.

The functions in this module are the following;

1. sim_tauLeap
2.


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.propagation.recalculate_globals \
    import get_globals, apply_rules, reserve_events_words


def sim_tauLeap(tvar, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                rand_seed, del_coef, implicit=False, rfile=""):
    """This functions performs the tau-leaping integration

    Args:
        tvar ([type]): [description]
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
        rand_seed ([type]): [description]
        rand_seed (float): random  seed value picked  at random for each
            trajectory. They have been sampled from the calling program.
        implicit (bool, optional): True means report in time intervals
            similar to the input time intervals even if actual step is
            more or less. Defaults to False.
        rfile (string, optional): name of topology file where some
            parameters or components are negative indicating  they  have
            to be estimated. Defaults to None.

    Returns:
        tuple : (time, trajectories)
    """
    get_globals(rfile)
    tmax = tvar[-1]
    np.random.seed(int(rand_seed * 100))
    tnew = []
    stch_var2 = stch_var**2
    stch_var = stch_var.T

    all_sp = list(sp_comp.keys())  # [z for z in sp_comp]
    spc = [z for z in sp_comp if z not in reserve_events_words]
    spc2 = [z for z in sp_comp if z in reserve_events_words]
    concz = {xvar: conc[xvar] for xvar in conc}
    yconc = {xvar: conc[xvar] for xvar in conc}
    apply_rules(concz, yconc)
    update_sp = [all_sp.index(z) for z in spc]
    zlist = [[concz[z] for z in spc]]
    t_c = 0
    tnew.append(t_c)

    g_i = []
    for i, _ in enumerate(spc):
        keys = sp_comp[spc[i]]
        maxrlen = 0
        rlen = 0
        for key in keys:
            if spc[i] in r_dict[key]:
                maxrlen = max(maxrlen, len(r_dict[key]))
            elif spc[i] in p_dict[key]:
                maxrlen = max(maxrlen, len(p_dict[key]))
        if maxrlen == 1:
            for key in keys:
                if spc[i] in r_dict[key]:
                    rlen = max(rlen, r_dict[key][spc[i]])
                elif spc[i] in p_dict[key]:
                    rlen = max(rlen, p_dict[key][spc[i]])
            if rlen == 1:
                g_i.append(lambda xvar: 1)
            else:
                g_i.append(lambda xvar: 2)
        elif maxrlen == 2:
            try:
                g_i.append(lambda xvar: 2 + 1 / (xvar - 1))
            except:
                g_i.append(lambda xvar: 2)
        else:
            g_i.append(lambda xvar: 1)

    if not implicit:
        while t_c < tmax:
            for xvar, _ in enumerate(spc):
                concz[spc[xvar]] = zlist[-1][xvar]
            prop_flux = propensity_vec(ks_dict, concz, r_dict, p_dict)
            alp = np.sum(prop_flux)
            r_1 = np.random.uniform()
            while r_1 == 0:
                r_1 = np.random.uniform()
            dt1 = (1 / alp) * (np.log(1 / r_1))

            uuj = np.matmul(stch_var, prop_flux)
            sig = np.matmul(stch_var2, prop_flux)
            exigi = np.array([zlist[-1][j] * (1 / g_i[j](zlist[-1][j]))
                              for j in range(len(spc))]) * 0.03 * del_coef
            exigi1 = np.maximum(exigi, np.full(len(spc), 1))
            dt2 = min(np.min(exigi1 / np.abs(uuj)),
                      np.min(exigi1 * exigi1 / np.abs(sig)))

            dtime = max(dt1, dt2)
            kvar = np.round(np.random.poisson(prop_flux * dtime))
            all_pos = True
            c_c = {}
            b_b = np.sum(kvar * stch_var[:, update_sp], 0)
            for xvar, _ in enumerate(spc):
                holder = zlist[-1][xvar] + b_b[xvar]
                if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                    c_c[spc[xvar]] = holder
                else:
                    all_pos = False
                    break
            if all_pos:
                for xvar, _ in enumerate(spc):
                    concz[spc[xvar]] = c_c[spc[xvar]]
                for xvar, _ in enumerate(spc2):
                    concz[spc2[xvar]] = concz[spc2[xvar]] + dtime
                apply_rules(concz, yconc)
                zlist.append([concz[xvar] for xvar in spc])
                t_c = t_c + dtime
                tnew.append(t_c)

    else:
        tindex = 1
        z_c = [zlist[-1]]
        while t_c < tmax:
            prop_flux = propensity_vec(ks_dict, concz, r_dict, p_dict)
            prop_flux[prop_flux < 0] = 0
            alp = np.sum(prop_flux)
            r_1 = np.random.uniform()
            while r_1 == 0:
                r_1 = np.random.uniform()
            dt1 = (1 / alp) * (np.log(1 / r_1))

            uuj = np.matmul(stch_var, prop_flux)
            sig = np.matmul(stch_var2, prop_flux)
            exigi = np.array([zlist[-1][j] * (1 / g_i[j](zlist[-1][j]))
                              for j in range(len(spc))]) * 0.03 * del_coef
            exigi1 = np.maximum(exigi, np.full(len(spc), 1))
            dt2 = min(np.min(exigi1 / np.abs(uuj)),
                      np.min(exigi1 * exigi1 / np.abs(sig)))

            dtime = max(dt1, dt2)
            if t_c + dtime > tvar[tindex]:
                dtime = max(0, tvar[tindex] - t_c)
                kvar = np.round(np.random.poisson(prop_flux * dtime))
                all_pos = True
                c_c = {}
                b_b = np.sum(kvar * stch_var[:, update_sp], 0)
                for xvar, _ in enumerate(spc):
                    holder = zlist[-1][xvar] + b_b[xvar]
                    if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                        c_c[spc[xvar]] = holder
                    else:
                        all_pos = False
                        break
                if all_pos:
                    for xvar, _ in enumerate(spc):
                        concz[spc[xvar]] = c_c[spc[xvar]]
                    for xvar, _ in enumerate(spc2):
                        concz[spc2[xvar]] = concz[spc2[xvar]] + dtime
                    apply_rules(concz, yconc)
                    zlist.append([concz[xvar] for xvar in spc])
                    z_c.append(zlist[-1])
                    t_c = t_c + dtime
                    tindex = tindex + 1
            else:
                dtime = max(0, dtime)
                kvar = np.round(np.random.poisson(prop_flux * dtime))
                all_pos = True
                c_c = {}
                b_b = np.sum(kvar * stch_var[:, update_sp], 0)
                for xvar, _ in enumerate(spc):
                    holder = zlist[-1][xvar] + b_b[xvar]
                    if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                        c_c[spc[xvar]] = holder
                    else:
                        all_pos = False
                        break
                if all_pos:
                    for xvar, _ in enumerate(spc):
                        concz[spc[xvar]] = c_c[spc[xvar]]
                    for xvar, _ in enumerate(spc2):
                        concz[spc2[xvar]] = concz[spc2[xvar]] + dtime
                    apply_rules(concz, yconc)
                    zlist.append([concz[xvar] for xvar in spc])
                    t_c = t_c + dtime

        if len(tvar) != len(z_c):
            while len(tvar) != len(z_c):
                z_c.append(zlist[-1])
        tnew = tvar
        zlist = z_c
    return (tnew, np.array(zlist))
