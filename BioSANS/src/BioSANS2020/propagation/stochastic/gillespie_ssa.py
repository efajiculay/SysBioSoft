"""

                 This module is the gillespie_ssa module

This module propagates the stochastic simulation algorithm also known as
Gillespie.

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec
from BioSANS2020.propagation.recalculate_globals \
    import get_globals, apply_rules, reserve_events_words
from BioSANS2020.myglobal import mglobals as globals2


def gillespie_ssa(tvar, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                  rand_seed, implicit=False, rfile=""):
    """[summary]

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
    tmax = tvar[-1]
    np.random.seed(int(rand_seed * 100))
    tnew = []
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
    if not implicit:
        while t_c < tmax:
            prop_flux = propensity_vec(ks_dict, concz, r_dict, p_dict)
            alp = np.sum(prop_flux)
            r_1 = np.random.uniform()
            while r_1 == 0:
                r_1 = np.random.uniform()
            p_r = np.cumsum([d / alp for d in prop_flux])
            d_t = (1 / alp) * (np.log(1 / r_1))
            if np.isnan(d_t) or np.isinf(d_t):
                break
            r_2 = np.random.uniform()
            for i, _ in enumerate(p_r):
                if r_2 <= p_r[i]:
                    all_pos = True
                    for xvar, _ in enumerate(spc):
                        holder = zlist[-1][xvar] + stch_var[i][update_sp[xvar]]
                        if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                            concz[spc[xvar]] = holder
                        else:
                            all_pos = True
                            break
                    if all_pos:
                        for xvar, _ in enumerate(spc2):
                            concz[spc2[xvar]] = concz[spc2[xvar]] + d_t
                        apply_rules(concz, yconc)
                        zlist.append([concz[xvar] for xvar in spc])
                        t_c = t_c + d_t
                        tnew.append(t_c)
                    else:
                        for xvar, _ in enumerate(spc):
                            concz[spc[xvar]] = zlist[-1][xvar]
                    break
    else:
        z_c = []
        tindex = 0
        index = 0
        tchlen = len(globals2.TCHECK)
        while t_c < tmax:
            prop_flux = propensity_vec(ks_dict, concz, r_dict, p_dict)
            alp = np.sum(prop_flux)
            r_1 = np.random.uniform()
            while r_1 == 0:
                r_1 = np.random.uniform()
            p_r = np.cumsum([d / alp for d in prop_flux])
            d_t = (1 / alp) * (np.log(1 / r_1))

            if index != tchlen:
                if t_c + d_t >= globals2.TCHECK[index]:
                    d_t = globals2.TCHECK[index] - t_c
                    index = index + 1

            if np.isnan(d_t) or np.isinf(d_t):
                z_c.append(zlist[-1])
                while tvar[tindex] != tmax:
                    z_c.append(zlist[-1])
                    tindex = tindex + 1
                break
            r_2 = np.random.uniform()
            for i, _ in enumerate(p_r):
                if r_2 <= p_r[i]:
                    all_pos = True
                    for xvar, _ in enumerate(spc):
                        holder = zlist[-1][xvar] + stch_var[i][update_sp[xvar]]
                        if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                            concz[spc[xvar]] = holder
                        else:
                            all_pos = False
                            break
                    if all_pos:
                        for xvar, _ in enumerate(spc2):
                            concz[spc2[xvar]] = concz[spc2[xvar]] + d_t

                        t_c = t_c + d_t
                        if "t" in spc2:
                            concz["t"] = t_c
                        elif "time" in spc2:
                            concz["time"] = t_c
                        else:
                            pass

                        apply_rules(concz, yconc)
                        zlist.append([concz[xvar] for xvar in spc])
                        # t_c = t_c + d_t
                        try:
                            if t_c == tvar[tindex]:
                                z_c.append(zlist[-1])
                                tindex = tindex + 1
                            else:
                                while t_c > tvar[tindex]:
                                    z_c.append(zlist[-2])
                                    tindex = tindex + 1
                        except:
                            pass
                        tnew.append(t_c)
                    else:
                        for xvar, _ in enumerate(spc):
                            concz[spc[xvar]] = zlist[-1][xvar]
                    break
        tnew = tvar
        zlist = z_c
    return (tnew, np.array(zlist))
