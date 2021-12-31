"""

                  This module is the tau_leaping2 module

The purpose of this module is to propagate stochastic trajectories using
the tau-leaping algorithm.

The functions in this module are the following;

1. tau_leaping2
2. step_3to5
3. ssa_support


"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.propagation.propensity import propensity_vec
from BioSANS2020.propagation.recalculate_globals \
    import get_globals, apply_rules, reserve_events_words


def tau_leaping2(tvar, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                 rand_seed, del_coef=1, implicit=False, rfile=""):
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
    dto = tvar[-1] - tvar[-2]
    tchlen = len(globals2.TCHECK)

    all_sp = list(sp_comp.keys())  # [z for z in sp_comp]
    spc = [z for z in sp_comp if z not in reserve_events_words]
    spc2 = [z for z in sp_comp if z in reserve_events_words]
    concz = {xvar: conc[xvar] for xvar in conc}
    yconc = {xvar: conc[xvar] for xvar in conc}
    apply_rules(concz, yconc)
    update_sp = [all_sp.index(z) for z in spc]

    vcri = np.where(stch_var < 0)
    vncr = np.array(list(range(len(stch_var[0]))))
    vncr = vncr[~np.isin(vncr, vcri[1])]

    zconc = [[concz[z] for z in spc]]
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
            prop_flux = propensity_vec(ks_dict, concz, r_dict, p_dict)
            if len(vcri[1]): # > 0:
                lcri = set()
                lncr = set()
                for xvar in range(len(vncr)):
                    lncr.add(xvar)

                for xvar in range(len(vcri[0])):
                    i, j = [vcri[0][xvar], vcri[1][xvar]]
                    if abs(zconc[-1][i] / stch_var[i, j]
                           ) < 10 and prop_flux[j] > 0:
                        lcri.add(j)
                    else:
                        lncr.add(j)
            else:
                lcri = set()

            if len(lcri) == len(stch_var[0]):
                alp = np.sum(prop_flux)
                d_t = (1 / alp) * (np.log(1 / np.random.uniform()))
            elif not lcri:  # len(lcri) == 0:
                alp = np.sum(prop_flux)
                uuj = np.matmul(stch_var, prop_flux)
                sig = np.matmul(stch_var2, prop_flux)
                exigi = np.array([zconc[-1][j] * (1 / g_i[j](zconc[-1][j]))
                                  for j in range(len(spc))]) * 0.03 * del_coef
                exigi1 = np.maximum(exigi, np.full(len(spc), 1))
                d_t = min(np.min(exigi1 / np.abs(uuj)),
                          np.min(exigi1 * exigi1 / np.abs(sig)))
            else:
                lcri = np.array(list(lcri))
                lncr = np.array(list(lncr))
                alpc = np.sum(prop_flux[lcri])
                dtc = (1 / alpc) * (np.log(1 / np.random.uniform()))

                alp = np.sum(prop_flux[lncr])
                uuj = np.matmul(stch_var[:, lncr], prop_flux[lncr])
                sig = np.matmul(stch_var2[:, lncr], prop_flux[lncr])
                exigi = np.array([zconc[-1][j] * (1 / g_i[j](zconc[-1][j]))
                                  for j in range(len(spc))]) * 0.03 * del_coef
                exigi1 = np.maximum(exigi, np.full(len(spc), 1))
                d_t = min(np.min(exigi1 / np.abs(uuj)),
                          np.min(exigi1 * exigi1 / np.abs(sig)), dtc)

            kvar = np.random.poisson(prop_flux * d_t)
            all_pos = True
            c_c = {}
            b_b = np.sum(kvar * stch_var[:, update_sp], 0)
            for xvar, _ in enumerate(spc):
                holder = zconc[-1][xvar] + b_b[xvar]
                if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                    c_c[spc[xvar]] = holder
                else:
                    all_pos = False
                    break
            if all_pos:
                for xvar, _ in enumerate(spc):
                    concz[spc[xvar]] = c_c[spc[xvar]]
                for xvar, _ in enumerate(spc2):
                    concz[spc2[xvar]] = concz[spc2[xvar]] + d_t
                apply_rules(concz, yconc)
                zconc.append([concz[xvar] for xvar in spc])
                t_c = t_c + d_t
                tnew.append(t_c)
    else:
        tindex = 1
        index = 0
        cvar = [concz[z] for z in spc]
        while tvar[tindex] < tmax:
            prop_flux = propensity_vec(ks_dict, concz, r_dict, p_dict)
            # step 1
            if len(vcri[1]):  # > 0:
                lcri = set()
                lncr = set()
                for xvar in range(len(vncr)):
                    lncr.add(xvar)

                for xvar in range(len(vcri[0])):
                    i, j = [vcri[0][xvar], vcri[1][xvar]]
                    if abs(zconc[-1][i] / stch_var[i, j]
                           ) < 10 and prop_flux[j] > 0:
                        lcri.add(j)
                    else:
                        lncr.add(j)
            else:
                lcri = set()

            # step 2
            epsilon = 0.03
            lcri = np.array(list(lcri))
            lncr = np.array(list(lncr))
            if not lncr:  # len(lncr) == 0:  # No non critical reactions
                dt1 = 1.0e+10
            else:
                alp = np.sum(prop_flux[lncr])
                uuj = np.matmul(stch_var[:, lncr], prop_flux[lncr])
                sig = np.matmul(stch_var2[:, lncr], prop_flux[lncr])
                exigi = np.array([cvar[j] * (1 / g_i[j](cvar[j]))
                                  for j in range(len(spc))]) \
                    * epsilon * del_coef
                exigi1 = np.maximum(exigi, np.full(len(spc), 1))
                dt1 = min(np.min(exigi1 / np.abs(uuj)),
                          np.min(exigi1 * exigi1 / np.abs(sig)))

            all_pos = False
            while not all_pos:
                kmul, d_t, do_ssa, _ = step_3to5(
                    prop_flux, lcri, dt1)  # step_3to5
                d_t = min(dto, d_t)
                if do_ssa:
                    t_c, tindex, index, break_now, all_pos, cvar = ssa_support(
                        tvar, ks_dict, r_dict, p_dict, stch_var, rfile,
                        tindex, index, t_c, zconc, spc, spc2, concz,
                        yconc, update_sp)
                    if break_now:
                        break
                else:
                    g_update = False
                    if index != tchlen:
                        if t_c + d_t > globals2.TCHECK[index]:
                            d_t = globals2.TCHECK[index] - t_c
                            index = index + 1
                            g_update = True

                    if t_c + d_t > tvar[tindex]:
                        d_t = max(0, tvar[tindex] - t_c)
                        kvar = np.round(
                            np.random.poisson(
                                prop_flux * d_t)) * kmul
                        all_pos = True
                        c_c = {}
                        b_b = np.sum(kvar * stch_var[:, update_sp], 0)
                        for xvar, _ in enumerate(spc):
                            holder = cvar[xvar] + b_b[xvar]
                            if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                                c_c[spc[xvar]] = holder
                            else:
                                all_pos = False
                                dt1 = dt1 / 2
                                break
                        if all_pos:
                            for xvar, _ in enumerate(spc):
                                concz[spc[xvar]] = c_c[spc[xvar]]
                            for xvar, _ in enumerate(spc2):
                                concz[spc2[xvar]] = concz[spc2[xvar]] + d_t
                            apply_rules(concz, yconc)
                            zconc.append([concz[xvar] for xvar in spc])
                            cvar = zconc[-1]
                            t_c = t_c + d_t
                            tindex = tindex + 1
                        else:
                            if g_update:
                                # pass
                                index = index - 1
                    else:
                        d_t = max(0, d_t)
                        kvar = np.round(
                            np.random.poisson(
                                prop_flux * d_t)) * kmul
                        all_pos = True
                        c_c = {}
                        b_b = np.sum(kvar * stch_var[:, update_sp], 0)
                        for xvar, _ in enumerate(spc):
                            holder = cvar[xvar] + b_b[xvar]
                            if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                                c_c[spc[xvar]] = holder
                            else:
                                all_pos = False
                                dt1 = dt1 / 2
                                break
                        if all_pos:
                            for xvar, _ in enumerate(spc):
                                concz[spc[xvar]] = c_c[spc[xvar]]
                            for xvar, _ in enumerate(spc2):
                                concz[spc2[xvar]] = concz[spc2[xvar]] + d_t
                            apply_rules(concz, yconc)
                            cvar = [concz[xvar] for xvar in spc]
                            t_c = t_c + d_t
                        else:
                            if g_update:
                                index = index - 1
                                # pass
            if tindex >= len(tvar):
                break
        if len(tvar) != len(zconc):
            while len(tvar) != len(zconc):
                zconc.append(zconc[-1])
        tnew = tvar
    return (tnew, np.array(zconc))


def step_3to5(prop_flux, lcri, dt1):
    """Additional steps in tau-leaping2"""
    # step 3
    kmul = np.ones((len(prop_flux), 1))
    r_1 = np.random.uniform()
    while r_1 == 0:
        r_1 = np.random.uniform()

    alpo = np.sum(prop_flux)
    do_ssa = False
    if dt1 < 10 * (1 / alpo):
        d_t = (1 / alpo) * (np.log(1 / r_1))
        do_ssa = True
    else:
        # step 4
        if lcri:  # > 0:
            alpc = np.sum(prop_flux[lcri])
            dt2 = (1 / alpc) * (np.log(1 / r_1))
        else:
            dt2 = 1.0e+5
            #dt2 = dt1

        # step 5
        if dt1 < dt2:
            d_t = dt1
            if lcri:  # > 0:
                kmul[lcri] = 0
        else:
            d_t = dt2
            if lcri:  # > 0:
                kmul[lcri] = 0
                pvar = np.cumsum([d / alpc for d in prop_flux[lcri]])
                r_2 = np.random.uniform()
                while r_2 == 0:
                    r_2 = np.random.uniform()
                for i, _ in enumerate(pvar):
                    if r_2 <= pvar[i]:
                        j_c = i
                        kmul[lcri[j_c]] = 1
                        break
    if d_t == 1.0e+5:
        do_ssa = True
    return [kmul, d_t, do_ssa, alpo]


def ssa_support(tvar, ks_dict, r_dict, p_dict, stch_var,
                rfile="", tindex=0, index=0, t_c=0, z_conc=[],
                spc=None, spc2=None, concz=None, yconc=None, update_sp=None):
    """SSA steps of tau-leaping2"""
    get_globals(rfile)
    stch_var = stch_var.T
    zconc = [[concz[z] for z in spc]]
    tmax = tvar[-1]

    break_now = False
    tchlen = len(globals2.TCHECK)
    all_pos = True
    for _ in range(100):
        if t_c < tmax:
            prop_flux = propensity_vec(ks_dict, concz, r_dict, p_dict)
            alp = np.sum(prop_flux)
            r_1 = np.random.uniform()
            while r_1 == 0:
                r_1 = np.random.uniform()
            pvar = np.cumsum([d / alp for d in prop_flux])
            d_t = (1 / alp) * (np.log(1 / r_1))

            #g_update = False
            if index != tchlen:
                if t_c + d_t >= globals2.TCHECK[index]:
                    d_t = globals2.TCHECK[index] - t_c
                    index = index + 1
                    #g_update = True

            if np.isnan(d_t) or np.isinf(d_t):
                z_conc.append(zconc[-1])
                while tvar[tindex] != tmax:
                    z_conc.append(zconc[-1])
                    tindex = tindex + 1
                t_c = tvar[-1]
                break_now = True
                break

            r_2 = np.random.uniform()
            for i, _ in enumerate(pvar):
                if r_2 <= pvar[i]:
                    all_pos = True
                    for xvar, _ in enumerate(spc):
                        holder = zconc[-1][xvar] + stch_var[i][update_sp[xvar]]
                        if holder >= 0 or spc[xvar] in globals2.MODIFIED:
                            concz[spc[xvar]] = holder
                        else:
                            all_pos = False
                            break
                    if all_pos:
                        for xvar, _ in enumerate(spc2):
                            concz[spc2[xvar]] = concz[spc2[xvar]] + d_t
                        apply_rules(concz, yconc)
                        zconc.append([concz[xvar] for xvar in spc])
                        t_c = t_c + d_t
                        try:
                            if t_c == tvar[tindex]:
                                z_conc.append(zconc[-1])
                                tindex = tindex + 1
                            else:
                                while t_c > tvar[tindex]:
                                    z_conc.append(zconc[-2])
                                    tindex = tindex + 1
                        except:
                            pass
                    else:
                        for xvar, _ in enumerate(spc):
                            concz[spc[xvar]] = zconc[-1][xvar]
                        # if g_update:
                            #index = index - 1
                    break
        else:
            break_now = True
            break
    return [t_c, tindex, index, break_now, all_pos, zconc[-1]]
