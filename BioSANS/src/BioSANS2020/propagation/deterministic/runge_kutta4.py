"""

                 This module is the runge_kutta4 module

This module serves in the integration of ODE trajectory using RK4 and
tau-adaptive RK4-algorithm.

The list of functions in this module are

1. rk4_model
2. runge_kutta_forth
3. rungek4_int
rkck
rkqs
rungek4a_int


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals, \
    apply_rules
from BioSANS2020.myglobal import mglobals as globals2


def rk4_model(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
              tvar, molar=False):
    """This function returns  the evaluated  value of the  derivative of
    each component with respwct to time at a particular instant based on
    the state of the system at that instant.

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
        ksn_dict (dict): dictionary of rate constant that appear in each
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
        stoch_var (numpy.ndarray): stoichiometric matrix. For example

            stoch_var = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        tvar (list): time stamp of trajectories i.e. [0, 0.1, 0.2, ...]
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.

    Returns:
        np.ndarray: value of dx/dt
    """
    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict)
    else:
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict)
    dxdt = np.matmul(stch_var, prop_flux).reshape(len(sp_comp))
    if globals2.CON_BOUNDARY:
        spc = [spi for spi in sp_comp]
        for xvar in globals2.CON_BOUNDARY:
            ind = spc.index(xvar)
            dxdt[ind] = 0
    return dxdt


def runge_kutta_forth(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                      tvar, delt, n_sp, molar=False):
    """this function prepare some stuffs for the integration

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
        ks_dict (dict): dictionary of rate constant that appears in each
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
        conc ([type]): [description]
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
        v_stoich (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        tvar (list): time stamp of trajectories i.e. [0, 0.1, 0.2, ...]
        delt (float): step size
        n_sp (dict): dictionary of keywords
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.

    Returns:
        list: [updated concn, updated time]
    """

    yconc = {xvar: conc[xvar] for xvar in n_sp}

    da_dt = rk4_model(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                      tvar, molar)
    new_a1 = da_dt * delt
    ind = 0
    for spi in sp_comp:
        yconc[spi] = conc[spi] + 0.5 * new_a1[ind]
        ind = ind + 1

    da_dt = rk4_model(sp_comp, ks_dict, yconc, r_dict, p_dict, stch_var,
                      tvar + 0.5 * delt, molar)
    new_a2 = da_dt * delt
    ind = 0
    for spi in sp_comp:
        yconc[spi] = conc[spi] + 0.5 * new_a2[ind]
        ind = ind + 1

    da_dt = rk4_model(sp_comp, ks_dict, yconc, r_dict, p_dict, stch_var,
                      tvar + 0.5 * delt, molar)
    new_a3 = da_dt * delt
    ind = 0
    for spi in sp_comp:
        yconc[spi] = conc[spi] + 0.5 * new_a3[ind]
        ind = ind + 1

    da_dt = rk4_model(sp_comp, ks_dict, yconc, r_dict, p_dict, stch_var,
                      tvar + delt, molar)
    new_a4 = da_dt * delt

    new = (new_a1 + 2.0 * new_a2 + 2.0 * new_a3 + new_a4) / 6.0

    new_a = yconc
    ind = 0
    for spi in sp_comp:
        new_a[spi] = conc[spi] + new[ind]
        ind = ind + 1

    return [new_a, tvar + delt]


def rungek4_int(conc, time, sp_comp, ks_dict, r_dict, p_dict, stch_var,
                molar=False, delx=1, rfile=""):
    """Peforms the RK4 integration

    Args:
        conc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
        tvar (list): time stamp of trajectories i.e. [0, 0.1, 0.2, ...]
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
        ks_dict (dict): dictionary of rate constant that appears in each
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
        stch_var (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.
        delx (float, optional): stepsize modifier. Defaults to 1.
        rfile (string, optional): name of topology file where some
            parameters or components are negative indicating  they  have
            to be estimated. Defaults to None.

    Returns:
        list: [time stamp, trajectory]
    """
    get_globals(rfile)
    tend = time[-1]
    div = max(1, int(1 / delx))
    delt = (time[-1] - time[-2]) / div
    yconc = {xvar: conc[xvar] for xvar in conc}
    slabels = [a for a in sp_comp]
    apply_rules(conc, yconc, [0], [conc[a] for a in sp_comp], slabels)
    zvar = [conc[a] for a in sp_comp]
    zlist = [zvar]
    tvar = 0.0
    t_c = [tvar]
    n_sp = [xvar for xvar in conc if xvar not in sp_comp]
    d_v = 0
    tc2 = [tvar]
    z_2 = [zvar]

    while abs(tvar - tend) > 1.0e-10:
        conc, tvar = runge_kutta_forth(
            sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, tvar,
            delt, n_sp, molar)
        # conc, yerr = rkck(delt,sp_comp,ks_dict,conc,r_dict,p_dict,
        # stch_var,tvar,n_sp,molar) #higher order runge-kutta
        # tvar = tvar + delt     # higher order
        # runge-kutta
        apply_rules(conc, yconc, tc2, z_2, slabels)
        z_2.append([conc[a] for a in sp_comp])
        tc2.append(tvar)
        if d_v == div - 1:
            zlist.append([conc[a] for a in sp_comp])
            t_c.append(tvar)
            d_v = 0
        else:
            d_v = d_v + 1
    zlist = np.array(zlist)
    return [t_c, zlist]


def rkck(hvar, sp_comp, ks_dict, conc, r_dict, p_dict,
         stch_var, tvar, n_sp, molar):
    """This function is a helper function for tau-adaptive RK4"""

    a_2 = 0.2
    a_3 = 0.3
    a_4 = 0.6
    a_5 = 1.0
    a_6 = 0.875

    b_21 = 0.2
    b_31 = 3.0 / 40.0
    b_32 = 9.0 / 40.0

    b_41 = 0.3
    b_42 = -0.9
    b_43 = 1.2

    b_51 = -11.9 / 54.0
    b_52 = 2.5
    b_53 = -70.0 / 27.0
    b_54 = 35.0 / 27.0

    b_61 = 1631.0 / 55296.0
    b_62 = 175.0 / 512.0
    b_63 = 575.0 / 13824.0
    b_64 = 44275.0 / 110592.0
    b_65 = 253.0 / 4096.0

    c_1 = 37.0 / 378.0
    c_3 = 250.0 / 621.0
    c_4 = 125.0 / 594.0
    c_6 = 512.0 / 1771.0

    dc_1 = c_1 - 2825.0 / 27648.0
    dc_3 = c_3 - 18575.0 / 48384.0
    dc_4 = c_4 - 13525.0 / 55296.0
    dc_5 = -277.0 / 14336.0
    dc_6 = c_6 - 0.25

    dydx = rk4_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                     stch_var, tvar, molar)

    ytemp = {xvar: conc[xvar] for xvar in n_sp}
    ind = 0
    for spi in sp_comp:
        ytemp[spi] = conc[spi] + b_21 * hvar * dydx[ind]
        ind = ind + 1
    ak2 = rk4_model(sp_comp, ks_dict, ytemp, r_dict, p_dict,
                    stch_var, tvar + a_2 * hvar, molar)

    ind = 0
    for spi in sp_comp:
        ytemp[spi] = conc[spi] + hvar * (b_31 * dydx[ind] + b_32 * ak2[ind])
        ind = ind + 1
    ak3 = rk4_model(sp_comp, ks_dict, ytemp, r_dict, p_dict,
                    stch_var, tvar + a_3 * hvar, molar)

    ind = 0
    for spi in sp_comp:
        ytemp[spi] = conc[spi] + hvar * \
            (b_41 * dydx[ind] + b_42 * ak2[ind] + b_43 * ak3[ind])
        ind = ind + 1
    ak4 = rk4_model(sp_comp, ks_dict, ytemp, r_dict, p_dict,
                    stch_var, tvar + a_4 * hvar, molar)

    ind = 0
    for spi in sp_comp:
        ytemp[spi] = conc[spi] + hvar * \
            (b_51 * dydx[ind] + b_52 * ak2[ind] +
             b_53 * ak3[ind] + b_54 * ak4[ind])
        ind = ind + 1
    ak5 = rk4_model(sp_comp, ks_dict, ytemp, r_dict, p_dict,
                    stch_var, tvar + a_5 * hvar, molar)

    ind = 0
    for spi in sp_comp:
        ytemp[spi] = conc[spi] + hvar * (b_61 * dydx[ind] + b_62
                                         * ak2[ind] + b_63 * ak3[ind]
                                         + b_64 * ak4[ind] + b_65 * ak5[ind])
        ind = ind + 1
    ak6 = rk4_model(sp_comp, ks_dict, ytemp, r_dict, p_dict,
                    stch_var, tvar + a_6 * hvar, molar)

    yout = ytemp
    yerr = []
    ind = 0
    for spi in sp_comp:
        yout[spi] = conc[spi] + hvar * \
            (c_1 * dydx[ind] + c_3 * ak3[ind] +
             c_4 * ak4[ind] + c_6 * ak6[ind])
        yerr.append(hvar * (dc_1 * dydx[ind] + dc_3 * ak3[ind] + dc_4 *
                            ak4[ind] + dc_5 * ak5[ind] + dc_6 * ak6[ind]))
        ind = ind + 1

    return [yout, yerr]


def rkqs(htry, eps, yscal, sp_comp, ks_dict, conc, r_dict, p_dict,
         stch_var, tvar, n_sp, molar):
    """This function is a helper function for tau-adaptive RK4"""
    safety = 0.9
    pgrow = -0.2
    pshrnk = -0.25
    errcon = 1.89e-4
    errmax = 1000
    hvar = htry
    while errmax > 1:
        ytemp, yerr = rkck(hvar, sp_comp, ks_dict, conc, r_dict, p_dict,
                           stch_var, tvar, n_sp, molar)
        errmax = max(0, np.max(np.array(yerr) / np.array(yscal))) / eps
        if errmax > 1:
            hvar = safety * hvar * (errmax**pshrnk)
            if hvar < 0.1 * hvar:
                hvar = 0.1 * hvar
        else:
            if errmax > errcon:
                hnext = safety * hvar * (errmax**pgrow)
            else:
                hnext = 5.0 * hvar
            hdid = hvar
            conc = ytemp

    return [conc, hdid, hnext, yerr]


def rungek4a_int(tvar, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                 yscal=10, molar=False, implicit=False, rfile=""):
    """[summary]

    Args:
        tvar (list): time stamp of trajectories i.e. [0, 0.1, 0.2, ...]
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
        ks_dict (dict): dictionary of rate constant that appears in each
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
        yscal (int, optional): [description]. Defaults to 10.
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.
        implicit (bool, optional): True means report in time intervals
            similar to the input time intervals even if actual step is
            more or less. Defaults to False.
        rfile (string, optional): name of topology file where some
            parameters or components are negative indicating  they  have
            to be estimated. Defaults to None.

    Returns:
        list: [time stamp, trajectory]
    """
    get_globals(rfile)
    tnew = []
    delt = tvar[-1] - tvar[-2]
    eps = 1.0e-8
    n_sp = [xvar for xvar in conc if xvar not in sp_comp]
    yconc = {xvar: conc[xvar] for xvar in sp_comp}
    slabels = [a for a in sp_comp]
    apply_rules(conc, yconc, [0], [conc[a] for a in sp_comp], slabels)
    yvar = [conc[a] for a in sp_comp]
    s_list = [yvar]
    if not implicit:
        tnow = tvar[0]
        tnew.append(tnow)
        tindex = 1

        while abs((tnew[-1] - tvar[-1]) / tvar[-1]) > 1.0e-5:

            y_old = conc
            conc, dtime, delt, _ = rkqs(
                delt, eps, yscal, sp_comp, ks_dict, conc, r_dict,
                p_dict, stch_var, tnow, n_sp, molar)
            tnow = tnow + dtime
            apply_rules(conc, yconc)
            if tnow > tvar[tindex]:
                tnow = tnow - dtime
                dtime = tvar[tindex] - tnow
                conc, tnow = runge_kutta_forth(
                    sp_comp, ks_dict, y_old, r_dict, p_dict, stch_var,
                    tnow, dtime, n_sp, molar)
                delt = 5 * dtime
                apply_rules(conc, yconc)
                tindex = tindex + 1
            s_list.append([conc[a] for a in sp_comp])
            tnew.append(tnow)
    else:
        tnow = tvar[0]
        tnew.append(tnow)
        tindex = 1

        while abs((tnew[-1] - tvar[-1]) / tvar[-1]) > 1.0e-5:

            y_old = conc
            conc, dtime, delt, _ = rkqs(
                delt, eps, yscal, sp_comp, ks_dict, conc, r_dict, p_dict,
                stch_var, tnow, n_sp, molar)
            tnow = tnow + dtime
            apply_rules(conc, yconc)
            if tnow > tvar[tindex]:
                tnow = tnow - dtime
                dtime = tvar[tindex] - tnow
                conc, tnow = runge_kutta_forth(
                    sp_comp, ks_dict, y_old, r_dict, p_dict, stch_var,
                    tnow, dtime, n_sp, molar)
                delt = 5 * dtime
                apply_rules(conc, yconc)
                s_list.append([conc[a] for a in sp_comp])
                tnew.append(tnow)
                tindex = tindex + 1
    return [tnew, np.array(s_list)]
