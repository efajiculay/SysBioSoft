"""

                    This is the euler_mod module

The main  task  this module  do is to  integrate  ordinary  differential
equation (ODE)  using the Euler method. In  this module, two implementa-
tion is provided. Both are tau-adaptive.

The list of functions in this module are as follows;

1. euler_model
2. euler_int
3. euler2_model
4. euler_wer_est
5. euler_help
6. euler2_int

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals, \
    apply_rules
from BioSANS2020.propagation.deterministic.lna_approx import \
    lna_steady_state
from BioSANS2020.myglobal import mglobals as globals2


def euler_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                stch_var, d_time, del_coef, molar=False):
    """This function returns the ODE model as a list of evaluated deri-
    vative of each components at a particular intance given by the state
    of the inputs.

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

            stch_var = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        d_time ([type]): [description]
        del_coef (float, optional): factor for modifying time steps used
            in the integration/propagation of ODE. Defaults to 10.
        molar : True if concentration or amount is in molar otherwise it
            is  False. Defaults to  False. If molar is True, macroscopic
            propensity is used. If it is False, microscopic propensity
            is used.

    Returns:
        list: 2 element list i.e. [value of dx/dt, d_time] where x are
            the components
    """

    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    fofx = np.matmul(stch_var, prop_flux)
    h1_div = np.abs(fofx) + 1.0e-30
    d_time = min(del_coef * np.min(1 / h1_div), d_time)
    dxdt = fofx * d_time
    if globals2.CON_BOUNDARY:
        spc = list(sp_comp.keys())  # [xvar for xvar in sp_comp]
        for xvar in globals2.CON_BOUNDARY:
            ind = spc.index(xvar)
            dxdt[ind] = 0
    return [dxdt.reshape(len(sp_comp)), d_time]


def euler_int(
        t_var, sp_comp, ks_dict, sconc, r_dict, p_dict, stch_var,
        del_coef=10, lna_solve=False, items=None, implicit=False,
        molar=False, rfile=""):
    """This function performs tau-adpative euler integration. The tau is
    adjusted in such a way that limit the change of the fastest reaction
    to del_coef amounts.

    Args:
        t_var (list): list of time points in the simulation
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
        stoch_var (numpy.ndarray): stoichiometric matrix. For example

            stoch_var = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        del_coef (float, optional): factor for modifying time steps used
            in the integration/propagation of ODE. Defaults to 10.
        lna_solve (bool, optional): if True, proceed to linar noise
            approximation calculations. Defaults to False.
        items (list, optional): 3 item list of [canvas, scroll_x,
            scroll_y], Defaults to None.
        implicit (bool, optional): True means report in time intervals
            similar to the input time intervals even if actual step is
            more or less. Defaults to False.
        molar (bool, optional): True if concentration or amount is in
            molar and you want to use macroscopic equations otherwise it
            should be False and using microscopic equations. Defaults to
            False.
        rfile (str): file name of BioSANS topology file.

    Returns:
        tuple: (time, trajectory)
    """
    get_globals(rfile)
    tnew = []
    d_time = t_var[-1] - t_var[-2]
    yconc = {xvar: sconc[xvar] for xvar in sconc}
    conc = {xvar: sconc[xvar] for xvar in sconc}
    apply_rules(conc, yconc)
    straj_list = [[conc[z] for z in sp_comp]]

    if not implicit:
        tnow = t_var[0]
        tnew.append(tnow)
        while tnow < t_var[-1]:
            mvar = np.nan_to_num(euler_model(
                sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                d_time, del_coef, molar))
            mm_var = mvar[0].reshape(1, len(mvar[0]))[0]
            tnow = tnow + mvar[1]
            ind = 0
            for spi in sp_comp:
                conc[spi] = straj_list[-1][ind] + mm_var[ind]
                ind = ind + 1
            apply_rules(conc, yconc)
            straj_list.append([conc[z] for z in sp_comp])
            tnew.append(tnow)
        if lna_solve:
            return lna_steady_state(
                t_var, sp_comp, ks_dict, conc, r_dict, p_dict,
                stch_var, items=items)
    else:
        tnow = t_var[0]
        tnew.append(tnow)
        tindex = 1
        conc_list = [conc[z] for z in sp_comp]
        while tnew[-1] < t_var[-1]:
            mvar = np.nan_to_num(euler_model(
                sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                d_time, del_coef, molar))
            tnow = tnow + mvar[1]
            if tnow > t_var[tindex]:
                tnow = tnow - mvar[1]
                dt2 = t_var[tindex] - tnow
                mvar = np.nan_to_num(euler_model(
                    sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                    dt2, del_coef, molar))
                mm_var = mvar[0].reshape(1, len(mvar[0]))[0]
                tnow = tnow + mvar[1]

                ind = 0
                for spi in sp_comp:
                    conc[spi] = conc_list[ind] + mm_var[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                conc_list = [conc[z] for z in sp_comp]
                for spi in sp_comp:
                    conc[spi] = max(0, conc[spi])
                straj_list.append(conc_list)

                tnew.append(tnow)
                tindex = tindex + 1
            else:
                mm_var = mvar[0].reshape(1, len(mvar[0]))[0]
                ind = 0
                for spi in sp_comp:
                    conc[spi] = conc_list[ind] + mm_var[ind]
                    ind = ind + 1
                apply_rules(conc, yconc)
                conc_list = [conc[z] for z in sp_comp]
                for spi in sp_comp:
                    conc[spi] = max(0, conc[spi])
    return (tnew, np.array(straj_list))


def euler2_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                 stch_var, molar=False):
    """[summary]

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
        stoch_var (numpy.ndarray): stoichiometric matrix. For example

            stoch_var = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        molar (bool, optional): True if concentration or amount is in
            molar and you want to use macroscopic equations otherwise it
            should be False and using microscopic equations. Defaults to
            False.

    Returns:
        np.ndarray: derivative of components with respect to time at a
        particular time or instant based on the current state of the
        inputs.
    """
    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(stch_var, prop_flux)
    if globals2.CON_BOUNDARY:
        spc = list(sp_comp.keys())  # [xvar for xvar in sp_comp]
        for xvar in globals2.CON_BOUNDARY:
            ind = spc.index(xvar)
            dxdt[ind] = 0
    return dxdt.reshape(len(sp_comp))


def euler_wer_est(h_var, sp_comp, ks_dict, conc, r_dict, p_dict,
                  stch_var, molar):
    """This function serves to do some part of the task in euler_help
    and euler2_int.
    """
    c_delta = h_var * euler2_model(sp_comp, ks_dict, conc, r_dict, p_dict,
                                   stch_var, molar)
    yconc = {}
    ind = 0
    for xvar in sp_comp:
        yconc[xvar] = conc[xvar] + c_delta[ind]
        ind = ind + 1
    err = 0.5 * h_var \
        * euler2_model(sp_comp, ks_dict, yconc, r_dict, p_dict,
                       stch_var, molar) - 0.5 * c_delta
    return [yconc, err]


def euler_help(htry, eps, yscal, sp_comp, ks_dict, conc, r_dict, p_dict,
               stch_var, molar):
    """This function serves to do some part of the task in euler2_int.
    """
    safety = 0.9
    pgrow = -0.2
    pshrnk = -0.25
    errcon = 1.89e-4
    errmax = 1000
    h_var = htry
    while errmax > 1:
        ytemp, yerr = euler_wer_est(
            h_var, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, molar)
        errmax = max(0, np.max(yerr / yscal)) / eps
        if errmax > 1:
            h_var = safety * h_var * (errmax**pshrnk)
            if h_var < 0.1 * h_var:
                h_var = 0.1 * h_var
        else:
            if errmax > errcon:
                hnext = safety * h_var * (errmax**pgrow)
            else:
                hnext = 5.0 * h_var
            hdid = h_var
            yvar = ytemp

    return [yvar, hdid, hnext, yerr]


def euler2_int(
        t_var, sp_comp, ks_dict, conc, r_dict, p_dict,
        stch_var, yscal=10, lna_solve=False, items=None, implicit=False,
        molar=False, rfile=""):
    """This function performs tau-adpative euler integration. The tau is
    adjusted to limit the error in each integration step as compared to
    what a second order runge-kutta would have predicted.

    Args:
        t_var (list): list of time points in the simulation
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
        stoch_var (numpy.ndarray): stoichiometric matrix. For example

            stoch_var = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        yscal (float, optional): factor for modifying time steps used
            in the integration/propagation of ODE. Defaults to 10.
        lna_solve (bool, optional): if True, proceed to linar noise
            approximation calculations. Defaults to False.
        items (list, optional): 3 item list of [canvas, scroll_x,
            scroll_y], Defaults to None.
        implicit (bool, optional): True means report in time intervals
            similar to the input time intervals even if actual step is
            more or less. Defaults to False.
        molar (bool, optional): True if concentration or amount is in
            molar and you want to use macroscopic equations otherwise it
            should be False and using microscopic equations. Defaults to
            False.
        rfile (str): file name of BioSANS topology file.

    Returns:
        tuple: (time, trajectory)
    """

    get_globals(rfile)
    tnew = []
    delt = t_var[-1] - t_var[-2]
    straj_list = [[conc[z] for z in sp_comp]]
    eps = 1.0e-4

    if not implicit:
        tnow = t_var[0]
        tnew.append(tnow)
        tindex = 1
        while tnew[-1] < t_var[-1]:
            conc_old = {xvar: conc[xvar] for xvar in conc}
            conc, d_time, delt, _ = euler_help(
                delt, eps, yscal, sp_comp, ks_dict, conc, r_dict, p_dict,
                stch_var, molar)
            tnow = tnow + d_time
            if tnow > t_var[tindex]:
                tnow = tnow - d_time
                d_time = t_var[tindex] - tnow
                conc, _ = euler_wer_est(
                    d_time, sp_comp, ks_dict, conc_old, r_dict, p_dict,
                    stch_var, molar)
                delt = 5 * d_time
                tindex = tindex + 1
                tnow = tnow + d_time
            straj_list.append([conc[z] for z in sp_comp])
            tnew.append(tnow)

        if lna_solve:
            return lna_steady_state(
                t_var, sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                items=items)
    else:
        tnow = t_var[0]
        tnew.append(tnow)
        tindex = 1
        while tnew[-1] < t_var[-1]:
            conc_old = {xvar: conc[xvar] for xvar in conc}
            conc, d_time, delt, _ = euler_help(
                delt, eps, yscal, sp_comp, ks_dict,
                conc, r_dict, p_dict, stch_var, molar)
            tnow = tnow + d_time
            if tnow > t_var[tindex]:
                tnow = tnow - d_time
                d_time = t_var[tindex] - tnow
                conc, _ = euler_wer_est(
                    d_time, sp_comp, ks_dict, conc_old, r_dict, p_dict,
                    stch_var, molar)
                delt = 5 * d_time
                straj_list.append([conc[z] for z in sp_comp])
                tnow = tnow + d_time
                tnew.append(tnow)
                tindex = tindex + 1
    return (tnew, np.array(straj_list))
