"""

                  This module in the ode_int module

This module uses the odeint integrator in python to propagate ODE trajec
tory.
"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from numpy import matmul as np_matmul
from scipy.integrate import odeint
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals
from BioSANS2020.myglobal import mglobals as globals2


def ode_model(zlist, tvar, sp_comp, ks_dict, r_dict, p_dict,
              stch_var, molar=False):
    """This function returns dx/dt where x are the components or species
    and t is time.

    Args:
        zlist (list): list of components or species amounts
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
        molar : True if conc. is in molar else False

    Returns:
        np.ndarray: dx/dt where x are the components and t is time
    """
    spc = [s for s in sp_comp]
    conc = {spc[xvar]: zlist[xvar] for xvar in range(len(spc))}
    if not molar:
        prop_flux = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        prop_flux = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np_matmul(stch_var, prop_flux).reshape(len(zlist))
    for xvar in globals2.CON_BOUNDARY:
        ind = spc.index(xvar)
        dxdt[ind] = 0

    return dxdt


def ode_int(conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
            stch_var, molar=False, rfile=""):
    """This function returns the inetegration result of odeint

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
        stoch_var (numpy.ndarray): stoichiometric matrix. For example

            stoch_var = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.
        rfile (str): file name of BioSANS topology file.

    Returns:
        np.ndarray: trajectories
    """
    get_globals(rfile)
    zlist = [conc[a] for a in sp_comp]
    return odeint(ode_model, zlist, tvar,
                  args=(sp_comp, ks_dict, r_dict, p_dict, stch_var, molar))
