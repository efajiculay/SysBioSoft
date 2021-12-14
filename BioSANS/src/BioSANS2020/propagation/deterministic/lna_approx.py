"""

                     This is the lna_approx module

     This module performs numerical linear noise approximation or LNA by
exploiting the following relationship.

                         AC + CA.T + BB = 0

where A is defined as d(S*f)/dx where S is the stoichiometric  matrix, f
are the propensities or fluxes, and x are the components or species. The
flux is a function of  x and rate constant k. BB is the diffusion matrix
equivalent  to S*diag(f)*S.T where diag(f) is a  square matrix with zero
non-diagonal elements and f[i] in each diagonal elements.

The following are the list of function for this module.

1. rem_rowcol_zero
2. lna_ss_jacobian
3. lna_model_ss
4. lna_steady_state


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from scipy import linalg as LA
from scipy.optimize import fsolve

from BioSANS2020.propagation.propensity import propensity_vec_molar
from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text


def rem_rowcol_zero(a_mat):
    """This function removes rows and columns without non-zero entries.

    Args:
        a_mat (np.ndarray): numpy matrix of A or d(S*f)/dx as described
            in the module docstring.

    Returns:
        np.ndarray: numpy matrix with no rows and columns without non-
            zero entries.
    """
    return a_mat[:, ~np.all(a_mat == 0, axis=0)][~np.all(a_mat == 0, axis=1)]


def lna_ss_jacobian(model, zlist, sp_comp, stch_var,
                    ks_dict, r_dict, p_dict):
    """This  function calculataes  the  jacobian  of  the model ODE with
    respect to the list of species concentration zlist.

    Args:
        model (function): the ODE model returning derivative of species
            or components as a function of time
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
        stch_var (np.ndarray): stoichiometric matrix
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

    Returns:
        np.ndarray: jacobian matrix or d(S*f)/dx or A
    """
    jacbn = np.zeros((len(zlist), len(zlist)))
    div = 1
    new = 1000
    old = 2000
    check = 1000
    while check > 1.0e-7:
        for i, _ in enumerate(zlist):
            h_var = abs(zlist[i] * 1.0e-8) / div
            h2_var = 2 * h_var
            zlist[i] = zlist[i] - h_var
            ini_val = model(zlist, sp_comp, ks_dict, r_dict,
                            p_dict, stch_var)
            zlist[i] = zlist[i] + 2 * h_var
            out_val = model(zlist, sp_comp, ks_dict, r_dict,
                            p_dict, stch_var)
            zlist[i] = zlist[i] - h_var
            for j, _ in enumerate(out_val):
                jacbn[j, i] = (out_val[j] - ini_val[j]) / h2_var
        div = div * 2
        old = new
        new = sum([jacbn[k, k] for k in range(len(zlist))])
        check = abs(new - old) / abs(new)
    return np.array(jacbn)


def lna_model_ss(zlist, sp_comp, ks_dict, r_dict, p_dict, stch_var):
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

    Returns:
        np.ndarray: derivative of species with respect to time.d(S*f)/dt
    """
    conc = {}
    ind = 0
    for spi in sp_comp:
        conc[spi] = zlist[ind]
        ind = ind + 1
    prop_flux = propensity_vec_molar(
        ks_dict, conc, r_dict, p_dict, True)
    fofx = np.matmul(stch_var, prop_flux)
    return fofx.reshape(len(sp_comp))


def lna_steady_state(t_var, sp_comp, ks_dict, conc, r_dict, p_dict,
                     stch_var, items=None):
    """[summary]

    Args:
        t_var (list): time stamp
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
        stch_var (np.ndarray): stoichiometric matrix

        items (tuplel): (canvas, scroll_x, scroll_y). Defaults to None.

    Returns:
        [type]: [description]
    """
    ind = 0
    zlist = []
    for spi in sp_comp:
        zlist.append(conc[spi])
    zlist = fsolve(
        lna_model_ss,
        tuple(zlist),
        xtol=1.0e-10,
        args=(sp_comp, ks_dict, r_dict, p_dict, stch_var))
    ind = 0
    for spi in sp_comp:
        conc[spi] = zlist[ind]
        ind = ind + 1
    dsf_dx = lna_ss_jacobian(
        lna_model_ss, zlist, sp_comp, stch_var, ks_dict, r_dict, p_dict)
    f_prop = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)
    bb_diffsn = np.matmul(
        np.matmul(
            stch_var,
            np.diag(
                f_prop.flatten())),
        stch_var.T)

    dsf_dx = np.nan_to_num(dsf_dx)
    bb_diffsn = np.nan_to_num(bb_diffsn)
    cov_mat = LA.solve_continuous_lyapunov(dsf_dx, -bb_diffsn)

    if items:
        text = prepare_scroll_text(items)

        def fprint(xvar):
            return text.insert('insert', " ".join([str(y) for y in xvar]))
    else:
        def fprint(xvar):
            return print(" ".join([str(y) for y in xvar]), end="")

    fprint(["\nConcentrations\n\n"])
    ind = 0
    for spi in sp_comp:
        if spi[0] != "-":
            fprint([spi, " = ", zlist[ind], "\n"])
            ind = ind + 1

    fprint(["\nCovariance\n\n"])
    i_ind = 0
    for spi in sp_comp:
        j = 0
        for spj in sp_comp:
            if j >= i_ind:
                val = cov_mat[i_ind, j]
                if str(val) not in {"None", "nan", "0.0"} \
                        and spi[0] != "-" and spj[0] != "-":
                    fprint([" ".join(["Covr", spi, spj])
                            .ljust(50), "=", val, "\n"])
            j = j + 1
        i_ind = i_ind + 1
    fprint(["\nCorrelation\n\n"])
    i_ind = 0
    for spi in sp_comp:
        j = 0
        for spj in sp_comp:
            if j >= i_ind:
                val = cov_mat[i_ind, j] / \
                    np.sqrt(np.abs(cov_mat[i_ind, i_ind]
                                   * cov_mat[j, j]))
                if str(val) not in {"None", "nan", "0.0"} \
                        and spi[0] != "-" and spj[0] != "-":
                    fprint([" ".join(["Corr", spi, spj]).ljust(50),
                            "=", val, "\n"])
            j = j + 1
        i_ind = i_ind + 1

    fprint(["\nFano Factor\n\n"])
    ind = 0
    for spi in sp_comp:
        val = cov_mat[ind, ind] / zlist[ind]
        if str(val) not in {"None", "nan", "0.0"} and spi[0] != "-":
            fprint([" ".join(["Fano Factor for", spi])
                    .ljust(50), "=", val, "\n"])
        ind = ind + 1

    return [cov_mat, zlist]
