"""

            This module is the law_of_localization module

The sole purpose of this module is to implement the law of localization


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from sympy import Matrix, Symbol, simplify, NonSquareMatrixError


def subs2(zvar, cval):
    """This function substitute actual values to sympy symbols.

    Args:
        zvar (sympy.core): sympy expression
        cval (disct): dictionary of Symbols : value

    Returns:
        (sympy.core): substituted expression
    """
    for xvar in cval:
        zvar = zvar.subs(xvar, cval[xvar])
    return zvar


def prepare_ffrrint(items):
    """This function instantiate how the output will be printed by crea-
    thing the function ffprint

    Args:
        text (Text): text area where the outputs are written
        items (tuple): 3 item list of (canvas, scroll_x, scroll_y)

    Returns:
        function: either print in console or in text area
    """
    if items:
        text = prepare_scroll_text(items)

        def ffprint(xvar):
            return text.insert('insert', " ".join([str(y) for y in xvar]))
    else:

        def ffprint(xvar):
            return print(" ".join([str(y) for y in xvar]), end="")
    return ffprint


def prepare_dict_list1(sp_comp, conc, ks_dict):
    """This function prepares cs_var, cso_var, time

    Args:
        sp_comp (dict): dictionary of components position in reaction
        conc (dict) : dictionary of initail concentration

    Returns:
        (dict/list): cs_var, cso_var, equivals, equi_ks
    """
    cs_var = {}
    cso_var = {}
    equivals = []
    equi_co = []
    equi_ks = []

    for xvar in sp_comp:
        cs_var[xvar] = Symbol(xvar, real=True, negative=False)
        cso_var[xvar] = Symbol(xvar + 'o', real=True, negative=False) * \
            (0 if conc[xvar] == 0 else 1)
        equivals.append((cso_var[xvar], conc[xvar]))
        equi_co.append((cso_var[xvar], conc[xvar]))

    kcs = []
    for i, _ in enumerate(ks_dict):
        row = []
        if len(ks_dict[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[0], ks_dict[i][0]))
            equi_ks.append((row[0], ks_dict[i][0]))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[0], ks_dict[i][0]))
            equi_ks.append((row[0], ks_dict[i][0]))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[1], ks_dict[i][1]))
            equi_ks.append((row[1], ks_dict[i][1]))
        kcs.append(row)
    return cs_var, equivals, equi_co, equi_ks, kcs


def prepare_list2(stch_var, slabels):
    """This function prepares two list from stoichiometric matrix

    Args:
        stch_var (np.ndarray): toichiometric matrix

    Returns:
        list: [non zero row values of  stch_var,
            nonzero rows index of stch_var]
    """
    ss_var = []
    nz_var = []
    for row in range(stch_var.shape[0]):
        if sum(abs(stch_var[row, :])) != 0 and slabels[row][0] != "-":
            ss_var.append(list(stch_var[row, :]))
            nz_var.append(row)
    ss_var = Matrix(ss_var)
    return ss_var, nz_var


def law_loc_symbolic(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                     items=None, molar=False, mode=None, numer=False):
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
        stch_var (np.ndarray): stoichiometric matrix
        items (tuple): 3 item list of (canvas, scroll_x, scroll_y)
        molar (bool, optional): [description]. Defaults to False.
        mode (str, optional): either "Numeric", "fofks", "fofCo".
        numer (bool, optional): If True, numeric substitution will be
            done. Defaults to False.

    Returns:
        list: [0, 0] - useless return
    """

    ffprint = prepare_ffrrint(items)
    cs_var, equivals, equi_co, equi_ks, kcs \
        = prepare_dict_list1(sp_comp, conc, ks_dict)

    if not molar:
        prop_flux = Matrix(propensity_vec(kcs, cs_var, r_dict, p_dict))
    else:
        prop_flux = Matrix(propensity_vec_molar(kcs, cs_var, r_dict, p_dict))

    if mode == "Numeric":
        prop_flux = prop_flux.subs(equivals)
    elif mode == "fofks":
        prop_flux = prop_flux.subs(equi_co)
    elif mode == "fofCo":
        prop_flux = prop_flux.subs(equi_ks)

    stch_var = Matrix(stch_var)
    # cs_var might have change after call
    for xvar in sp_comp:
        cs_var[xvar] = Symbol(xvar, real=True, negative=False)
    slabels = list(cs_var.keys())

    ss_var, nz_var = prepare_list2(stch_var, slabels)

    ccs = [cs_var[xvar] for xvar in cs_var]
    js_name = [ccs[xvar] for xvar in nz_var]
    a_jac = prop_flux.jacobian(js_name)
    ker_v = ss_var.nullspace()

    as_jac = a_jac.col_insert(a_jac.shape[0], -ker_v[0])
    for ih_index in range(1, len(ker_v)):
        as_jac = as_jac.col_insert(as_jac.shape[0], -ker_v[ih_index])

    knames = []
    for i, _ in enumerate(kcs):
        if len(kcs[i]) == 1:
            knames.append(kcs[i][0])
        else:
            knames.append(kcs[i][0])
            knames.append(kcs[i][1])

    ffprint([
        "For inhibition or downregulation : \n\nWhen k decrease, the \
            sign of the sensitivity tells what happen to species activity\n"
    ])
    ffprint([
        "When k increase, the reverse of the sign of the sensitivity \
            tells what happen to species activity\n\n"
    ])

    if numer:
        as_jac = as_jac.subs(equivals)

    try:
        stch_var = simplify(as_jac.inv()).T
        if numer:
            stch_var = stch_var.subs(equivals)
        for ih_index in range(stch_var.shape[0]):
            for ij_index in range(stch_var.shape[1]):
                try:
                    ffprint([knames[ih_index], "\t", js_name[ij_index],
                             "\t", stch_var[ih_index, ij_index], "\n"])
                except BaseException:
                    pass
            ffprint(["\n"])
    except NonSquareMatrixError:
        ffprint(["Non square matrix error"])
    finally:
        pass

    return [0, 0]
