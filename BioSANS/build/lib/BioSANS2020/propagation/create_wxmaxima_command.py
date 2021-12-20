"""

            This module is the create_wxmaxima_command module

This module attempts to transforms topology file into wxmaxima ODE input
that a user can copy and paste into wxmaxima to grab a symbolic solution


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from sympy import Function, Symbol, Matrix

from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.propagation.propensity import propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals
# from BioSANS2020.myglobal import mglobals as globals2


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


def prepare_dict1(sp_comp):
    """This function prepares cs_var, cso_var, time

    Args:
        sp_comp (dict): dictionary of components position in reaction

    Returns:
        multitype: cs_var, cso_var, time
    """
    cs_var = {}
    cso_var = {}
    time = Symbol('t', real=True, positive=True)
    for xvar in sp_comp:
        cs_var[xvar] = Function(xvar, positive=True)
        cs_var[xvar] = cs_var[xvar](time)
        cso_var[xvar] = Symbol(xvar + "o", real=True, negative=False)
    return cs_var, cso_var, time


def prepare_list1(ks_dict):
    """This function prepares a list of symbol for rate constant

    Args:
        ks_dict (dict): dictionary of rate constant positions

    Returns:
        list: rate constants symbols
    """
    kcs = []
    for i, _ in enumerate(ks_dict):
        row = []
        if len(ks_dict[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
        kcs.append(row)
    return kcs


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


def for_wxmaxima(sp_comp, ks_dict, conc, r_dict, p_dict,
                 stch_var, items=None, rfile=""):
    """This function  transforms topology file  into wxmaxima  ODE input
    that a user can copy and paste into wxmaxima to grab a symbolic
    solution.

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
        items (list): list of [canvas, scroll_x, scroll_y]
        rfile (string, optional): name of topology file where some
            parameters or components are negative indicating  they  have
            to be estimated. Defaults to None.

    Returns:
        list: [0, 0] - contains nothing
    """
    get_globals(rfile)
    ffprint = prepare_ffrrint(items)

    ffprint(["/*Copy and paste to wxmaxima and run the cell*/\n\n"])
    cs_var, cso_var, time = prepare_dict1(sp_comp)
    kcs = prepare_list1(ks_dict)

    prop_flux = Matrix(propensity_vec_molar(kcs, cs_var, r_dict, p_dict, True))
    stch_var = Matrix(stch_var)
    # cs_var might have change after call
    for xvar in sp_comp:
        cs_var[xvar] = Function(xvar, positive=True)
        cs_var[xvar] = cs_var[xvar](time)
    slabels = list(cs_var.keys())  # [xvar for xvar in cs_var]

    ss_var, nz_var = prepare_list2(stch_var, slabels)

    da_dt = ss_var * prop_flux
    ccs = [cs_var[xvar] for xvar in cs_var]
    ccso = [cso_var[xvar] for xvar in cso_var]
    js_var = [ccs[xvar] for xvar in nz_var]
    jso = [ccso[xvar] for xvar in nz_var]

    ds_dt2 = []
    for xvar in js_var:
        ds_dt2.append(xvar.diff(time))

    fe_var = []
    for i, _ in enumerate(ds_dt2):
        ffprint([
            "f" + str(i + 1) + ":" + str(ds_dt2[i])
            .replace("Derivative", "diff") + " = " + str(da_dt[i]) + ";",
            "\n"])

        fe_var.append("f" + str(i + 1))
    ffprint(["\n"])
    for xvar, _ in enumerate(js_var):
        ffprint([
            "atvalue(" + str(js_var[xvar]) + ",t=0," + str(jso[xvar]) + ");",
            "\n"])

    ffprint(["\ndesolve(" + str(fe_var).replace('"', '').replace("'", '') +
             "," + str(js_var) + ");"])
    ffprint(["\n"])
    ffprint(["\n\n/*Copy and paste to wxmaxima and run the cell*/\n"])
    ffprint(["\n"])

    for xvar in sp_comp:
        cso_var[xvar] = conc[xvar]

    prop_flux = Matrix(
        propensity_vec_molar(ks_dict, cs_var, r_dict, p_dict, True))
    stch_var = Matrix(stch_var)
    # cs_var might have change after call
    for xvar in sp_comp:
        cs_var[xvar] = Function(xvar, positive=True)
        cs_var[xvar] = cs_var[xvar](time)
    slabels = list(cs_var.keys())  # [xvar for xvar in cs_var]
    ss_var = []
    nz_var = []
    for row in range(stch_var.shape[0]):
        if sum(abs(stch_var[row, :])) != 0 and slabels[row][0] != "-":
            ss_var.append(list(stch_var[row, :]))
            nz_var.append(row)
    ss_var = Matrix(ss_var)

    da_dt = ss_var * prop_flux
    ccs = [cs_var[xvar] for xvar in cs_var]
    ccso = [cso_var[xvar] for xvar in cso_var]
    js_var = [ccs[xvar] for xvar in nz_var]
    jso = [ccso[xvar] for xvar in nz_var]

    ds_dt2 = []
    for xvar in js_var:
        ds_dt2.append(xvar.diff(time))

    fe_var = []
    for i, _ in enumerate(ds_dt2):
        ffprint([
            "f" + str(i + 1) + ":" + str(ds_dt2[i])
            .replace("Derivative", "diff") + " = " + str(da_dt[i]) + ";",
            "\n"])
        fe_var.append("f" + str(i + 1))
    ffprint(["\n"])
    for xvar, _ in enumerate(js_var):
        ffprint([
            "atvalue(" + str(js_var[xvar]) + ",t=0," + str(jso[xvar]) + ");",
            "\n"])

    ffprint(["\ndesolve(" + str(fe_var).replace('"', '').replace("'", '') +
             "," + str(js_var) + ");"])

    return [0, 0]
