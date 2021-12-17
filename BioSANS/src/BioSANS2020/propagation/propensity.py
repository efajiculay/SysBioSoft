"""

                     This is the propensity module

This module prepares the propensity vector or fluxes vector


"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.myglobal import mglobals as globals2
# from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT


def propensity_vec(ks_dict, conc, r_dict, p_dict, odeint=False):
    """Returns propensity vector using microscopic equations

    Args:
        ks_dict (dict): check details in process module
        conc (dict): check details in process module
        r_dict (dict): check details in process module
        p_dict (dict): check details in process module
        odeint (bool, optional): [description]. Defaults to False.

    Returns:
        np.ndarray: propensity vector
    """
    # this is for microscopic
    prop_flux = []
    rxn = len(ks_dict)

    if odeint:
        for xvar in globals2.MODIFIED:
            spvar, pfunc = globals2.MODIFIED[xvar][0]
            try:
                sprv = []
                for cvar, _ in enumerate(spvar):
                    sprv.append(conc[spvar[cvar].strip()])
                conc[xvar] = pfunc(*sprv)
            except BaseException:
                conc[xvar] = pfunc()

    for rvar in range(rxn):
        key = "Prop_" + str(rvar)
        if key in globals2.PROP_MODIFIED:
            row_prop = globals2.PROP_MODIFIED[key]
            for row in row_prop:
                try:
                    spvar, pfunc = row
                    sprv = []
                    for cvar, _ in enumerate(spvar):
                        sprv.append(conc[spvar[cvar].strip()])
                    prop_flux.append(pfunc(*sprv))
                except BaseException:
                    prop_flux.append(pfunc())
        else:
            if len(r_dict[rvar]) == 1:
                for xvar in r_dict[rvar]:
                    if r_dict[rvar][xvar] == 1:
                        prop_flux.append(ks_dict[rvar][0] * conc[xvar])
                    elif r_dict[rvar][xvar] == 2:
                        prop_flux.append(
                            ks_dict[rvar][0]
                            * max(conc[xvar] * (conc[xvar] - 1), 0) / 2)
                    elif r_dict[rvar][xvar] == 0:
                        prop_flux.append(ks_dict[rvar][0])
            elif len(r_dict[rvar]) == 2:
                pvar = ks_dict[rvar][0]
                for xvar in r_dict[rvar]:
                    pvar = pvar * conc[xvar]
                prop_flux.append(pvar)
            if len(ks_dict[rvar]) == 2:
                if len(p_dict[rvar]) == 1:
                    for xvar in p_dict[rvar]:
                        if p_dict[rvar][xvar] == 1:
                            prop_flux.append(ks_dict[rvar][1] * conc[xvar])
                        elif p_dict[rvar][xvar] == 2:
                            prop_flux.append(
                                ks_dict[rvar][1]
                                * max(conc[xvar] * (conc[xvar] - 1), 0) / 2)
                        elif p_dict[rvar][xvar] == 0:
                            prop_flux.append(ks_dict[rvar][1])
                elif len(p_dict[rvar]) == 2:
                    pvar = ks_dict[rvar][1]
                    for xvar in p_dict[rvar]:
                        if p_dict[rvar][xvar] == 1:
                            pvar = pvar * conc[xvar]
                        else:
                            pvar = pvar * max(conc[xvar]
                                              * (conc[xvar] - 1), 0) / 2
                    prop_flux.append(pvar)
    try:
        return np.array(prop_flux).reshape(len(prop_flux), 1).astype(float)
    except BaseException:
        return np.array(prop_flux).reshape(len(prop_flux), 1)


def propensity_vec_molar(ks_dict, conc, r_dict, p_dict, odeint=False):
    """Returns propensity vector using macroscopic equations

    Args:
        ks_dict (dict): check details in process module
        conc (dict): check details in process module
        r_dict (dict): check details in process module
        p_dict (dict): check details in process module
        odeint (bool, optional): [description]. Defaults to False.

    Returns:
        np.ndarray: propensity vector
    """
    # this is for macroscopic
    prop_flux = []
    rxn = len(ks_dict)

    if odeint:
        for xvar in globals2.MODIFIED:
            spvar, pfunc = globals2.MODIFIED[xvar][0]
            try:
                sprv = []
                for cvar, _ in enumerate(spvar):
                    sprv.append(conc[spvar[cvar].strip()])
                conc[xvar] = pfunc(*sprv)
            except BaseException:
                suby = pfunc()
                conc[xvar] = suby

    for rvar in range(rxn):
        key = "Prop_" + str(rvar)
        if key in globals2.PROP_MODIFIED:
            row_prop = globals2.PROP_MODIFIED[key]
            for row in row_prop:
                try:
                    spvar, pfunc = row
                    sprv = []
                    for cvar, _ in enumerate(spvar):
                        sprv.append(conc[spvar[cvar].strip()])
                    prop_flux.append(pfunc(*sprv))
                except BaseException:
                    prop_flux.append(pfunc())
        else:
            if len(r_dict[rvar]) == 1:
                for xvar in r_dict[rvar]:
                    if xvar not in conc:
                        conc[xvar] = 1
                    prop_flux.append(ks_dict[rvar][0]
                                     * conc[xvar]**r_dict[rvar][xvar])
            elif len(r_dict[rvar]) == 2:
                pvar = ks_dict[rvar][0]
                for xvar in r_dict[rvar]:
                    if xvar not in conc:
                        conc[xvar] = 1
                    pvar = pvar * conc[xvar]**r_dict[rvar][xvar]
                prop_flux.append(pvar)

            if len(ks_dict[rvar]) == 2:
                if len(p_dict[rvar]) == 1:
                    for xvar in p_dict[rvar]:
                        prop_flux.append(ks_dict[rvar][1]
                                         * conc[xvar]**p_dict[rvar][xvar])
                elif len(p_dict[rvar]) == 2:
                    pvar = ks_dict[rvar][1]
                    for xvar in p_dict[rvar]:
                        if xvar not in conc:
                            conc[xvar] = 1
                        pvar = pvar * conc[xvar]**p_dict[rvar][xvar]
                    prop_flux.append(pvar)

    try:
        return np.array(prop_flux).reshape(len(prop_flux), 1).astype(float)
    except BaseException:
        return np.array(prop_flux).reshape(len(prop_flux), 1)
