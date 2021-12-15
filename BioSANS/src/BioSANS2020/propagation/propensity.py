# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.myglobal import mglobals as globals2
# from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT


def propensity_vec(ks_dict, conc, r_dict, p_dict, odeint=False):
    # this is for microscopic
    prop_flux = []
    rxn = len(ks_dict)

    if odeint:
        for x in globals2.MODIFIED:
            spvar, pfunc = globals2.MODIFIED[x][0]
            try:
                sprv = []
                for c in range(len(spvar)):
                    sprv.append(conc[spvar[c].strip()])
                conc[x] = pfunc(*sprv)
            except BaseException:
                conc[x] = pfunc()

    for r in range(rxn):
        key = "Prop_" + str(r)
        if key in globals2.PROP_MODIFIED:
            rowProp = globals2.PROP_MODIFIED[key]
            for row in rowProp:
                try:
                    spvar, pfunc = row
                    sprv = []
                    for c in range(len(spvar)):
                        sprv.append(conc[spvar[c].strip()])
                    prop_flux.append(pfunc(*sprv))
                except BaseException:
                    prop_flux.append(pfunc())
        else:
            if len(r_dict[r]) == 1:
                for x in r_dict[r]:
                    if r_dict[r][x] == 1:
                        prop_flux.append(ks_dict[r][0] * conc[x])
                    elif r_dict[r][x] == 2:
                        prop_flux.append(
                            ks_dict[r][0]
                            * max(conc[x] * (conc[x] - 1), 0) / 2)
                    elif r_dict[r][x] == 0:
                        prop_flux.append(ks_dict[r][0])
            elif len(r_dict[r]) == 2:
                p = ks_dict[r][0]
                for x in r_dict[r]:
                    p = p * conc[x]
                prop_flux.append(p)
            if len(ks_dict[r]) == 2:
                if len(p_dict[r]) == 1:
                    for x in p_dict[r]:
                        if p_dict[r][x] == 1:
                            prop_flux.append(ks_dict[r][1] * conc[x])
                        elif p_dict[r][x] == 2:
                            prop_flux.append(
                                ks_dict[r][1]
                                * max(conc[x] * (conc[x] - 1), 0) / 2)
                        elif p_dict[r][x] == 0:
                            prop_flux.append(ks_dict[r][1])
                elif len(p_dict[r]) == 2:
                    p = ks_dict[r][1]
                    for x in p_dict[r]:
                        if p_dict[r][x] == 1:
                            p = p * conc[x]
                        else:
                            p = p * max(conc[x] * (conc[x] - 1), 0) / 2
                    prop_flux.append(p)
    try:
        return np.array(prop_flux).reshape(len(prop_flux), 1).astype(float)
    except BaseException:
        return np.array(prop_flux).reshape(len(prop_flux), 1)


def propensity_vec_molar(ks_dict, conc, r_dict, p_dict, odeint=False):
    # this is for macroscopic
    prop_flux = []
    rxn = len(ks_dict)

    if odeint:
        for x in globals2.MODIFIED:
            spvar, pfunc = globals2.MODIFIED[x][0]
            try:
                sprv = []
                for c in range(len(spvar)):
                    sprv.append(conc[spvar[c].strip()])
                conc[x] = pfunc(*sprv)
            except BaseException:
                suby = pfunc()
                conc[x] = suby

    for r in range(rxn):
        key = "Prop_" + str(r)
        if key in globals2.PROP_MODIFIED:
            rowProp = globals2.PROP_MODIFIED[key]
            for row in rowProp:
                try:
                    spvar, pfunc = row
                    sprv = []
                    for c in range(len(spvar)):
                        sprv.append(conc[spvar[c].strip()])
                    prop_flux.append(pfunc(*sprv))
                except BaseException:
                    prop_flux.append(pfunc())
        else:
            if len(r_dict[r]) == 1:
                for x in r_dict[r]:
                    if x not in conc:
                        conc[x] = 1
                    prop_flux.append(ks_dict[r][0] * conc[x]**r_dict[r][x])
            elif len(r_dict[r]) == 2:
                p = ks_dict[r][0]
                for x in r_dict[r]:
                    if x not in conc:
                        conc[x] = 1
                    p = p * conc[x]**r_dict[r][x]
                prop_flux.append(p)

            if len(ks_dict[r]) == 2:
                if len(p_dict[r]) == 1:
                    for x in p_dict[r]:
                        prop_flux.append(ks_dict[r][1] * conc[x]**p_dict[r][x])
                elif len(p_dict[r]) == 2:
                    p = ks_dict[r][1]
                    for x in p_dict[r]:
                        if x not in conc:
                            conc[x] = 1
                        p = p * conc[x]**p_dict[r][x]
                    prop_flux.append(p)

    try:
        return np.array(prop_flux).reshape(len(prop_flux), 1).astype(float)
    except BaseException:
        return np.array(prop_flux).reshape(len(prop_flux), 1)
