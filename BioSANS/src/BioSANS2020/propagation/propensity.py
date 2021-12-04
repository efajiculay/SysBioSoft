#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.math_functs.sbmlMath import *


def propensity_vec(Ks, conc, Rr, Rp, odeint=False):  # this is for microscopic
    D = []
    Rxn = len(Ks)

    if odeint:
        for x in globals2.modified:
            spvar, pfunc = globals2.modified[x][0]
            try:
                sprv = []
                for c in range(len(spvar)):
                    sprv.append(conc[spvar[c].strip()])
                conc[x] = pfunc(*sprv)
            except:
                conc[x] = pfunc()

    for r in range(Rxn):
        key = "Prop_"+str(r)
        if key in globals2.PropModified:
            rowProp = globals2.PropModified[key]
            for row in rowProp:
                try:
                    spvar, pfunc = row
                    sprv = []
                    for c in range(len(spvar)):
                        sprv.append(conc[spvar[c].strip()])
                    D.append(pfunc(*sprv))
                except:
                    D.append(pfunc())
        else:
            if len(Rr[r]) == 1:
                for x in Rr[r]:
                    if Rr[r][x] == 1:
                        D.append(Ks[r][0]*conc[x])
                    elif Rr[r][x] == 2:
                        D.append(Ks[r][0]*max(conc[x]*(conc[x]-1), 0)/2)
                    elif Rr[r][x] == 0:
                        D.append(Ks[r][0])
            elif len(Rr[r]) == 2:
                p = Ks[r][0]
                for x in Rr[r]:
                    p = p*conc[x]
                D.append(p)
            if len(Ks[r]) == 2:
                if len(Rp[r]) == 1:
                    for x in Rp[r]:
                        if Rp[r][x] == 1:
                            D.append(Ks[r][1]*conc[x])
                        elif Rp[r][x] == 2:
                            D.append(Ks[r][1]*max(conc[x]*(conc[x]-1), 0)/2)
                        elif Rp[r][x] == 0:
                            D.append(Ks[r][1])
                elif len(Rp[r]) == 2:
                    p = Ks[r][1]
                    for x in Rp[r]:
                        if Rp[r][x] == 1:
                            p = p*conc[x]
                        else:
                            p = p*max(conc[x]*(conc[x]-1), 0)/2
                    D.append(p)
    try:
        return np.array(D).reshape(len(D), 1).astype(float)
    except:
        return np.array(D).reshape(len(D), 1)


def propensity_vec_molar(Ks, conc, Rr, Rp, odeint=False):  # this is for macroscopic
    D = []
    Rxn = len(Ks)

    if odeint:
        for x in globals2.modified:
            spvar, pfunc = globals2.modified[x][0]
            try:
                sprv = []
                for c in range(len(spvar)):
                    sprv.append(conc[spvar[c].strip()])
                conc[x] = pfunc(*sprv)
            except:
                suby = pfunc()
                conc[x] = suby

    for r in range(Rxn):
        key = "Prop_"+str(r)
        if key in globals2.PropModified:
            rowProp = globals2.PropModified[key]
            for row in rowProp:
                try:
                    spvar, pfunc = row
                    sprv = []
                    for c in range(len(spvar)):
                        sprv.append(conc[spvar[c].strip()])
                    D.append(pfunc(*sprv))
                except:
                    D.append(pfunc())
        else:
            if len(Rr[r]) == 1:
                for x in Rr[r]:
                    if x not in conc:
                        conc[x] = 1
                    D.append(Ks[r][0]*conc[x]**Rr[r][x])
            elif len(Rr[r]) == 2:
                p = Ks[r][0]
                for x in Rr[r]:
                    if x not in conc:
                        conc[x] = 1
                    p = p*conc[x]**Rr[r][x]
                D.append(p)

            if len(Ks[r]) == 2:
                if len(Rp[r]) == 1:
                    for x in Rp[r]:
                        D.append(Ks[r][1]*conc[x]**Rp[r][x])
                elif len(Rp[r]) == 2:
                    p = Ks[r][1]
                    for x in Rp[r]:
                        if x not in conc:
                            conc[x] = 1
                        p = p*conc[x]**Rp[r][x]
                    D.append(p)

    try:
        return np.array(D).reshape(len(D), 1).astype(float)
    except:
        return np.array(D).reshape(len(D), 1)
