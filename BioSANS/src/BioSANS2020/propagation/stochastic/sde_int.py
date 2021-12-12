#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from scipy.integrate import odeint
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
import sdeint
from BioSANS2020.myglobal import mglobals as globals2


def sde_fx(z, t, Sp, ks_dict, r_dict, p_dict, V, molar=False):
    Spc = [s for s in Sp]
    conc = {Spc[x]: z[x] for x in range(len(Spc))}
    if not molar:
        D = propensity_vec(ks_dict, conc, r_dict, p_dict)
    else:
        D = propensity_vec_molar(ks_dict, conc, r_dict, p_dict)

    dxdt = np.matmul(V, D)
    for x in globals2.CON_BOUNDARY:
        ind = Spc.index(x)
        dxdt[ind] = 0
    return dxdt.T[0]


def sde_gx(z, t, Sp, ks_dict, r_dict, p_dict, V, molar=False):
    Spc = [s for s in Sp]
    conc = {Spc[x]: z[x] for x in range(len(Spc))}
    if not molar:
        D = propensity_vec(ks_dict, conc, r_dict, p_dict)
    else:
        D = propensity_vec_molar(ks_dict, conc, r_dict, p_dict)
    G = np.sqrt(D)
    #nlen = np.random.randn(len(D)).reshape(len(D),1)
    dxdt = np.matmul(V, G)
    for x in globals2.CON_BOUNDARY:
        ind = Spc.index(x)
        dxdt[ind] = 0
    return dxdt


def SDE_int(conc, t, Sp, ks_dict, r_dict, p_dict, V, molar=False, ito=True):
    z = [conc[a] for a in Sp]
    if ito:
        return sdeint.itoint(
            lambda z, t: sde_fx(z, t, Sp, ks_dict, r_dict, p_dict, V, molar),
            lambda z, t: sde_gx(z, t, Sp, ks_dict, r_dict, p_dict, V, molar), z, t
        )
    else:
        return sdeint.stratint(
            lambda z, t: sde_fx(z, t, Sp, ks_dict, r_dict, p_dict, V, molar),
            lambda z, t: sde_gx(z, t, Sp, ks_dict, r_dict, p_dict, V, molar), z, t
        )
