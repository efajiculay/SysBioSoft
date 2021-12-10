#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from scipy.integrate import odeint
from BioSANS2020.propagation.propensity import *
from BioSANS2020.propagation.recalculate_globals import *
from BioSANS2020.myglobal import mglobals as globals2


def ODE_model(z, t, Sp, ks_dict, r_dict, p_dict, V, molar=False):
    Spc = [s for s in Sp]
    conc = {Spc[x]: z[x] for x in range(len(Spc))}
    if not molar:
        D = propensity_vec(ks_dict, conc, r_dict, p_dict, True)
    else:
        D = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)

    dxdt = np.matmul(V, D).reshape(len(z))
    for x in globals2.CON_BOUNDARY:
        ind = Spc.index(x)
        dxdt[ind] = 0
    return dxdt


def ODE_int(conc, t, Sp, ks_dict, r_dict, p_dict, V, molar=False, rfile=""):
    get_globals(rfile)
    z = [conc[a] for a in Sp]
    return odeint(ODE_model, z, t, args=(Sp, ks_dict, r_dict, p_dict, V, molar))
