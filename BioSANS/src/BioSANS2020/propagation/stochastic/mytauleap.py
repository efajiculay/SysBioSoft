#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.propagation.recalculate_globals import get_globals


def Sim_TauLeap(t, Sp, ks_dict, conc, r_dict, p_dict, V, rr,
                del_coef, implicit=False, rfile=""):
    get_globals(rfile)
    tmax = t[-1]
    np.random.seed(int(rr * 100))
    tnew = []
    V2 = V**2
    stch_var = V.T

    AllSp = [z for z in Sp]
    Spc = [z for z in Sp if z not in reserve_events_words]
    Spc2 = [z for z in Sp if z in reserve_events_words]
    concz = {x: conc[x] for x in conc}
    yconc = {x: conc[x] for x in conc}
    apply_rules(concz, yconc)
    UpdateSp = [AllSp.index(z) for z in Spc]
    Z = [[concz[z] for z in Spc]]
    tc = 0
    tnew.append(tc)

    gi = []
    for i in range(len(Spc)):
        keys = Sp[Spc[i]]
        maxrlen = 0
        rlen = 0
        for key in keys:
            if Spc[i] in r_dict[key]:
                maxrlen = max(maxrlen, len(r_dict[key]))
            elif Spc[i] in p_dict[key]:
                maxrlen = max(maxrlen, len(p_dict[key]))
        if maxrlen == 1:
            for key in keys:
                if Spc[i] in r_dict[key]:
                    rlen = max(rlen, r_dict[key][Spc[i]])
                elif Spc[i] in p_dict[key]:
                    rlen = max(rlen, p_dict[key][Spc[i]])
            if rlen == 1:
                gi.append(lambda x: 1)
            else:
                gi.append(lambda x: 2)
        elif maxrlen == 2:
            try:
                gi.append(lambda x: 2 + 1 / (x - 1))
            except BaseException:
                gi.append(lambda x: 2)
        else:
            gi.append(lambda x: 1)

    if not implicit:
        while tc < tmax:
            for x in range(len(Spc)):
                concz[Spc[x]] = Z[-1][x]
            D = propensity_vec(ks_dict, concz, r_dict, p_dict)
            alp = np.sum(D)
            r1 = np.random.uniform()
            while r1 == 0:
                r1 = np.random.uniform()
            dt1 = (1 / alp) * (np.log(1 / r1))

            uuj = np.matmul(V, D)
            sig = np.matmul(V2, D)
            exigi = np.array([Z[-1][j] * (1 / gi[j](Z[-1][j]))
                              for j in range(len(Spc))]) * 0.03 * del_coef
            exigi1 = np.maximum(exigi, np.full(len(Spc), 1))
            dt2 = min(np.min(exigi1 / np.abs(uuj)),
                      np.min(exigi1 * exigi1 / np.abs(sig)))

            dt = max(dt1, dt2)
            K = np.round(np.random.poisson(D * dt))
            Allpos = True
            cc = {}
            bb = np.sum(K * stch_var[:, UpdateSp], 0)
            for x in range(len(Spc)):
                holder = Z[-1][x] + bb[x]
                if holder >= 0 or Spc[x] in globals2.MODIFIED:
                    cc[Spc[x]] = holder
                else:
                    Allpos = False
                    break
            if Allpos:
                for x in range(len(Spc)):
                    concz[Spc[x]] = cc[Spc[x]]
                for x in range(len(Spc2)):
                    concz[Spc2[x]] = concz[Spc2[x]] + dt
                apply_rules(concz, yconc)
                Z.append([concz[x] for x in Spc])
                tc = tc + dt
                tnew.append(tc)

    else:
        tindex = 1
        Zc = [Z[-1]]
        while tc < tmax:
            D = propensity_vec(ks_dict, concz, r_dict, p_dict)
            D[D < 0] = 0
            alp = np.sum(D)
            r1 = np.random.uniform()
            while r1 == 0:
                r1 = np.random.uniform()
            dt1 = (1 / alp) * (np.log(1 / r1))

            uuj = np.matmul(V, D)
            sig = np.matmul(V2, D)
            exigi = np.array([Z[-1][j] * (1 / gi[j](Z[-1][j]))
                              for j in range(len(Spc))]) * 0.03 * del_coef
            exigi1 = np.maximum(exigi, np.full(len(Spc), 1))
            dt2 = min(np.min(exigi1 / np.abs(uuj)),
                      np.min(exigi1 * exigi1 / np.abs(sig)))

            dt = max(dt1, dt2)
            if tc + dt > t[tindex]:
                dt = max(0, t[tindex] - tc)
                K = np.round(np.random.poisson(D * dt))
                Allpos = True
                cc = {}
                bb = np.sum(K * stch_var[:, UpdateSp], 0)
                for x in range(len(Spc)):
                    holder = Z[-1][x] + bb[x]
                    if holder >= 0 or Spc[x] in globals2.MODIFIED:
                        cc[Spc[x]] = holder
                    else:
                        Allpos = False
                        break
                if Allpos:
                    for x in range(len(Spc)):
                        concz[Spc[x]] = cc[Spc[x]]
                    for x in range(len(Spc2)):
                        concz[Spc2[x]] = concz[Spc2[x]] + dt
                    apply_rules(concz, yconc)
                    Z.append([concz[x] for x in Spc])
                    Zc.append(Z[-1])
                    tc = tc + dt
                    tindex = tindex + 1
            else:
                dt = max(0, dt)
                K = np.round(np.random.poisson(D * dt))
                Allpos = True
                cc = {}
                bb = np.sum(K * stch_var[:, UpdateSp], 0)
                for x in range(len(Spc)):
                    holder = Z[-1][x] + bb[x]
                    if holder >= 0 or Spc[x] in globals2.MODIFIED:
                        cc[Spc[x]] = holder
                    else:
                        Allpos = False
                        break
                if Allpos:
                    for x in range(len(Spc)):
                        concz[Spc[x]] = cc[Spc[x]]
                    for x in range(len(Spc2)):
                        concz[Spc2[x]] = concz[Spc2[x]] + dt
                    apply_rules(concz, yconc)
                    Z.append([concz[x] for x in Spc])
                    tc = tc + dt

        if len(t) != len(Zc):
            while len(t) != len(Zc):
                Zc.append(Z[-1])
        tnew = t
        Z = Zc
    return (tnew, np.array(Z))
