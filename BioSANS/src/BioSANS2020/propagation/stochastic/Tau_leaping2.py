#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propagation.propensity import *
from BioSANS2020.propagation.recalculate_globals import *


def Tau_leaping2(t, Sp, ks_dict, conc, r_dict, p_dict, V, rr,
                 del_coef=1, implicit=False, rfile=""):
    get_globals(rfile)
    tmax = t[-1]
    np.random.seed(int(rr * 100))
    tnew = []
    V2 = V**2
    stch_var = V.T
    dto = t[-1] - t[-2]
    tchlen = len(globals2.tCheck)

    AllSp = [z for z in Sp]
    Spc = [z for z in Sp if z not in reserve_events_words]
    Spc2 = [z for z in Sp if z in reserve_events_words]
    concz = {x: conc[x] for x in conc}
    yconc = {x: conc[x] for x in conc}
    apply_rules(concz, yconc)
    UpdateSp = [AllSp.index(z) for z in Spc]

    Vcri = np.where(V < 0)
    Vncr = np.array(list(range(len(V[0]))))
    Vncr = Vncr[~np.isin(Vncr, Vcri[1])]

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
            D = propensity_vec(ks_dict, concz, r_dict, p_dict)
            if len(Vcri[1]) > 0:
                Lcri = set()
                Lncr = set()
                for x in range(len(Vncr)):
                    Lncr.add(x)

                for x in range(len(Vcri[0])):
                    i, j = [Vcri[0][x], Vcri[1][x]]
                    if abs(Z[-1][i] / V[i, j]) < 10 and D[j] > 0:
                        Lcri.add(j)
                    else:
                        Lncr.add(j)
            else:
                Lcri = set()

            if len(Lcri) == len(V[0]):
                alp = np.sum(D)
                dt = (1 / alp) * (np.log(1 / np.random.uniform()))
            elif len(Lcri) == 0:
                alp = np.sum(D)
                uuj = np.matmul(V, D)
                sig = np.matmul(V2, D)
                exigi = np.array([Z[-1][j] * (1 / gi[j](Z[-1][j]))
                                  for j in range(len(Spc))]) * 0.03 * del_coef
                exigi1 = np.maximum(exigi, np.full(len(Spc), 1))
                dt = min(np.min(exigi1 / np.abs(uuj)),
                         np.min(exigi1 * exigi1 / np.abs(sig)))
            else:
                Lcri = np.array(list(Lcri))
                Lncr = np.array(list(Lncr))
                alpc = np.sum(D[Lcri])
                dtc = (1 / alpc) * (np.log(1 / np.random.uniform()))

                alp = np.sum(D[Lncr])
                uuj = np.matmul(V[:, Lncr], D[Lncr])
                sig = np.matmul(V2[:, Lncr], D[Lncr])
                exigi = np.array([Z[-1][j] * (1 / gi[j](Z[-1][j]))
                                  for j in range(len(Spc))]) * 0.03 * del_coef
                exigi1 = np.maximum(exigi, np.full(len(Spc), 1))
                dt = min(np.min(exigi1 / np.abs(uuj)),
                         np.min(exigi1 * exigi1 / np.abs(sig)), dtc)

            K = np.random.poisson(D * dt)
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
        index = 0
        C = [concz[z] for z in Spc]
        while t[tindex] < tmax:
            D = propensity_vec(ks_dict, concz, r_dict, p_dict)
            # step 1
            if len(Vcri[1]) > 0:
                Lcri = set()
                Lncr = set()
                for x in range(len(Vncr)):
                    Lncr.add(x)

                for x in range(len(Vcri[0])):
                    i, j = [Vcri[0][x], Vcri[1][x]]
                    if abs(Z[-1][i] / V[i, j]) < 10 and D[j] > 0:
                        Lcri.add(j)
                    else:
                        Lncr.add(j)
            else:
                Lcri = set()

            # step 2
            epsilon = 0.03
            Lcri = np.array(list(Lcri))
            Lncr = np.array(list(Lncr))
            if len(Lncr) == 0:  # No non critical reactions
                dt1 = 1.0e+10
            else:
                alp = np.sum(D[Lncr])
                uuj = np.matmul(V[:, Lncr], D[Lncr])
                sig = np.matmul(V2[:, Lncr], D[Lncr])
                exigi = np.array([C[j] * (1 / gi[j](C[j]))
                                  for j in range(len(Spc))]) * epsilon * del_coef
                exigi1 = np.maximum(exigi, np.full(len(Spc), 1))
                dt1 = min(np.min(exigi1 / np.abs(uuj)),
                          np.min(exigi1 * exigi1 / np.abs(sig)))

            Allpos = False
            while Allpos == False:
                Kmul, dt, doSSA, alpo = step_3to5(D, Lcri, dt1)  # step_3to5
                dt = min(dto, dt)
                if doSSA:
                    tc, tindex, index, break_now, Allpos, C = SSA_support(
                        t, Sp, ks_dict, r_dict, p_dict, V, rfile, tindex, index, tc, Z, Spc, Spc2, concz, yconc, UpdateSp)
                    if break_now:
                        break
                else:
                    Gupdate = False
                    if index != tchlen:
                        if tc + dt > globals2.tCheck[index]:
                            dt = globals2.tCheck[index] - tc
                            index = index + 1
                            Gupdate = True

                    if tc + dt > t[tindex]:
                        dt = max(0, t[tindex] - tc)
                        K = np.round(np.random.poisson(D * dt)) * Kmul
                        Allpos = True
                        cc = {}
                        bb = np.sum(K * stch_var[:, UpdateSp], 0)
                        for x in range(len(Spc)):
                            holder = C[x] + bb[x]
                            if holder >= 0 or Spc[x] in globals2.MODIFIED:
                                cc[Spc[x]] = holder
                            else:
                                Allpos = False
                                dt1 = dt1 / 2
                                break
                        if Allpos:
                            for x in range(len(Spc)):
                                concz[Spc[x]] = cc[Spc[x]]
                            for x in range(len(Spc2)):
                                concz[Spc2[x]] = concz[Spc2[x]] + dt
                            apply_rules(concz, yconc)
                            Z.append([concz[x] for x in Spc])
                            C = Z[-1]
                            tc = tc + dt
                            tindex = tindex + 1
                        else:
                            if Gupdate:
                                # pass
                                index = index - 1
                    else:
                        dt = max(0, dt)
                        K = np.round(np.random.poisson(D * dt)) * Kmul
                        Allpos = True
                        cc = {}
                        bb = np.sum(K * stch_var[:, UpdateSp], 0)
                        for x in range(len(Spc)):
                            holder = C[x] + bb[x]
                            if holder >= 0 or Spc[x] in globals2.MODIFIED:
                                cc[Spc[x]] = holder
                            else:
                                Allpos = False
                                dt1 = dt1 / 2
                                break
                        if Allpos:
                            for x in range(len(Spc)):
                                concz[Spc[x]] = cc[Spc[x]]
                            for x in range(len(Spc2)):
                                concz[Spc2[x]] = concz[Spc2[x]] + dt
                            apply_rules(concz, yconc)
                            C = [concz[x] for x in Spc]
                            tc = tc + dt
                        else:
                            if Gupdate:
                                index = index - 1
                                # pass
            if tindex >= len(t):
                break
        if len(t) != len(Z):
            while len(t) != len(Z):
                Z.append(Z[-1])
        tnew = t
    return (tnew, np.array(Z))


def step_3to5(D, Lcri, dt1):
    # step 3
    Kmul = np.ones((len(D), 1))
    r1 = np.random.uniform()
    while r1 == 0:
        r1 = np.random.uniform()

    alpo = np.sum(D)
    doSSA = False
    if dt1 < 10 * (1 / alpo):
        dt = (1 / alpo) * (np.log(1 / r1))
        doSSA = True
    else:
        # step 4
        if len(Lcri) > 0:
            alpc = np.sum(D[Lcri])
            dt2 = (1 / alpc) * (np.log(1 / r1))
        else:
            dt2 = 1.0e+5
            #dt2 = dt1

        # step 5
        if dt1 < dt2:
            dt = dt1
            if len(Lcri) > 0:
                Kmul[Lcri] = 0
        else:
            dt = dt2
            if len(Lcri) > 0:
                Kmul[Lcri] = 0
                P = np.cumsum([d / alpc for d in D[Lcri]])
                r2 = np.random.uniform()
                while r2 == 0:
                    r2 = np.random.uniform()
                for i in range(len(P)):
                    if r2 <= P[i]:
                        Jc = i
                        Kmul[Lcri[Jc]] = 1
                        break
    if dt == 1.0e+5:
        doSSA = True
    return [Kmul, dt, doSSA, alpo]


def SSA_support(t, Sp, ks_dict, r_dict, p_dict, V, rfile="", tindex=0, index=0, tc=0, Zc=[
], Spc=None, Spc2=None, concz=None, yconc=None, UpdateSp=None):
    get_globals(rfile)
    stch_var = V.T
    Z = [[concz[z] for z in Spc]]
    tmax = t[-1]

    break_now = False
    tchlen = len(globals2.tCheck)
    Allpos = True
    for ssa in range(100):
        if tc < tmax:
            D = propensity_vec(ks_dict, concz, r_dict, p_dict)
            alp = np.sum(D)
            r1 = np.random.uniform()
            while r1 == 0:
                r1 = np.random.uniform()
            P = np.cumsum([d / alp for d in D])
            dt = (1 / alp) * (np.log(1 / r1))

            #Gupdate = False
            if index != tchlen:
                if tc + dt >= globals2.tCheck[index]:
                    dt = globals2.tCheck[index] - tc
                    index = index + 1
                    #Gupdate = True

            if np.isnan(dt) or np.isinf(dt):
                Zc.append(Z[-1])
                while t[tindex] != tmax:
                    Zc.append(Z[-1])
                    tindex = tindex + 1
                tc = t[-1]
                break_now = True
                break

            r2 = np.random.uniform()
            for i in range(len(P)):
                if r2 <= P[i]:
                    Allpos = True
                    for x in range(len(Spc)):
                        holder = Z[-1][x] + stch_var[i][UpdateSp[x]]
                        if holder >= 0 or Spc[x] in globals2.MODIFIED:
                            concz[Spc[x]] = holder
                        else:
                            Allpos = False
                            break
                    if Allpos:
                        for x in range(len(Spc2)):
                            concz[Spc2[x]] = concz[Spc2[x]] + dt
                        apply_rules(concz, yconc)
                        Z.append([concz[x] for x in Spc])
                        tc = tc + dt
                        try:
                            if tc == t[tindex]:
                                Zc.append(Z[-1])
                                tindex = tindex + 1
                            else:
                                while tc > t[tindex]:
                                    Zc.append(Z[-2])
                                    tindex = tindex + 1
                        except BaseException:
                            pass
                    else:
                        for x in range(len(Spc)):
                            concz[Spc[x]] = Z[-1][x]
                        # if Gupdate:
                            #index = index - 1
                    break
        else:
            break_now = True
            break
    return [tc, tindex, index, break_now, Allpos, Z[-1]]
