#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.gui_functs.scrollable_text import *
from BioSANS2020.propagation.propensity import *
from sympy import *


def subs2(Z, cval):
    for x in cval:
        Z = Z.subs(x, cval[x])
    return Z


def LNA_symbolic(Sp, Ks, conc, Rr, Rp, V, items=None, molar=False, mode=None):

    Cs = {}
    Cso = {}
    equivals = []
    equiCo = []
    equiKs = []

    for x in Sp:
        Cs[x] = Symbol(x, real=True, negative=False)
        Cso[x] = Symbol(x + 'o', real=True, negative=False) * \
            (0 if conc[x] == 0 else 1)
        equivals.append((Cso[x], conc[x]))
        equiCo.append((Cso[x], conc[x]))

    KCs = []
    for i in range(len(Ks)):
        row = []
        if len(Ks[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[0], Ks[i][0]))
            equiKs.append((row[0], Ks[i][0]))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[0], Ks[i][0]))
            equiKs.append((row[0], Ks[i][0]))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[1], Ks[i][1]))
            equiKs.append((row[1], Ks[i][1]))
        KCs.append(row)

    if not molar:
        f = Matrix(propensity_vec(KCs, Cs, Rr, Rp))
    else:
        f = Matrix(propensity_vec_molar(KCs, Cs, Rr, Rp))

    if mode == "Numeric":
        f = f.subs(equivals)
    elif mode == "fofks":
        f = f.subs(equiCo)
    elif mode == "fofCo":
        f = f.subs(equiKs)

    stch_var = Matrix(V)
    # Cs might have change after call
    for x in Sp:
        Cs[x] = Symbol(x, real=True, negative=False)
    slabels = [x for x in Cs]

    Ss = []
    nz = []
    for row in range(stch_var.shape[0]):
        if sum(abs(stch_var[row, :])) != 0 and slabels[row][0] != "-":
            Ss.append(list(stch_var[row, :]))
            nz.append(row)
    Ss = Matrix(Ss)

    dA_dt = Ss * f
    # print(dA_dt)
    ccs = [Cs[x] for x in Cs]
    js = [ccs[x] for x in nz]
    A = dA_dt.jacobian(js)
    F = [[0] * len(f) for x in range(len(f))]
    for i in range(len(f)):
        F[i][i] = F[i][i] + f[i]
    F = Matrix(F)

    BBT = Ss * F * Ss.T

    cov = []
    for i in nz:
        row = []
        for j in nz:
            key = "C" + str(min(i + 1, j + 1)) + "_" + str(max(i + 1, j + 1))
            if i == j:
                row.append(Symbol(key, real=True, negative=False))
            else:
                row.append(Symbol(key, real=True))
        cov.append(row)
    cov = Matrix(cov)
    # print(Ss)
    same = {}
    for i in range(Ss.shape[0] - 1):
        for j in range(i + 1, Ss.shape[0]):
            if Ss[i, :] == Ss[j, :]:
                same[(i, j)] = 1
            elif Ss[i, :] == -Ss[j, :]:
                same[(i, j)] = -1
            else:
                pass

    L = A * cov + cov * A.T + BBT

    xs = list(dict.fromkeys(flatten(cov)))

    if items:
        text = prepare_scroll_text(items)
        def ffprint(x): return text.insert(
            INSERT, " ".join([str(y) for y in x]))
    else:
        def ffprint(x): return print(" ".join([str(y) for y in x]), end="")

    Sps = [x for x in Sp]
    CSps = {}
    for i in nz:
        for j in nz:
            if j >= i:
                key = "C" + str(min(i + 1, j + 1)) + "_" + \
                    str(max(i + 1, j + 1))
                ffprint([key, " = ", str(Sps[i]) + "_" + str(Sps[j]), "\n"])
                CSps[key] = Sps[i] + "_" + Sps[j]

    eqs = []
    ffprint(["\nEquations to solve\n\n"])
    for i in range(len(nz)):
        for j in range(i, len(nz)):
            if L[i, j] != 0:
                row = collect(L[i, j], xs)
                eqs.append(row)
                ffprint([row, "\n"])

    for en in same:
        i, j = en
        for k in range(j, Ss.shape[0]):
            eqs.append(cov[j, k] - same[en] * cov[i, k])

    post_proc = False
    multiple_cval = False

    val = solve(dA_dt, js)

    cval = {}
    if val:
        for x in js:
            key = x
            if key in val:
                cval[key] = val[key]
            else:
                post_proc = True
    else:
        post_proc = True

    if post_proc:
        ffprint(
            ["\nSteady state concentrations not properly calculated due to incomplete constraint\n\n"])
        ffprint(["\n", dA_dt, "\n"])
        ffprint(["\nAttempting the use of boundary condition", "\n"])

        tosum = set()
        sets = [tosum]
        used = set()
        for ih in range(len(Rr)):
            used2 = set()
            if len(Rr[ih]) == 1:
                for s in sets:
                    key = list(Rr[ih].keys())[0]
                    if key not in used and Rr[ih][key] != 0:
                        s.add(key)
                        used2.add(key)
            else:
                for j in range(len(sets)):
                    s = sets[j]
                    sets.append(s.copy())
                    key = list(Rr[ih].keys())[0]
                    if key not in used and Rr[ih][key] != 0:
                        s.add(key)
                        used2.add(key)
                    key = list(Rr[ih].keys())[1]
                    if key not in used and Rr[ih][key] != 0:
                        sets[-1].add(key)
                        used2.add(key)

            if len(Rp[ih]) == 1:
                for s in sets:
                    key = list(Rp[ih].keys())[0]
                    if key not in used and Rp[ih][key] != 0:
                        s.add(key)
                        used2.add(key)
            else:
                for j in range(len(sets)):
                    s = sets[j]
                    sets.append(s.copy())
                    key = list(Rp[ih].keys())[0]
                    if key not in used and Rp[ih][key] != 0:
                        s.add(key)
                        used2.add(key)
                    key = list(Rp[ih].keys())[1]
                    if key not in used and Rp[ih][key] != 0:
                        sets[-1].add(key)
                        used2.add(key)
            for zs in used2:
                used.add(zs)

        Fe = [dA_dt]
        for s in sets:
            Fe.append(sum([Cs[x] - Cso[x] for x in s]))

        if mode == "Numeric":
            Fe = [entry.subs(equivals) for entry in Fe]
        elif mode == "fofks":
            Fe = [entry.subs(equiCo) for entry in Fe]
        elif mode == "fofCo":
            Fe = [entry.subs(equiKs) for entry in Fe]

        val2 = solve(Fe, {x for x in js})
        cval = {}
        for x in js:
            key = x
            if key in val2:
                cval[key] = val2[key]
            else:
                multiple_cval = True
                break
        # print(Fe)

    ffprint(["\nUsing Algebraic manipulation of AC + CA.T + BT = 0\n"])
    if not multiple_cval:
        print(1)
        ffprint(["\nSteady state concentrations\n\n"])
        for x in cval:
            ffprint([x, " = ", cval[x], "\n\n"])
        L = simplify(subs2(L, cval))
        eqs = []
        ffprint(["\nEquations to solve\n\n"])
        for i in range(len(nz)):
            for j in range(i, len(nz)):
                if L[i, j] != 0:
                    row = collect(L[i, j], xs)
                    eqs.append(row)
                    ffprint([row, "\n\n"])

        for en in same:
            i, j = en
            for k in range(j, Ss.shape[0]):
                eqs.append(cov[j, k] - same[en] * cov[i, k])

        sol = solve(eqs, xs)

        if not sol:
            #LA = Matrix(eqs)
            #LHS = LA.jacobian(xs)
            #RHS = LA-LHS*Matrix(xs)
            #sol = simplify(-(LHS**-1)*RHS)
            sol = list(nonlinsolve(eqs, xs))
            sol = {xs[x]: sol[x] for x in range(len(xs))}
        ffprint(["\nCovariance\n\n"])
        for x in sol:
            Ch = factor(simplify(subs2(sol[x], cval)))
            #Ch = sol[x]
            ffprint(["Cov(" + str(CSps[str(x)]) + ")", " = ", Ch, "\n\n"])
    else:
        print(2)
        ffprint(["\nMultiple solutions detected"])
        valset = 0
        for val in val2:
            ffprint(["\nValues ", valset, "\n"])
            cval = {}
            for x in js:
                key = x
                if key in val:
                    cval[key] = val[key]
            ffprint(["\nSteady state concentrations\n\n"])
            for x in cval:
                ffprint([x, " = ", cval[x], "\n\n"])
            L = simplify(subs2(L, cval))
            eqs = []
            ffprint(["\nEquations to solve\n\n"])
            for i in range(len(nz)):
                for j in range(i, len(nz)):
                    if L[i, j] != 0:
                        row = collect(L[i, j], xs)
                        eqs.append(row)
                        ffprint([row, "\n\n"])

            for en in same:
                i, j = en
                for k in range(j, Ss.shape[0]):
                    eqs.append(cov[i, k] - same[en] * cov[j, k])

            sol = solve(eqs, xs)
            if not sol:
                sol = list(nonlinsolve(eqs, xs))
                sol = {xs[x]: sol[0][x] for x in range(len(xs))}
            ffprint(["\nCovariance\n\n"])
            fact = 1
            for x in sol:
                u, w = str(x).split("_")
                Ch = factor(simplify(subs2(sol[x], cval))).subs(equivals)
                if u == 'C' + w and Ch < 0:
                    fact = -1
                    break

            for x in sol:
                Ch = factor(simplify(subs2(sol[x], cval)))
                #Ch = sol[x]
                ffprint(["Cov(" + str(CSps[str(x)]) + ")",
                         " = ", Ch * fact, "\n\n"])
            valset = valset + 1

    return [0, 0]
