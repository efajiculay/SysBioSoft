#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from scipy.integrate import odeint
from BioSANS2020.propagation.deterministic.ode_int import *
from BioSANS2020.model.param_est.my_mcem import *
from tkinter import filedialog
from BioSANS2020.gui_functs.scrollable_text import *
from scipy import optimize
from BioSANS2020.myglobal import proc_global as proc_global


def load_data(file=None):
    if file is None:
        file = filedialog.askopenfilename(title="Select file")
    with open(file, "r") as f:
        data = []
        ddvar = []
        row1 = str(f.readline()).strip()
        slabels = row1.split("\t")[1:]
        for row in f:
            cols = [float(x) for x in row.split("\t")]
            if cols[0] == 0.0 and len(ddvar) > 0:
                data.append(np.array(ddvar))
                ddvar = []
            ddvar.append(cols)
        data.append(np.array(ddvar))
        return (data, slabels)


def custom_likelihood(ks, args=None):
    data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss, molar, rfile = args
    ind = 0
    slabels, data = data
    for row in Kmiss:
        i, j = row
        Ks[i][j] = ks[ind]
        ind = ind + 1
    for s in Cmiss:
        conc[s] = ks[ind]
        ind = ind + 1
    # print(conc)
    z = ODE_int(conc, t, Sp, Ks, Rr, Rp, V, molar, rfile)
    spc1 = np.array([s for s in Sp])
    # print(z[0,np.isin(spc1,Cmiss)])
    spc1 = ~np.isin(spc1, Cmiss)
    spc2 = np.array(slabels)
    spc2 = ~np.isin(spc2, Cmiss)
    return -np.sum((data[:, spc2] - z[:, spc1])**2)


def ave_abs_dev(ks, args=None):
    data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss, molar, rfile = args
    ind = 0
    slabels, data = data
    for row in Kmiss:
        i, j = row
        Ks[i][j] = ks[ind]
        ind = ind + 1
    for s in Cmiss:
        conc[s] = ks[ind]
        ind = ind + 1
    z = ODE_int(conc, t, Sp, Ks, Rr, Rp, V, molar, rfile)
    spc1 = np.array([s for s in Sp])
    spc1 = ~np.isin(spc1, Cmiss)
    spc2 = np.array(slabels)
    spc2 = ~np.isin(spc2, Cmiss)
    denom = data[:, spc2]
    denom[np.abs(denom) < 1.0e-10] = 1
    return np.mean(np.abs((data[:, spc2] - z[:, spc1]) / denom))


def labelParam(Kmiss, Cmiss, ks):
    pars = {}
    count = 0
    L = ["kf", "kb"]
    for x in Kmiss:
        i, j = x
        pars[L[j] + str(i + 1)] = ks[count]
        count = count + 1

    for x in Cmiss:
        pars[x + "o"] = ks[count]
        count = count + 1
    return pars


def param_estimate(conc, t, Sp, Ks, Rr, Rp, V, items,
                   molar=False, mode="MCEM", TrueDataFil=None, rfile=None):

    data2 = load_data(TrueDataFil)
    slabels = data2[1]
    Spc = [s for s in Sp]
    z = [conc[a] for a in Sp]
    Cmiss = []
    for i in range(len(z)):
        if z[i] < 0:
            Cmiss.append(Spc[i])
    Kmiss = []
    for i in range(len(Ks)):
        if len(Ks[i]) == 1:
            if Ks[i][0] < 0:
                Kmiss.append((i, 0))
        elif len(Ks[i]) == 2:
            if Ks[i][0] < 0:
                Kmiss.append((i, 0))
            if Ks[i][1] < 0:
                Kmiss.append((i, 1))
    npar = len(Cmiss) + len(Kmiss)

    t = data2[0][0][:, 0]
    data = (slabels, data2[0][0][:, 1:])
    # print(data)

    if items:
        text = prepare_scroll_text(items)
        def ffprint(x): return text.insert(
            INSERT, " ".join([str(y) for y in x]))
    else:
        def ffprint(x): return print(" ".join([str(y) for y in x]))

    param_res = {}
    if mode == "MCEM":
        f = 1
        k = 10
        ks, er_min = run_mcem(f, npar, f, k * f, positive_only=True, likelihood=custom_likelihood,
                              arg=(data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss, molar, rfile))
        er = ave_abs_dev(ks, (data, conc, t, Sp, Ks, Rr,
                              Rp, V, Cmiss, Kmiss, molar, rfile))
        count = 0
        best_er = 1e+100
        best_ks = ks

        rands = [(x + 1) * np.random.uniform(0, 1) for x in range(200000)]
        ind = 0
        while er > 1.0e-6 and k <= 10000:
            ks, er_min = run_mcem(min(f, 10), npar, f, k * f, positive_only=True, likelihood=custom_likelihood,
                                  arg=(data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss, molar, rfile), rr=rands[ind])
            er = ave_abs_dev(ks, (data, conc, t, Sp, Ks, Rr,
                                  Rp, V, Cmiss, Kmiss, molar, rfile))
            if er < best_er:
                best_er = er
                best_ks = ks

                param_res = labelParam(Kmiss, Cmiss, best_ks)
                ffprint(["\n\npartial result =\n", "errors = ", er, er_min])
                ffprint(["\n"])
                for x in param_res:
                    ffprint(["\n" + str(x) + " = " + str(param_res[x])])

            f = f + 1
            count = count + 1
            if count == 10:
                k = k * 10
                f = 1
                count = 0
            ind = ind + 1
        param_res = labelParam(Kmiss, Cmiss, best_ks)
        ffprint(["\n\nFinal result =\n", best_ks, er, er_min])
        ffprint(["\n"])
        for x in param_res:
            ffprint(["\n" + str(x) + " = " + str(param_res[x])])
    elif mode == "DEvol":
        xo = []
        for row in Kmiss:
            xo.append(np.random.uniform())
        for s in Cmiss:
            xo.append(np.random.uniform())
        ddvar = optimize.minimize(lambda ks: -custom_likelihood(ks, (data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss,
                                                                     molar, rfile)), xo, method='Nelder-Mead', tol=1e-10, options={'maxiter': 10000, 'adaptive': True})

        bounds = []
        counts = 0
        for row in Kmiss:
            val = ddvar.x[counts]
            bounds.append((max(0, 0.1 * val), max(0, 10 * val)))
            counts = counts + 1
        for s in Cmiss:
            val = ddvar.x[counts]
            bounds.append((max(0, 0.1 * val), max(0, 10 * val)))
            counts = counts + 1

        ddvar = optimize.differential_evolution(lambda ks: -custom_likelihood(
            ks, (data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss, molar, rfile)), bounds, maxiter=100, popsize=10, tol=1.0e-10)
        ddvar = optimize.minimize(lambda ks: -custom_likelihood(ks, (data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss,
                                                                     molar, rfile)), ddvar.x, method='Nelder-Mead', tol=1e-10, options={'maxiter': 10000, 'adaptive': True})
        param_res = labelParam(Kmiss, Cmiss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for x in param_res:
            ffprint(["\n" + str(x) + " = " + str(param_res[x])])
    elif mode == "NeldMead":
        xo = []
        for row in Kmiss:
            xo.append(np.random.uniform())
        for s in Cmiss:
            xo.append(np.random.uniform())
        ddvar = optimize.minimize(lambda ks: -custom_likelihood(ks, (data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss,
                                                                     molar, rfile)), xo, method='Nelder-Mead', tol=1e-10, options={'maxiter': 100000, 'adaptive': True})
        param_res = labelParam(Kmiss, Cmiss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for x in param_res:
            ffprint(["\n" + str(x) + " = " + str(param_res[x])])
    elif mode == "Powell":
        xo = []
        for row in Kmiss:
            xo.append(np.random.uniform())
        for s in Cmiss:
            xo.append(np.random.uniform())
        ddvar = optimize.minimize(lambda ks: -custom_likelihood(ks, (data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss,
                                                                     molar, rfile)), xo, method='Powell', tol=1e-10, options={'maxiter': 100000, 'ftol': 1e-10, 'xtol': 1e-10})
        param_res = labelParam(Kmiss, Cmiss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for x in param_res:
            ffprint(["\n" + str(x) + " = " + str(param_res[x])])
    elif mode == "L-BFGS-B":
        xo = []
        for row in Kmiss:
            xo.append(np.random.uniform())
        for s in Cmiss:
            xo.append(np.random.uniform())
        ddvar = optimize.minimize(lambda ks: -custom_likelihood(ks, (data, conc, t, Sp, Ks, Rr, Rp, V, Cmiss, Kmiss, molar, rfile)),
                                  xo, method='L-BFGS-B', tol=1e-10, options={'maxiter': 100000, 'ftol': 1e-10, 'xtol': 1e-10})
        param_res = labelParam(Kmiss, Cmiss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for x in param_res:
            ffprint(["\n" + str(x) + " = " + str(param_res[x])])
    return (0, param_res)
