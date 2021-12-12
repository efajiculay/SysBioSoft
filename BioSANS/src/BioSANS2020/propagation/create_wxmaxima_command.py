#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.gui_functs.scrollable_text import *
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from sympy import *
from BioSANS2020.propagation.recalculate_globals import get_globals
#from BioSANS2020.myglobal import mglobals as globals2


def for_wxmaxima(Sp, ks_dict, conc, r_dict, p_dict, V, items=None, rfile=""):
    get_globals(rfile)
    if items:
        text = prepare_scroll_text(items)
        def ffprint(x): return text.insert(
            INSERT, " ".join([str(y) for y in x]))
    else:
        def ffprint(x): return print(" ".join([str(y) for y in x]), end="")

    ffprint(["/*Copy and paste to wxmaxima and run the cell*/\n\n"])
    Cs = {}
    Cso = {}
    t = Symbol('t', real=True, positive=True)
    for x in Sp:
        Cs[x] = Function(x, positive=True)(t)
        Cso[x] = Symbol(x + "o", real=True, negative=False)

    KCs = []
    for i in range(len(ks_dict)):
        row = []
        if len(ks_dict[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
        KCs.append(row)

    f = Matrix(propensity_vec_molar(KCs, Cs, r_dict, p_dict, True))
    stch_var = Matrix(V)
    # Cs might have change after call
    for x in Sp:
        Cs[x] = Function(x, positive=True)(t)
    slabels = [x for x in Cs]
    Ss = []
    nz = []
    for row in range(stch_var.shape[0]):
        if sum(abs(stch_var[row, :])) != 0 and slabels[row][0] != "-":
            Ss.append(list(stch_var[row, :]))
            nz.append(row)
    Ss = Matrix(Ss)

    dA_dt = Ss * f
    ccs = [Cs[x] for x in Cs]
    ccso = [Cso[x] for x in Cso]
    js = [ccs[x] for x in nz]
    jso = [ccso[x] for x in nz]

    dA_dt2 = []
    for x in js:
        dA_dt2.append(x.diff(t))

    Fe = []
    for i in range(len(dA_dt2)):
        ffprint(["f" + str(i + 1) + ":" + str(dA_dt2[i]).replace("Derivative",
                                                                 "diff") + " = " + str(dA_dt[i]) + ";", "\n"])
        Fe.append("f" + str(i + 1))
    ffprint(["\n"])
    for x in range(len(js)):
        ffprint(["atvalue(" + str(js[x]) + ",t=0," + str(jso[x]) + ");", "\n"])

    ffprint(["\ndesolve(" + str(Fe).replace('"', '').replace("'", '') +
             "," + str([x for x in js]) + ");"])
    ffprint(["\n"])
    ffprint(["\n\n/*Copy and paste to wxmaxima and run the cell*/\n"])
    ffprint(["\n"])

    for x in Sp:
        Cso[x] = conc[x]

    f = Matrix(propensity_vec_molar(ks_dict, Cs, r_dict, p_dict, True))
    stch_var = Matrix(V)
    # Cs might have change after call
    for x in Sp:
        Cs[x] = Function(x, positive=True)(t)
    slabels = [x for x in Cs]
    Ss = []
    nz = []
    for row in range(stch_var.shape[0]):
        if sum(abs(stch_var[row, :])) != 0 and slabels[row][0] != "-":
            Ss.append(list(stch_var[row, :]))
            nz.append(row)
    Ss = Matrix(Ss)

    dA_dt = Ss * f
    ccs = [Cs[x] for x in Cs]
    ccso = [Cso[x] for x in Cso]
    js = [ccs[x] for x in nz]
    jso = [ccso[x] for x in nz]

    dA_dt2 = []
    for x in js:
        dA_dt2.append(x.diff(t))

    Fe = []
    for i in range(len(dA_dt2)):
        ffprint(["f" + str(i + 1) + ":" + str(dA_dt2[i]).replace("Derivative",
                                                                 "diff") + " = " + str(dA_dt[i]) + ";", "\n"])
        Fe.append("f" + str(i + 1))
    ffprint(["\n"])
    for x in range(len(js)):
        ffprint(["atvalue(" + str(js[x]) + ",t=0," + str(jso[x]) + ");", "\n"])

    ffprint(["\ndesolve(" + str(Fe).replace('"', '').replace("'", '') +
             "," + str([x for x in js]) + ");"])

    return [0, 0]
