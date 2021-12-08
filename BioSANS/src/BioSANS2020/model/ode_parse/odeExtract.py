from sympy import sympify, expand, Matrix
import numpy as np
from BioSANS2020.gui_functs.scrollable_text import INSERT
from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text


DONE_PARSING = set()


def process(xvar):
    global DONE_PARSING
    xxvar = xvar.strip("*").strip().replace("*/", "/")
    if xxvar[0] == "/":
        xxvar = "1" + xxvar
    if xxvar[-1] == "/":
        xxvar = xxvar.strip("/")
    val = sympify(xxvar.strip("*"))
    if val in DONE_PARSING:
        return None
    else:
        DONE_PARSING.add(val)
    return val


def propExt(expr, prop):
    ex = str(expr)
    open = 0
    close = 0
    diff = 0

    collect = ""
    for v in ex:
        if v in ["+", "-"] and diff == 0:
            if len(collect) > 0:
                d = process(collect)
                if d is not None:
                    prop.append(d)
                collect = ""
                open = 0
                close = 0
        elif v == "(":
            open = open + 1
            collect = collect + v
        elif v == ")":
            close = close + 1
            collect = collect + v
        elif (v.isnumeric() or v == ".") and (diff == 0):
            if len(collect) >= 3:
                if collect[-2:] == "**" or collect[-1] not in ["*", "/"]:
                    collect = collect + v
            elif len(collect) >= 1:
                if collect[-1] not in ["*", "/", " "]:
                    collect = collect + v
            else:
                pass
        else:
            collect = collect + v
        diff = open - close

    d = process(collect)
    if d is not None:
        prop.append(d)


def termExt(expr):
    term = []
    ex = str(expr)
    open = 0
    close = 0
    diff = 0

    collect = ""
    last_sign = ""
    for v in ex:
        if v in ["+", "-"] and diff == 0:
            if len(collect) > 0:
                term.append(last_sign + collect)
                collect = ""
                open = 0
                close = 0
            last_sign = v
        elif v == "(":
            open = open + 1
            collect = collect + v
        elif v == ")":
            close = close + 1
            collect = collect + v
        else:
            collect = collect + v
        diff = open - close

    term.append(last_sign + collect)
    return term


def get_prop_stoich(dxdt):
    prop = []
    dAdt = []
    for expro in dxdt:
        expr = expand(sympify(expro))
        propExt(expr, prop)
        dAdt.append(expr)

    w_var = prop
    v_stoich = [[0 for xvar in range(len(w_var))] for y in range(len(dAdt))]

    for i in range(len(dAdt)):
        for j in range(len(w_var)):
            s = 0
            for xxvar in termExt(dAdt[i]):
                xvar = sympify(xxvar) / w_var[j]
                try:
                    s = s + float(xvar)
                except BaseException:
                    pass
            v_stoich[i][j] = s

    return Matrix(v_stoich), Matrix(w_var)


def print_stoich_prop(dxdt):
    global DONE_PARSING
    DONE_PARSING = set()
    print()
    v_stoich, w_var = get_prop_stoich(dxdt)
    for tvar in v_stoich * w_var:
        print(tvar)
    print()

    for c in np.array(v_stoich):
        print([round(y, 4) for y in c])
    print()

    for c in np.array(w_var):
        print(c)


def transform_to_rxn(xvar, dxdt, xo, ks, items):
    global DONE_PARSING

    if items:
        text = prepare_scroll_text(items)
        def ffrint(xvar): return text.insert(INSERT, xvar + "\n")
    else:
        def ffrint(xvar): return print(" ".join([str(y) for y in xvar]),
                                       end="")

    DONE_PARSING = set()
    v_stoich, w_var = get_prop_stoich(dxdt)
    S = np.around(np.array(v_stoich).astype(float), 3)
    Rxn = []
    Ksn = set()
    ind = 0

    for col in S.T:
        R = ""
        P = ""
        for i in range(len(col)):
            if col[i] != 0:
                if col[i] < 0:
                    R = R + str(abs(col[i])) + " " + xvar[i] + " " + "+ "
                else:
                    P = P + str(abs(col[i])) + " " + xvar[i] + " " + "+ "
        if R.strip() == "":
            R = "0 NONE"
        if P.strip() == "":
            P = "0 NONE"

        inSp = []
        for s in w_var[ind].free_symbols:
            ss = str(s)
            if ss in xvar:
                inSp.append(ss)
            else:
                Ksn.add(ss)
        inSp = ",".join(inSp)
        Rxn.append(R.strip("+ ") + " => " + P.strip("+ ") +
                   ", 1 ::::: lambda " + inSp + " : " + str(w_var[ind]))
        ind = ind + 1

    ffrint("Function_Definitions:")
    for k in Ksn:
        ffrint(k.strip() + " = type actual value") if k.strip(
        ) not in ks else ffrint(k.strip() + " = " + ks[k.strip()])
    for c in xvar:
        ffrint(c.strip() + "_ini = type actual value") if c.strip(
        ) not in xo else ffrint(c.strip() + "_ini = " + xo[c.strip()])

    ffrint("")
    ffrint("#REACTIONS")
    for rx in Rxn:
        ffrint(rx)

    ffrint("")
    ffrint("@CONCENTRATION")
    for c in xvar:
        ffrint(c + " , " + c.strip() + "_ini")
    ffrint("NONE, 1")
    return text


def odedxdt_to_topo(mfile, items):
    ddvar = open(mfile, "r")
    print()
    xvar = []
    xo = {}
    ks = {}
    dxdt = []
    last = ""
    for xxvar in ddvar:
        if last == "ODE_DECLARATIONS" and xxvar.strip(
        ) not in ["INI_CONCENTRATIONS:", "RATE_CONSTANTS:"]:
            if xxvar.strip() != "":
                row = xxvar.split("=")
                xvar.append(row[0].strip())
                dxdt.append(row[1])
        elif last == "INI_CONCENTRATIONS" and xxvar.strip() \
                not in ["ODE_DECLARATIONS:", "RATE_CONSTANTS:"]:
            if xxvar.strip() != "":
                row = xxvar.split("=")
                xo[row[0].strip()] = row[1].strip()
        elif last == "RATE_CONSTANTS" and xxvar.strip() \
                not in ["ODE_DECLARATIONS:", "INI_CONCENTRATIONS:"]:
            if xxvar.strip() != "":
                row = xxvar.split("=")
                ks[row[0].strip()] = row[1].strip()
        elif xxvar.strip() == "ODE_DECLARATIONS:":
            last = "ODE_DECLARATIONS"
        elif xxvar.strip() == "INI_CONCENTRATIONS:":
            last = "INI_CONCENTRATIONS"
        elif xxvar.strip() == "RATE_CONSTANTS:":
            last = "RATE_CONSTANTS"

    # print_stoich_prop(dxdt)
    return transform_to_rxn(xvar, dxdt, xo, ks, items)

# odedxdt_to_topo("NewText.txt")
