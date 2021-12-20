"""

                This module is the analytical_sol module

This module handles analytical expression derivation

The following are the functions in this module

1. solve_with_timeout
2. get_sets
3. grab_steady_state
4. cs_to_csr
5. analyt_soln

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from func_timeout import func_timeout, FunctionTimedOut
from sympy import simplify, Matrix, solve, Add, integrate, Function, \
    dsolve, Symbol, exp, Eq, zeros
# import time

from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.propagation.propensity import propensity_vec_molar
# from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.propagation.recalculate_globals import get_globals


def solve_with_timeout(f_e, j_s):
    """resturn solution to the system of equation in f_e with variables
    symbos in the list j_s

    Args:
        f_e (list): equations
        j_s (list): Symbols

    Returns:
        dict: sometimes list but normally {variable : value}
    """
    sol = None
    sol = solve(f_e, j_s)
    return sol


def get_sets(r_dict, p_dict):
    """[summary]

    Args:
        r_dict (dict): dictionary of reactants
        p_dict (dict): dictionary of products

    Returns:
        set: sets to sum that will be used as constraints
    """
    tosum = set()
    sets = [tosum]
    used = set()
    for ih_var, _ in enumerate(r_dict):
        used2 = set()
        if len(r_dict[ih_var]) == 1:
            for svar in sets:
                key = list(r_dict[ih_var].keys())[0]
                if key not in used and r_dict[ih_var][key] != 0:
                    svar.add(key)
                    used2.add(key)
        else:
            for j, _ in enumerate(sets):
                svar = sets[j]
                sets.append(svar.copy())
                key = list(r_dict[ih_var].keys())[0]
                if key not in used and r_dict[ih_var][key] != 0:
                    svar.add(key)
                    used2.add(key)
                key = list(r_dict[ih_var].keys())[1]
                if key not in used and r_dict[ih_var][key] != 0:
                    sets[-1].add(key)
                    used2.add(key)

        if len(p_dict[ih_var]) == 1:
            for svar in sets:
                key = list(p_dict[ih_var].keys())[0]
                if key not in used and p_dict[ih_var][key] != 0:
                    svar.add(key)
                    used2.add(key)
        else:
            for j, _ in enumerate(sets):
                svar = sets[j]
                sets.append(svar.copy())
                key = list(p_dict[ih_var].keys())[0]
                if key not in used and p_dict[ih_var][key] != 0:
                    svar.add(key)
                    used2.add(key)
                key = list(p_dict[ih_var].keys())[1]
                if key not in used and p_dict[ih_var][key] != 0:
                    sets[-1].add(key)
                    used2.add(key)
        for z_s in used2:
            used.add(z_s)
    return sets


def grab_steady_state(r_dict, p_dict, c_s, c_so, da_dt, j_s):
    """This function return steady state amounts of components

    Args:
        r_dict (dict): reactants dictionary
        p_dict (dict): products dictionary
        c_s (dict): dictionary of concentration Symbols/values
        c_so (dict): dictionary of initial concentration Symbols/values
        da_dt ([type]): ode expression
        j_s ([type]): list of components

    Returns:
        Matrix: steady state result
    """
    sets = get_sets(r_dict, p_dict)
    fez = [da_dt]
    for svar in sets:
        fez.append(sum([c_s[xvar] - c_so[xvar] for xvar in svar]))
    val2 = solve(fez, list(j_s))
    x_s = Matrix([val2[j_s[ih_var]] for ih_var in range(len(j_s))])
    return x_s


def cs_to_csr(f_x, c_s, cs_r, not_semi):
    """[summary]

    Args:
        f_x (expression): sympy expression
        c_s (Symbol): sympy symbols
        cs_r (Symbol): equivalent of c_s
        not_semi (boolean): if True, use pure symbolic

    Returns:
        expression: substituted expression
    """
    if not_semi:
        r_r = []
        for xvar in c_s:
            r_r.append((c_s[xvar], cs_r[xvar]))
        return f_x.subs(r_r)
    return f_x


def analyt_soln(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var, items=None,
                rfile="", not_semi=True, mode=None):
    """This function finds an analytical solution to a system of ODE.

    Args:
        sp_comp (dict): dictionary of components
        ks_dict (dict):  dictionary of rate constant
        conc (dict): dictionary of concentrations
        r_dict (dict): dictionary of reactants
        p_dict (dict): dictionary of products
        stch_var (Matrix): stoichiomentric matrix
        items (tuple, optional): (canvas, scroll_x, scrolly).
        rfile (str, optional): BioSANS topology file.
        not_semi (bool, optional): if True, use pure symbolic solution;
                which is a function of time, initial concentration, and
                rate constant. Defaults to True.
        mode (str, optional): mode of computation;
                ftxo - function of time and initial condition
                ftks - function of time and rate constant
            Defaults to None.

    Returns:
        dict: dictionary of solution
    """
    get_globals(rfile)
    my_equations = {}

    if items:
        text = prepare_scroll_text(items)

        def ffprint(xvar): return text.insert(
            'insert', " ".join([str(yvar) for yvar in xvar]))
    elif items == 0:
        def ffprint(xvar):
            return print("", end="")
    else:
        def ffprint(xvar):
            return print(" ".join([str(yvar) for yvar in xvar]), end="")

    if not_semi:
        ffprint(["\n\nComplex Analytical expressions\n"])
        ffprint(
            ["\n\nThe complex expression is because sympy do not know how",
             "you want to simplify the expression\n\n"])
    else:
        ffprint(["Simple semi-analytical expression\n\n"])

    c_s = {}
    c_so = {}
    tvar = Symbol('t', negative=False, real=True)

    equivals = [(tvar, 100)]
    cs_r = {}
    for xvar in sp_comp:
        c_s[xvar] = Function(xvar, negative=False, real=True)
        c_s[xvar] = c_s[xvar](tvar)
        cs_r[xvar] = Symbol(xvar, negative=False, real=True)
        if mode == "ftxo":
            c_so[xvar] = Symbol(xvar + "o", negative=False, real=True)
        elif mode == "ftks":
            c_so[xvar] = conc[xvar]
        else:
            c_so[xvar] = Symbol(xvar + "o", negative=False,
                                real=True) if not_semi else conc[xvar]
        equivals.append((c_so[xvar], conc[xvar]))

    kcs = []
    for i, _ in enumerate(ks_dict):
        row = []
        if len(ks_dict[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, negative=False, real=True))
            equivals.append((row[0], ks_dict[i][0]))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, negative=False, real=True))
            equivals.append((row[0], ks_dict[i][0]))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, negative=False, real=True))
            equivals.append((row[1], ks_dict[i][1]))
        kcs.append(row)

    if mode == "ftks":
        fvar = Matrix(propensity_vec_molar(kcs, c_s, r_dict, p_dict, True))
    elif mode == "ftxo":
        fvar = Matrix(propensity_vec_molar(ks_dict, c_s, r_dict, p_dict, True))
    else:
        fvar = Matrix(propensity_vec_molar(kcs, c_s, r_dict, p_dict, True)) \
            if not_semi else Matrix(
                propensity_vec_molar(ks_dict, c_s, r_dict, p_dict, True))

    stch_var = Matrix(stch_var)
    # c_s might have change after call
    for xvar in sp_comp:
        c_s[xvar] = Function(xvar, negative=False, real=True)
        c_s[xvar] = c_s[xvar](tvar)
    slabels = list(c_s.keys())  # [xvar for xvar in c_s]

    def stoich(xvar):
        return max(1, abs(stch_var[slabels.index(xvar)][0]))

    ss_var = []
    nz_var = []
    for row in range(stch_var.shape[0]):
        if sum(abs(stch_var[row, :])) != 0 and slabels[row][0] != "-":
            ss_var.append(list(stch_var[row, :]))
            nz_var.append(row)
    ss_var = Matrix(ss_var)

    da_dt = ss_var * fvar
    ccs = [c_s[xvar] for xvar in c_s]
    ccso = [c_so[xvar] for xvar in c_so]
    j_s = [ccs[xvar] for xvar in nz_var]
    jso = [ccso[xvar] for xvar in nz_var]

    da_dt2 = []
    for xvar in j_s:
        da_dt2.append(xvar.diff(tvar))

    fes = [Eq(da_dt2[i], da_dt[i]) for i in range(len(da_dt2))]

    left_hand = da_dt.jacobian(j_s)
    # print(left_hand)
    right_hand = da_dt - left_hand * Matrix(j_s)
    # print(right_hand)
    # print(fes)
    # print(da_dt)
    # print(j_s)

    check = True
    work = False
    for i in range(left_hand.shape[0]):
        if sum(abs(left_hand.col(i))) == 0:
            check = False

    if check:
        print("Inside Check")
        left_hand_j = np.array(left_hand)
        # print(left_hand_j)
        left_hand_j = Matrix(left_hand_j.flatten()).jacobian(j_s)
        # print(left_hand_j)
        if left_hand_j != zeros(left_hand_j.shape[0], left_hand_j.shape[1]):
            sets = get_sets(r_dict, p_dict)
            fe2 = []
            for svar in sets:
                fe2.append(
                    Eq(sum([c_s[xvar] / stoich(xvar) - c_so[xvar]
                            / stoich(xvar) for xvar in svar]), 0))
            xret = solve(fe2, j_s)
            # print(xret)
            # print(fes)

            for svar in xret:
                # print(fes[-1].atoms(Function))
                fes[-1] = Eq(fes[-1].lhs, fes[-1].rhs.subs(svar, xret[svar]))
            denom = cs_to_csr(fes[-1].rhs, c_s, cs_r, not_semi)
            itvar = cs_to_csr(j_s[-1], c_s, cs_r, not_semi)
            # pot_an = integrate(1/denom,(itvar,jso[-1],j_s[-1]))-tvar
            pot_an = integrate(
                1 / fes[-1].rhs, (j_s[-1], jso[-1], j_s[-1])) - tvar
            fe2.append(pot_an)
            try:
                sol = func_timeout(20, solve, args=(fe2, j_s))
                for c_c in sol:
                    rval = Matrix(c_c).subs(equivals)
                    if np.sum(np.array(rval) < 0) > 0:
                        pass
                    else:
                        x_sol = [c_c[xvar] for xvar in range(len(j_s))]
                        ind = 0
                        for xvar in x_sol:
                            var_sp = str(j_s[ind])
                            answer = str(simplify(xvar))
                            ffprint([var_sp, " = ", answer, "\n\n"])
                            ind = ind + 1
                            if var_sp not in my_equations:
                                my_equations[var_sp] = answer
                        break
                work = True
            except FunctionTimedOut:
                work = False
        else:
            print("Attempting Matrix method")
            svar = Symbol('s', negative=False, real=True)
            x_o = Matrix(jso)
            if sum(abs(right_hand)) == 0:
                # print("llll",1)
                try:
                    x_s = grab_steady_state(
                        r_dict, p_dict, c_s, c_so, da_dt, j_s)
                    try:
                        exp_d = func_timeout(
                            200, exp, args=(left_hand * tvar,))
                        x_sol = x_s + exp_d * (x_o - x_s)
                        ind = 0
                        for xvar in x_sol:
                            var_sp = str(j_s[ind])
                            answer = str(simplify(xvar))
                            ffprint([var_sp, " = ", answer, "\n\n"])
                            ind = ind + 1
                            if var_sp not in my_equations:
                                my_equations[var_sp] = answer
                        work = True
                    except FunctionTimedOut:
                        work = False
                except:
                    pass
            else:
                # print("llll",2)
                try:
                    try:
                        x_s = simplify(-left_hand**-1 * right_hand)
                        x_sol = x_s + exp(left_hand * tvar) * (x_o - x_s)
                    except:
                        x_sol = simplify(
                            exp(left_hand * tvar) * x_o
                            + integrate(exp(left_hand * tvar)
                                        * right_hand, (svar, 0, tvar)))
                    ind = 0
                    for xvar in x_sol:
                        var_sp = str(j_s[ind])
                        answer = str(simplify(xvar))
                        ffprint([var_sp, " = ", answer, "\n\n"])
                        ind = ind + 1
                        if var_sp not in my_equations:
                            my_equations[var_sp] = answer
                    work = True
                except:
                    pass

    if not work:
        used3 = {}
        try:
            print(1)
            try:
                # sol = dsolve(fes, j_s, ics={j_s[j].subs(tvar,0):jso[j]
                #                             for j in range(len(j_s))})
                sol = func_timeout(
                    200, lambda: dsolve(
                        fes, j_s, ics={j_s[j].subs(tvar, 0): jso[j]
                                       for j in range(len(j_s))}))
                for xvar, _ in enumerate(sol):
                    var_sp = str(sol[xvar].lhs)
                    answer = str(simplify(sol[xvar].rhs))
                    ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                    if var_sp not in my_equations:
                        my_equations[var_sp] = answer
                    used3[var_sp] = True

            except FunctionTimedOut:
                float("jjj")
        except:
            try:
                print(2)
                ind = 0
                not_f = []
                good = []
                for fex in fes:
                    ics = {j_s[ind].subs(tvar, 0): jso[ind]}
                    spc = list(fex.lhs.atoms(Function))
                    spr = list(fex.rhs.atoms(Function))
                    if len(spr) == 1 and spc[0] == spr[0]:
                        # print(3)
                        sol = [dsolve(fex, ics=ics)]
                        good.append([sol[0].lhs, sol[0].rhs])
                        if str(sol[0].lhs) not in used3:
                            var_sp = str(sol[0].lhs)
                            answer = str(simplify(sol[0].rhs))
                            ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                            if var_sp not in my_equations:
                                my_equations[var_sp] = answer
                            used3[var_sp] = True
                    else:
                        # print(4)
                        # print([fex,ics])
                        not_f.append([fex, ics])
                    ind = ind + 1

                # print(good)
                # print(not_f)
                ffprint(["\n\n", "Attempting to evaluate", "\n\n"])
                flen = 2 * len(not_f)
                c_c = 0

                sets = get_sets(r_dict, p_dict)
                fez = []
                # print(slabels,stch_var)
                for svar in sets:
                    fez.append(
                        Eq(sum([c_s[xvar] / stoich(xvar) - c_so[xvar]
                                / stoich(xvar) for xvar in svar]), 0))
                xret = solve(fez, j_s)
                # print(xret)

                while len(not_f) > 0 and c_c < flen:
                    xvar = not_f.pop(0)
                    # print(xvar,100000000000000000000000000000000000)
                    spc = list(xvar[0].lhs.atoms(Function))
                    spr = list(xvar[0].rhs.atoms(Function))
                    for svar in xret:
                        # print(svar not in spc)
                        if svar not in spc:
                            subx = Eq(Symbol("Something_not_needed"),
                                      xret[svar]).rhs.atoms(Function)
                            if spc[0] in subx:
                                # print(xvar[0],xvar[1],111)
                                xvar[0] = xvar[0].subs(svar, xret[svar])
                                # print(xvar[0],xvar[1],222)
                    # print(xvar,100000000000000000000000000000000000)
                    spr = list(xvar[0].rhs.atoms(Function))
                    terms = len(Add.make_args(xvar[0].rhs))
                    if len(spr) == 1 and spc[0] == spr[0]:
                        try:
                            # print(xvar[0],xvar[1],333)
                            sol = [dsolve(xvar[0], ics=xvar[1])]
                            fez.append(Eq(sol[0].lhs[0] - sol[0].rhs, 0))
                            good.append([sol[0].lhs, sol[0].rhs])
                            if str(sol[0].lhs) not in used3:
                                var_sp = str(sol[0].lhs)
                                answer = str(simplify(sol[0].rhs))
                                ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                                if var_sp not in my_equations:
                                    my_equations[var_sp] = answer
                                used3[var_sp] = True
                        except:
                            denom = cs_to_csr(xvar[0].rhs, c_s, cs_r, not_semi)
                            itvar = cs_to_csr(spc[0], c_s, cs_r, not_semi)
                            val = integrate(
                                1 / denom, (itvar,
                                            xvar[1][spc[0].subs(tvar, 0)],
                                            spc[0])) - tvar
                            fez.append(val)
                            sol = solve(val, spc[0])
                            if isinstance(sol, dict):
                                sol = simplify(sol[spc[0]])
                            else:
                                sol = simplify(sol[0])
                            good.append([spc[0], sol])
                            var_sp = str(spc[0])
                            answer = str(sol)
                            ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                            if var_sp not in my_equations:
                                my_equations[var_sp] = answer
                            used3[var_sp] = True
                    elif len(spr) == 2 and terms == 1 \
                        and (len(r_dict) == 1 or len(p_dict) == 1) \
                            and not_semi:
                        subt = Eq(
                            spr[0] - spr[1], jso[j_s.index(spr[0])]
                            - jso[j_s.index(spr[1])])
                        if spr[0] != spc[0]:
                            vari = spr[0]
                        else:
                            vari = spr[1]
                        sol = solve(subt, vari)
                        if isinstance(sol, dict):
                            sol = sol[vari]
                        else:
                            sol = sol[0]
                        subx = Eq(Symbol("Something_not_needed"),
                                  sol).rhs.atoms(Function)
                        if spc[0] in subx:
                            # print(xvar[0],xvar[1],444)
                            xvar[0] = xvar[0].subs(vari, sol)
                            # print(xvar[0],xvar[1],555)
                        not_f.append([xvar[0], xvar[1]])
                    else:
                        not_f.append([xvar[0], xvar[1]])
                    c_c = c_c + 1

                if len(not_f) > 0 and len(used3) < len(j_s):
                    ffprint(
                        ["\n", "The following are difficult to solve",
                         "using dsolve", "\n\n"])
                    for xvar in not_f:
                        ffprint(["\n", xvar[0], "\n\n"])

                    ffprint(
                        ["\n", "The following expression may be wrong,",
                         "please compare with ode_int result :", "\n\n"])
                    c_c = 0
                    while len(not_f) > 0 and c_c < flen:
                        xvar = not_f.pop(0)
                        # print(xvar[0])
                        for svar in xret:
                            if svar not in xvar[0].lhs.atoms(Function):
                                xvar[0] = xvar[0].subs(svar, xret[svar])
                        # print(xvar[0])

                        # for yvar in good:
                            # if yvar[0] not in xvar[0].lhs.atoms(Function):
                                # xvar[0] = xvar[0].subs(yvar[0],yvar[1])
                        try:
                            spc = list(xvar[0].lhs.atoms(Function))
                            spr = list(xvar[0].rhs.atoms(Function))
                            if len(spr) == 1 and spc[0] == spr[0]:
                                spc = spc[0]
                                try:
                                    denom = cs_to_csr(
                                        xvar[0].rhs, c_s, cs_r, not_semi)
                                    itvar = cs_to_csr(spc, c_s, cs_r, not_semi)
                                    val = integrate(
                                        1 / denom, (itvar,
                                                    xvar[1][spc.subs(tvar, 0)],
                                                    spc)) - tvar
                                except:
                                    denom = cs_to_csr(
                                        xvar[0].rhs, c_s, cs_r, False)
                                    itvar = cs_to_csr(spc, c_s, cs_r, False)
                                    val = integrate(
                                        1 / denom, (itvar,
                                                    xvar[1][spc.subs(tvar, 0)],
                                                    spc)) - tvar

                                ffprint(
                                    ["\n", "Not simplified solution", "\n\n"])
                                ffprint(["need to solve for ", spc,
                                         " to simplify :\n", val, "\n\n"])
                                fez.append(val)
                                sol = solve(val, spc)
                                if isinstance(sol, dict):
                                    sol = sol[spc]
                                else:
                                    sol = simplify(sol[0])
                                good.append([spc, sol])
                                # var_sp = str(spc)
                                # answer = str(sol)
                                # ffprint(["\n",var_sp," = ",answer,"\n\n"])
                                # if var_sp not in my_equations:
                                # my_equations[var_sp] = answer
                            else:
                                not_f.append([xvar[0], xvar[1]])
                        except:
                            not_f.append([xvar[0], xvar[1]])
                        c_c = c_c + 1

                    ffprint(["\n", "Attempting to simplify", "\n\n"])
                    # print(fez)
                    sol = solve(fez, j_s)
                    for c_c in sol:
                        rval = Matrix(c_c).subs(equivals)
                        if np.sum(np.array(rval) < 0) > 0:
                            pass
                        else:
                            x_sol = [c_c[xvar] for xvar in range(len(j_s))]
                            ind = 0
                            for xvar in x_sol:
                                var_sp = str(j_s[ind])
                                answer = str(simplify(xvar))
                                ffprint([var_sp, " = ", answer, "\n\n"])
                                if var_sp not in my_equations:
                                    my_equations[var_sp] = answer
                                ind = ind + 1
                            break
            except:
                print(3)
                ind = 0
                not_f = []
                good = []
                for fex in fes:
                    ics = {j_s[ind].subs(tvar, 0): jso[ind]}
                    spc = list(fex.lhs.atoms(Function))
                    spr = list(fex.rhs.atoms(Function))
                    if len(spr) == 1 and spc[0] == spr[0]:
                        # print(3)
                        sol = [dsolve(fex, ics=ics)]
                        good.append([sol[0].lhs, sol[0].rhs])
                        if str(sol[0].lhs) not in used3:
                            var_sp = str(sol[0].lhs)
                            answer = str(sol[0].rhs)
                            ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                            if var_sp not in my_equations:
                                my_equations[var_sp] = answer
                            used3[var_sp] = True
                    else:
                        # print(4)
                        # print([fex,ics])
                        not_f.append([fex, ics])
                    ind = ind + 1

                # print(not_f,1111111111111111111)
                ffprint(["\n\n", "Attempting to evaluate", "\n\n"])
                flen = 2 * len(not_f)
                c_c = 0

                sets = get_sets(r_dict, p_dict)
                fez = []
                # print(slabels,stch_var)

                while len(not_f) > 0 and c_c < flen:
                    xvar = not_f.pop(0)
                    # print(xvar,100000000000000000000000000000000000)
                    spc = list(xvar[0].lhs.atoms(Function))
                    spr = list(xvar[0].rhs.atoms(Function))
                    terms = len(Add.make_args(xvar[0].rhs))
                    right_ok = True
                    for s_s in spr:
                        if s_s in j_s:
                            right_ok = False
                            break
                    # print(xvar,100000000000000000000000000000000000)
                    if (len(spr) == 1 and spc[0] == spr[0]) or right_ok:
                        try:
                            # print(xvar[0],xvar[1],333)
                            sol = [dsolve(xvar[0], ics=xvar[1])]
                            fez.append(Eq(sol[0].lhs[0] - sol[0].rhs, 0))
                            good.append([sol[0].lhs, sol[0].rhs])
                            if str(sol[0].lhs) not in used3:
                                var_sp = str(sol[0].lhs)
                                answer = str(sol[0].rhs)
                                ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                                if var_sp not in my_equations:
                                    my_equations[var_sp] = answer
                                used3[var_sp] = True
                        except:
                            denom = cs_to_csr(xvar[0].rhs, c_s, cs_r, not_semi)
                            itvar = cs_to_csr(spc[0], c_s, cs_r, not_semi)
                            val = integrate(
                                1 / denom, (itvar,
                                            xvar[1][spc[0].subs(tvar, 0)],
                                            spc[0])) - tvar
                            fez.append(val)
                            sol = solve(val, spc[0])
                            if isinstance(sol, dict):
                                sol = simplify(sol[spc[0]])
                            else:
                                sol = simplify(sol[0])
                            good.append([spc[0], sol])
                            var_sp = str(spc[0])
                            answer = str(sol)
                            ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                            if var_sp not in my_equations:
                                my_equations[var_sp] = answer
                            used3[var_sp] = True
                    else:
                        not_f.append([xvar[0], xvar[1]])
                    c_c = c_c + 1

                # print(not_f)
                if len(not_f) > 0 and len(used3) < len(j_s):
                    ffprint(
                        ["\n", "The following are difficult to solve",
                         "using dsolve", "\n\n"])
                    for xvar in not_f:
                        ffprint(["\n", xvar[0], "\n\n"])

                    ffprint(
                        ["\n", "The following expression may be wrong",
                         "please compare with ode_int result :", "\n\n"])
                    c_c = 0
                    while len(not_f) > 0 and c_c < flen:
                        xvar = not_f.pop(0)
                        # print(xvar[0],1111)
                        for yvar in good:
                            if yvar[0] not in xvar[0].lhs.atoms(Function):
                                xvar[0] = xvar[0].subs(yvar[0], yvar[1])
                        # print(xvar[0],2222)
                        try:
                            spc = list(xvar[0].lhs.atoms(Function))
                            spr = list(xvar[0].rhs.atoms(Function))
                            right_ok = True
                            for s_s in spr:
                                if s_s in j_s:
                                    right_ok = False
                                    break

                            # print((len(spr) == 1 and spc[0] == spr[0]) or \
                            # right_ok, (len(spr) == 1 and spc[0] == spr[0]), \
                            # right_ok)
                            if (len(spr) == 1 and spc[0]
                                    == spr[0]) or right_ok:
                                # print(111111111)
                                spc = spc[0]
                                try:
                                    denom = cs_to_csr(
                                        xvar[0].rhs, c_s, cs_r, not_semi)
                                    itvar = cs_to_csr(spc, c_s, cs_r, not_semi)
                                    # print(
                                    #     1/denom,
                                    #     (itvar,xvar[1][spc.subs(tvar,0)],
                                    #      spc),1111)
                                    val = integrate(
                                        1 / denom,
                                        (itvar, xvar[1][spc.subs(tvar, 0)],
                                         spc)) - tvar
                                except:
                                    denom = cs_to_csr(
                                        xvar[0].rhs, c_s, cs_r, False)
                                    itvar = cs_to_csr(spc, c_s, cs_r, False)
                                    # print(1/denom,
                                    #       (itvar,xvar[1][spc.subs(tvar,0)],
                                    #        spc),2222)
                                    val = integrate(
                                        1 / denom,
                                        (itvar, xvar[1][spc.subs(tvar, 0)],
                                         spc)) - tvar

                                ffprint(
                                    ["\n", "Not simplified solution", "\n\n"])
                                ffprint(["need to solve for ", spc,
                                         " to simplify :\n", val, "\n\n"])
                                fez.append(val)
                                sol = solve(val, spc)
                                if isinstance(sol, dict):
                                    sol = sol[spc]
                                else:
                                    sol = simplify(sol[0])
                                good.append([spc, sol])
                                var_sp = str(spc)
                                answer = str(sol)
                                ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                                if var_sp not in my_equations:
                                    my_equations[var_sp] = answer
                            else:
                                not_f.append([xvar[0], xvar[1]])
                        except:
                            not_f.append([xvar[0], xvar[1]])
                        c_c = c_c + 1

                    ffprint(["\n", "Attempting to simplify", "\n\n"])
                    # print(fez)
                    sol = solve(fez, j_s)
                    if isinstance(sol, dict):
                        for c_c in sol:
                            var_sp = str(c_c)
                            answer = str(sol[c_c])
                            ffprint(["\n", var_sp, " = ", answer, "\n\n"])
                            if var_sp not in my_equations:
                                my_equations[var_sp] = answer
                    else:
                        for c_c in sol:
                            rval = Matrix(c_c).subs(equivals)
                            if np.sum(np.array(rval) < 0) > 0:
                                pass
                            else:
                                x_sol = [c_c[xvar] for xvar in range(len(j_s))]
                                ind = 0
                                for xvar in x_sol:
                                    var_sp = str(j_s[ind])
                                    answer = str(simplify(xvar))
                                    ffprint(
                                        ["\n", var_sp, " = ", answer, "\n\n"])
                                    if var_sp not in my_equations:
                                        my_equations[var_sp] = answer
                                    ind = ind + 1
                                break
    # print(slabels)
    return my_equations
