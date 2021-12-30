"""

                   This module is the lna_approx2 module

This module handles symbolic linear noise approximation.

The following are the list of function in this module;

1. subs2
2. lna_symbolic2

"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from sympy import simplify, Matrix, solve, Symbol, flatten, nonlinsolve, \
    collect, factor
from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.propagation.propensity \
    import propensity_vec, propensity_vec_molar


def subs2(zvar, cval):
    """This function helps in the substitution of value to sympy Symbols

    Args:
        zvar (Symbol): sympy expression
        cval (dict): dictionary of values

    Returns:
        (Symbol): substituted expression
    """
    for xvar in cval:
        zvar = zvar.subs(xvar, cval[xvar])
    return zvar


def lna_symbolic2(sp_comp, ks_dict, conc, r_dict, p_dict, stch_var,
                  items=None, molar=False, mode=None):
    """This function facilitates in the symbolic LNA computation.

    Args:
        sp_comp (dict): dictionary of appearance or position of species
            or component in the reaction tag of BioSANS topology file.

            For example;

                #REACTIONS
                A => B, kf1
                B => C, kf2

            The value of sp_comp is

                sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}

                A appears in first reaction with index 0
                B appears in first and second reaction with index 0, 1
                C appears in second reaction with index 1
        ks_dict (dict): dictionary of rate constant that appears in each
            reactions.

            For example;

                #REACTIONS
                A => B , 0.3        # first reaction
                B <=> C, 0.1, 0.2   # second reaction

            The value of ks_dict is

                ks_dict = {
                    0 : [0.3],      # first reaction
                    1 : [0.1, 0.2]  # second reaction
                }
        conc ([type]): [description]
        r_dict (dict): dictionary of reactant stoichiometry. For example

            r_dict = {
                0: {'A': 1},  # first reaction, coefficient of A is 1
                1: {'B': 1}   # second reaction, coefficient of B is 1
            }
        p_dict (dict): dictionary of product stoichiometry. For example

            p_dict = {
                0: {'B': 1},  # first reaction, coefficient of B is 1
                1: {'C': 1}   # second reaction, coefficient of C is 1
            }
        stch_var (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        items (list): list of [canvas, scroll_x, scroll_y]
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.
        mode (str, optional): method keywords : Numeric, fofks, fofCo

    Returns:
        list: [0, 0] - not used
    """

    c_s = {}
    c_so = {}
    equivals = []
    equi_co = []
    equi_ks = []

    for xvar in sp_comp:
        c_s[xvar] = Symbol(xvar, real=True, negative=False)
        c_so[xvar] = Symbol(xvar + 'o', real=True, negative=False) * \
            (0 if conc[xvar] == 0 else 1)
        equivals.append((c_so[xvar], conc[xvar]))
        equi_co.append((c_so[xvar], conc[xvar]))

    kcs = []
    for i, _ in enumerate(ks_dict):
        row = []
        if len(ks_dict[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[0], ks_dict[i][0]))
            equi_ks.append((row[0], ks_dict[i][0]))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[0], ks_dict[i][0]))
            equi_ks.append((row[0], ks_dict[i][0]))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, real=True, negative=False))
            equivals.append((row[1], ks_dict[i][1]))
            equi_ks.append((row[1], ks_dict[i][1]))
        kcs.append(row)

    if not molar:
        prop_flux = Matrix(propensity_vec(kcs, c_s, r_dict, p_dict))
    else:
        prop_flux = Matrix(propensity_vec_molar(kcs, c_s, r_dict, p_dict))

    if mode == "Numeric":
        prop_flux = prop_flux.subs(equivals)
    elif mode == "fofks":
        prop_flux = prop_flux.subs(equi_co)
    elif mode == "fofCo":
        prop_flux = prop_flux.subs(equi_ks)

    stch_var = Matrix(stch_var)
    # c_s might have change after call
    for xvar in sp_comp:
        c_s[xvar] = Symbol(xvar, real=True, negative=False)
    slabels = [xvar for xvar in c_s]

    s_s = []
    n_z = []
    for row in range(stch_var.shape[0]):
        if sum(abs(stch_var[row, :])) != 0 and slabels[row][0] != "-":
            s_s.append(list(stch_var[row, :]))
            n_z.append(row)
    s_s = Matrix(s_s)

    da_dt = s_s * prop_flux
    # print(da_dt)
    ccs = [c_s[xvar] for xvar in c_s]
    j_s = [ccs[xvar] for xvar in n_z]
    a_jac = da_dt.jacobian(j_s)
    diag_prof = [[0] * len(prop_flux) for xvar in range(len(prop_flux))]
    for i, _ in enumerate(prop_flux):
        diag_prof[i][i] = diag_prof[i][i] + prop_flux[i]
    diag_prof = Matrix(diag_prof)

    bbt = s_s * diag_prof * s_s.T

    cov = []
    for i in n_z:
        row = []
        for j in n_z:
            key = "C" + str(min(i + 1, j + 1)) + "_" + str(max(i + 1, j + 1))
            if i == j:
                row.append(Symbol(key, real=True, negative=False))
            else:
                row.append(Symbol(key, real=True))
        cov.append(row)
    cov = Matrix(cov)
    # print(s_s)
    same = {}
    for i in range(s_s.shape[0] - 1):
        for j in range(i + 1, s_s.shape[0]):
            if s_s[i, :] == s_s[j, :]:
                same[(i, j)] = 1
            elif s_s[i, :] == -s_s[j, :]:
                same[(i, j)] = -1
            else:
                pass

    dcov_dt = a_jac * cov + cov * a_jac.T + bbt

    x_s = list(dict.fromkeys(flatten(cov)))

    if items:
        text = prepare_scroll_text(items)

        def ffprint(xvar):
            return text.insert('insert', " ".join([str(y) for y in xvar]))
    else:

        def ffprint(xvar):
            return print(" ".join([str(y) for y in xvar]), end="")

    sps = [xvar for xvar in sp_comp]
    csps = {}
    for i in n_z:
        for j in n_z:
            if j >= i:
                key = "C" + str(min(i + 1, j + 1)) + "_" + \
                    str(max(i + 1, j + 1))
                ffprint([key, " = ", str(sps[i]) + "_" + str(sps[j]), "\n"])
                csps[key] = sps[i] + "_" + sps[j]

    eqs = []
    ffprint(["\nEquations to solve\n\n"])
    for i in range(len(n_z)):
        for j in range(i, len(n_z)):
            if dcov_dt[i, j] != 0:
                row = collect(dcov_dt[i, j], x_s)
                eqs.append(row)
                ffprint([row, "\n"])

    for e_n in same:
        i, j = e_n
        for k in range(j, s_s.shape[0]):
            eqs.append(cov[j, k] - same[e_n] * cov[i, k])

    post_proc = False
    multiple_cval = False

    val = solve(da_dt, j_s)

    cval = {}
    if val:
        for xvar in j_s:
            key = xvar
            if key in val:
                cval[key] = val[key]
            else:
                post_proc = True
    else:
        post_proc = True

    if post_proc:
        ffprint(
            ["\nSteady state concentrations not properly calculated",
             "due to incomplete constraint\n\n"])
        ffprint(["\n", da_dt, "\n"])
        ffprint(["\nAttempting the use of boundary condition", "\n"])

        tosum = set()
        sets = [tosum]
        used = set()
        for ihvar, _ in enumerate(r_dict):
            used2 = set()
            if len(r_dict[ihvar]) == 1:
                for svar in sets:
                    key = list(r_dict[ihvar].keys())[0]
                    if key not in used and r_dict[ihvar][key] != 0:
                        svar.add(key)
                        used2.add(key)
            else:
                for j, _ in enumerate(sets):
                    svar = sets[j]
                    sets.append(svar.copy())
                    key = list(r_dict[ihvar].keys())[0]
                    if key not in used and r_dict[ihvar][key] != 0:
                        svar.add(key)
                        used2.add(key)
                    key = list(r_dict[ihvar].keys())[1]
                    if key not in used and r_dict[ihvar][key] != 0:
                        sets[-1].add(key)
                        used2.add(key)

            if len(p_dict[ihvar]) == 1:
                for svar in sets:
                    key = list(p_dict[ihvar].keys())[0]
                    if key not in used and p_dict[ihvar][key] != 0:
                        svar.add(key)
                        used2.add(key)
            else:
                for j, _ in enumerate(sets):
                    svar = sets[j]
                    sets.append(svar.copy())
                    key = list(p_dict[ihvar].keys())[0]
                    if key not in used and p_dict[ihvar][key] != 0:
                        svar.add(key)
                        used2.add(key)
                    key = list(p_dict[ihvar].keys())[1]
                    if key not in used and p_dict[ihvar][key] != 0:
                        sets[-1].add(key)
                        used2.add(key)
            for z_s in used2:
                used.add(z_s)

        f_e = [da_dt]
        for svar in sets:
            f_e.append(sum([c_s[xvar] - c_so[xvar] for xvar in svar]))

        if mode == "Numeric":
            f_e = [entry.subs(equivals) for entry in f_e]
        elif mode == "fofks":
            f_e = [entry.subs(equi_co) for entry in f_e]
        elif mode == "fofCo":
            f_e = [entry.subs(equi_ks) for entry in f_e]

        val2 = solve(f_e, {xvar for xvar in j_s})
        cval = {}
        for xvar in j_s:
            key = xvar
            if key in val2:
                cval[key] = val2[key]
            else:
                multiple_cval = True
                break
        # print(f_e)

    ffprint(["\nUsing Algebraic manipulation of AC + CA.T + BT = 0\n"])
    if not multiple_cval:
        print(1)
        ffprint(["\nSteady state concentrations\n\n"])
        for xvar in cval:
            ffprint([xvar, " = ", cval[xvar], "\n\n"])
        dcov_dt = simplify(subs2(dcov_dt, cval))
        eqs = []
        ffprint(["\nEquations to solve\n\n"])
        for i in range(len(n_z)):
            for j in range(i, len(n_z)):
                if dcov_dt[i, j] != 0:
                    row = collect(dcov_dt[i, j], x_s)
                    eqs.append(row)
                    ffprint([row, "\n\n"])

        for e_n in same:
            i, j = e_n
            for k in range(j, s_s.shape[0]):
                eqs.append(cov[j, k] - same[e_n] * cov[i, k])

        sol = solve(eqs, x_s)

        if not sol:
            # LA = Matrix(eqs)
            # LHS = LA.jacobian(x_s)
            # RHS = LA-LHS*Matrix(x_s)
            # sol = simplify(-(LHS**-1)*RHS)
            sol = list(nonlinsolve(eqs, x_s))
            sol = {x_s[xvar]: sol[xvar] for xvar in range(len(x_s))}
        ffprint(["\nCovariance\n\n"])
        for xvar in sol:
            c_h = factor(simplify(subs2(sol[xvar], cval)))
            # c_h = sol[xvar]
            ffprint(["Cov(" + str(csps[str(xvar)]) + ")", " = ", c_h, "\n\n"])
    else:
        print(2)
        ffprint(["\nMultiple solutions detected"])
        valset = 0
        for val in val2:
            ffprint(["\nValues ", valset, "\n"])
            cval = {}
            for xvar in j_s:
                key = xvar
                if key in val:
                    cval[key] = val[key]
            ffprint(["\nSteady state concentrations\n\n"])
            for xvar in cval:
                ffprint([xvar, " = ", cval[xvar], "\n\n"])
            dcov_dt = simplify(subs2(dcov_dt, cval))
            eqs = []
            ffprint(["\nEquations to solve\n\n"])
            for i in range(len(n_z)):
                for j in range(i, len(n_z)):
                    if dcov_dt[i, j] != 0:
                        row = collect(dcov_dt[i, j], x_s)
                        eqs.append(row)
                        ffprint([row, "\n\n"])

            for e_n in same:
                i, j = e_n
                for k in range(j, s_s.shape[0]):
                    eqs.append(cov[i, k] - same[e_n] * cov[j, k])

            sol = solve(eqs, x_s)
            if not sol:
                sol = list(nonlinsolve(eqs, x_s))
                sol = {x_s[xvar]: sol[0][xvar] for xvar in range(len(x_s))}
            ffprint(["\nCovariance\n\n"])
            fact = 1
            for xvar in sol:
                uvar, wvar = str(xvar).split("_")
                c_h = factor(simplify(subs2(sol[xvar], cval))).subs(equivals)
                if uvar == 'C' + wvar and c_h < 0:
                    fact = -1
                    break

            for xvar in sol:
                c_h = factor(simplify(subs2(sol[xvar], cval)))
                # c_h = sol[xvar]
                ffprint(["Cov(" + str(csps[str(xvar)]) + ")",
                         " = ", c_h * fact, "\n\n"])
            valset = valset + 1

    return [0, 0]
