# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from scipy import linalg as LA
from scipy.optimize import fsolve

from BioSANS2020.propagation.propensity import propensity_vec_molar
from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text


def rem_rowcol_zero(a_mat):
    return a_mat[:, ~np.all(a_mat == 0, axis=0)][~np.all(a_mat == 0, axis=1)]


def lna_ss_jacobian(model, zlist, sp_comp, stch_var,
                    ks_dict, r_dict, p_dict):
    jacbn = np.zeros((len(zlist), len(zlist)))
    div = 1
    new = 1000
    old = 2000
    check = 1000
    while check > 1.0e-7:
        for i, _ in enumerate(zlist):
            h_var = abs(zlist[i] * 1.0e-8) / div
            h2_var = 2 * h_var
            zlist[i] = zlist[i] - h_var
            ini_val = model(zlist, sp_comp, ks_dict, r_dict,
                            p_dict, stch_var)
            zlist[i] = zlist[i] + 2 * h_var
            out_val = model(zlist, sp_comp, ks_dict, r_dict,
                            p_dict, stch_var)
            zlist[i] = zlist[i] - h_var
            for j, _ in enumerate(out_val):
                jacbn[j, i] = (out_val[j] - ini_val[j]) / h2_var
        div = div * 2
        old = new
        new = sum([jacbn[k, k] for k in range(len(zlist))])
        check = abs(new - old) / abs(new)
    return np.array(jacbn)


def lna_model_ss(zlist, sp_comp, ks_dict, r_dict, p_dict, stch_var):
    conc = {}
    ind = 0
    for spi in sp_comp:
        conc[spi] = zlist[ind]
        ind = ind + 1
    prop_flux = propensity_vec_molar(
        ks_dict, conc, r_dict, p_dict, True)
    fofx = np.matmul(stch_var, prop_flux)
    return fofx.reshape(len(sp_comp))


def lna_steady_state(t_var, sp_comp, ks_dict, conc, r_dict, p_dict,
                     stch_var, items=None):
    ind = 0
    zlist = []
    for spi in sp_comp:
        zlist.append(conc[spi])
    zlist = fsolve(
        lna_model_ss,
        tuple(zlist),
        xtol=1.0e-10,
        args=(sp_comp, ks_dict, r_dict, p_dict, stch_var))
    ind = 0
    for spi in sp_comp:
        conc[spi] = zlist[ind]
        ind = ind + 1
    dsf_dx = lna_ss_jacobian(
        lna_model_ss, zlist, sp_comp, stch_var, ks_dict, r_dict, p_dict)
    f_prop = propensity_vec_molar(ks_dict, conc, r_dict, p_dict, True)
    bb_diffsn = np.matmul(
        np.matmul(
            stch_var,
            np.diag(
                f_prop.flatten())),
        stch_var.T)

    dsf_dx = np.nan_to_num(dsf_dx)
    bb_diffsn = np.nan_to_num(bb_diffsn)
    cov_mat = LA.solve_continuous_lyapunov(dsf_dx, -bb_diffsn)

    if items:
        text = prepare_scroll_text(items)

        def fprint(xvar):
            return text.insert('insert', " ".join([str(y) for y in xvar]))
    else:
        def fprint(xvar):
            return print(" ".join([str(y) for y in xvar]), end="")

    fprint(["\nConcentrations\n\n"])
    ind = 0
    for spi in sp_comp:
        if spi[0] != "-":
            fprint([spi, " = ", zlist[ind], "\n"])
            ind = ind + 1

    fprint(["\nCovariance\n\n"])
    i_ind = 0
    for spi in sp_comp:
        j = 0
        for spj in sp_comp:
            if j >= i_ind:
                val = cov_mat[i_ind, j]
                if str(val) not in {"None", "nan", "0.0"} \
                        and spi[0] != "-" and spj[0] != "-":
                    fprint([" ".join(["Covr", spi, spj])
                            .ljust(50), "=", val, "\n"])
            j = j + 1
        i_ind = i_ind + 1
    fprint(["\nCorrelation\n\n"])
    i_ind = 0
    for spi in sp_comp:
        j = 0
        for spj in sp_comp:
            if j >= i_ind:
                val = cov_mat[i_ind, j] / \
                    np.sqrt(np.abs(cov_mat[i_ind, i_ind]
                                   * cov_mat[j, j]))
                if str(val) not in {"None", "nan", "0.0"} \
                        and spi[0] != "-" and spj[0] != "-":
                    fprint([" ".join(["Corr", spi, spj]).ljust(50),
                            "=", val, "\n"])
            j = j + 1
        i_ind = i_ind + 1

    fprint(["\nFano Factor\n\n"])
    ind = 0
    for spi in sp_comp:
        val = cov_mat[ind, ind] / zlist[ind]
        if str(val) not in {"None", "nan", "0.0"} and spi[0] != "-":
            fprint([" ".join(["Fano Factor for", spi])
                    .ljust(50), "=", val, "\n"])
        ind = ind + 1

    return [cov_mat, zlist]
