from numpy import array as np_array
from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT
from BioSANS2020.propagation.propensity \
    import propensity_vec, propensity_vec_molar
from sympy import Symbol, Function
    
def eval_dict(to_eval, loc_dict):
    """This function takes a string expression and return the evaluated
    expression using SBML_FUNCT_DICT and the locals() dictionary where
    eval_dict is called.

    Args:
        to_eval (str): the expression to evaluate
        loc_dict (dict): local dictionary from the calling function

    Returns:
        multitype: result of eval command
    """
    return eval(to_eval, loc_dict, SBML_FUNCT_DICT)


def tofloat(val, loc_dict):
    """This function attempts to convert the input val into float

    Args:
        val (str): the expression to evaluate
        loc_dict (dict): local dictionary from the calling function

    Returns:
        float: float equivalent of val
    """
    try:
        return float(val)
    except:
        return float(eval_dict(val, loc_dict))


def is_number(xvar):
    """This function checks if a string xvar is float

    Args:
        xvar (str): input string expression or number

    Returns:
        bool: True if xvar cal be converted to float otherwise False
    """
    try:
        float(xvar)
        return True
    except:
        return False
        
        
def read_topo_data(rfile):
    rows = []
    conc = {}
    with open(rfile, "r") as file:
        last = ""
        for row in file:
            row = row.strip()+" "
            if last == "Function_Definitions":
                if row.strip() != "" and row[0] != "#":
                    exec(row.strip(), locals(), SBML_FUNCT_DICT)
                elif row[0] == "#":
                    last = "#"
            elif last == "#":
                if row.strip() != "" and row[0] != "@":
                    rows.append(row)
                elif row.strip() != "" and row[0] != "#":
                    last = "@"
            elif last == "@":
                if row.strip() != "" and row[0] != "@":
                    cvar = row.split(",")
                    conc[cvar[0].strip()] = tofloat(cvar[1], locals())
            elif row[0] == "#":
                last = "#"
            elif row[0] == "@":
                last = "@"
            elif row.strip().upper() == "FUNCTION_DEFINITIONS:":
                last = "Function_Definitions"
    return rows, conc


def symbol_conv(sp_comp, ks_dict, r_dict, p_dict):

    c_s = {}
    tvar = Symbol('t', negative=False, real=True)

    cs_r = {}
    for xvar in sp_comp:
        c_s[xvar] = Function(xvar, negative=False, real=True)
        c_s[xvar] = c_s[xvar](tvar)
        cs_r[xvar] = Symbol(xvar, negative=False, real=True)

    kcs = []
    for i, _ in enumerate(ks_dict):
        row = []
        if len(ks_dict[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, negative=False, real=True))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, negative=False, real=True))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, negative=False, real=True))
        kcs.append(row)

    return kcs, cs_r


def dict_from_rows(rows, conc, algeb=False):
    ks_dict = {}
    r_dict = {}
    p_dict = {}
    sp_comp = {}
    rxn_rows = len(rows)
    for ih_ind in range(rxn_rows):
        r_dict[ih_ind] = {}
        p_dict[ih_ind] = {}
        col_row = rows[ih_ind].split(":::::")
        row = col_row[0].strip().split(",")

        if len(row) == 3:
            ks_dict[ih_ind] = [
                tofloat(row[1], locals()),
                tofloat(row[2], locals())]
        else:
            ks_dict[ih_ind] = [tofloat(row[1], locals())]
        col_var = row[0].split("<=>")
        if len(col_var) == 1:
            col_var = row[0].split("=>")

        sp_c = col_var[0]
        svar = sp_c.strip().split()
        if len(svar) > 1:
            last = 1
            for xvar in svar:
                if not is_number(xvar) and xvar != "+":
                    r_dict[ih_ind][xvar] = last
                    last = 1
                    if xvar in sp_comp:
                        sp_comp[xvar].add(ih_ind)
                    else:
                        sp_comp[xvar] = {ih_ind}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = svar[0]
            r_dict[ih_ind][xvar] = 1
            if xvar in sp_comp:
                sp_comp[xvar].add(ih_ind)
            else:
                sp_comp[xvar] = {ih_ind}

        sp_c = col_var[1]
        svar = sp_c.strip().split()
        if len(svar) > 1:
            last = 1
            for xvar in svar:
                if (not is_number(xvar) or xvar.lower() == "e") \
                        and xvar != "+":
                    p_dict[ih_ind][xvar] = last
                    last = 1
                    if xvar in sp_comp:
                        sp_comp[xvar].add(ih_ind)
                    else:
                        sp_comp[xvar] = {ih_ind}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = svar[0]
            p_dict[ih_ind][xvar] = 1
            if xvar in sp_comp:
                sp_comp[xvar].add(ih_ind)
            else:
                sp_comp[xvar] = {ih_ind}

    stoch_var = []

    for sp_c in sp_comp:
        row = []
        for r_ind in range(rxn_rows):
            prod = p_dict[r_ind][sp_c] if sp_c in p_dict[r_ind] else 0
            rect = r_dict[r_ind][sp_c] if sp_c in r_dict[r_ind] else 0
            row.append(prod - rect)
            if len(ks_dict[r_ind]) == 2:
                row.append(-row[-1])
        stoch_var.append(row)
    stoch_var = np_array(stoch_var)
        
    if algeb:
        ks_dict, conc = symbol_conv(sp_comp, ks_dict, r_dict, p_dict)
    return ks_dict, r_dict, p_dict, sp_comp, stoch_var, conc