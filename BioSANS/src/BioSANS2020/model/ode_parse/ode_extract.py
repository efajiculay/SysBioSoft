"""

                     This is the ode_extract module

This module grabs stoichiometric matrix, propensity vector, and chemical
reaction from a list of strings containing python expression

List of functions:

1. process
2. prop_extr
3. term_ext
4. get_prop_stoich
5. print_stoich_prop
6. grab_rxn_ksn
7. transform_to_rxn
8. odedxdt_to_topo

For a given set of chemical species and corresponding ODE expression,
this script can be used as follows;

1.	Declare the list of ODE and variables or grab the list from a file.

    x = ["m","P"]
    dx_dt = ["km - m*rm","kp*m - P*rp"]

2.	Use "transform_to_rxn" command to get reactions that match the ODE.

    transform_to_rxn(x, dx_dt)

The result of this command is as follows;

    0 NONE => 1.0 m, 1 ::::: lambda   : km
    1.0 m => 0 NONE, 1 ::::: lambda m : m*rm
    1.0 P => 0 NONE, 1 ::::: lambda P : P*rp
    0 NONE => 1.0 P, 1 ::::: lambda m : kp*m

The first line corresponds to the formation of mRNA m with a rate
constant of km. The second line is degradation of m with degradation
constant rm. The third line is degradation of protein with degradation
constant rp. The last line is the formation of P with rate constant kp.
The lambda expression after the delimiter ':::::' tells us the
propensity for each reaction which is provided after ':'. The 0 NONE
means there is no defined reactant and or degradation product for this
reaction.

3.	Use "print_stoich_prop" to print stoichiometric matrix and
    propensity vector.

        print_stoich_prop(dx_dt)

    The result is as follows;

    1.0*km - 1.0*m*rm
    -1.0*P*rp + 1.0*kp*m

    [1.00000000000000, -1.00000000000000, 0, 0]
    [0, 0, -1.00000000000000, 1.00000000000000]

    [km]
    [m*rm]
    [P*rp]
    [kp*m]

The first two line here is the declared ODE. After the space, is the
stoichiometric matrix. The rows corresponds to the species declared in
step 1 i.e. first row in the matrix corresponds to x[0] or m and the
second row corresponds to x[1] or P. The columns correspond to the terms
in the propensity as listed after the stoichiometric matrix.


"""

from sympy import sympify, expand, Matrix
import numpy as np
from BioSANS2020.gui_functs.scrollable_text import INSERT
from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text


DONE_PARSING = set()


def process(xvar):
    """This function formats the expression xvar into a form that can be
    easily processed by stripping and replacing invalid operators.
    Args:
        xvar : mathematical expresson i.e. 'A*B*ka/(C**2 + 1)'
    Returns:
        val : same as xvar but with invalid operator removed
    """
    # global DONE_PARSING

    # stripping leading and trailing * operator
    # replacing */ with /
    xxvar = xvar.strip("*").strip().replace("*/", "/")
    if xxvar[0] == "/":
        xxvar = "1" + xxvar
    if xxvar[-1] == "/":
        xxvar = xxvar.strip("/")
    val = sympify(xxvar.strip("*"))
    if val in DONE_PARSING:
        return None
    # else:
    DONE_PARSING.add(val)
    return val


def prop_extr(expr, prop):
    """This function extracts the propensity terms from the expression
    expr which is a string of python mathematical expression. Propensity
    terms are the terms that can be separated by + or - without the
    numerical coefficient or stoichiometric multiplier.
    Args:
        expr : mathematical expresson i.e. 'A*B*ka/(C**2 + 1)'
        prop : propensity list to append extracted propensity term
    """
    ex = str(expr)  # transform to string to ensure it is string
    open_p = 0      # number of opening parenthesis
    close_p = 0     # number of closing parenthesis
    diff = 0        # difference between open_p and close_p

    collect = ""    # variable to collect propensity terms
    # iteration in each character of the string exprression ex
    for vstr in ex:
        if vstr in ["+", "-"] and diff == 0:
            if len(collect) > 0:
                dvar = process(collect)
                if dvar is not None:
                    prop.append(dvar)
                collect = ""
                open_p = 0
                close_p = 0
        elif vstr == "(":
            open_p = open_p + 1
            collect = collect + vstr
        elif vstr == ")":
            close_p = close_p + 1
            collect = collect + vstr
        elif (vstr.isnumeric() or vstr == ".") and (diff == 0):
            if len(collect) >= 3:
                if collect[-2:] == "**" or collect[-1] not in ["*", "/"]:
                    collect = collect + vstr
            elif len(collect) >= 1:
                if collect[-1] not in ["*", "/", " "]:
                    collect = collect + vstr
            else:
                pass
        else:
            collect = collect + vstr
        diff = open_p - close_p

    dvar = process(collect)
    if dvar is not None:
        prop.append(dvar)


def term_ext(expr):
    """This function extracts the terms from the expression expr which
    is a string of python mathematical expression. The terms are part of
    expr which can be separated by + or - with the numerical coefficient
    or stoichiometric multiplier.
    Args:
        expr : mathematical expresson i.e. 'A*B*ka/(C**2 + 1)'
    Returns:
        term : the list of extracted terms
    """
    term = []
    ex = str(expr)
    open_p = 0
    close_p = 0
    diff = 0

    collect = ""
    last_sign = ""
    for vstr in ex:
        if vstr in ["+", "-"] and diff == 0:
            if len(collect) > 0:
                term.append(last_sign + collect)
                collect = ""
                open_p = 0
                close_p = 0
            last_sign = vstr
        elif vstr == "(":
            open_p = open_p + 1
            collect = collect + vstr
        elif vstr == ")":
            close_p = close_p + 1
            collect = collect + vstr
        else:
            collect = collect + vstr
        diff = open_p - close_p

    term.append(last_sign + collect)
    return term


def get_prop_stoich(dxdt):
    """This function extracts the propensity vector and stoichiometric
    matrix from the list of ordinary differential equation.
    Args:
        dxdt : list of strings of mathematical expression
    Returns :
        v_stoich : stoichiometric matrix (sympy Matrix)
        w_var : propensity vector ( sympy type Matrix )
    """
    prop = []
    da_dt = []
    for expro in dxdt:
        expr = expand(sympify(expro))
        prop_extr(expr, prop)
        da_dt.append(expr)

    w_var = prop
    v_stoich = [[0 for xvar in range(len(w_var))] for y in range(len(da_dt))]

    for i, _ in enumerate(da_dt):
        for j, _ in enumerate(w_var):
            s_var = 0
            for xxvar in term_ext(da_dt[i]):
                xvar = sympify(xxvar) / w_var[j]
                try:
                    s_var = s_var + float(xvar)
                except BaseException:
                    pass
            v_stoich[i][j] = s_var

    return Matrix(v_stoich), Matrix(w_var)


def print_stoich_prop(dxdt):
    """This function extracts the propensity vector and stoichiometric
    matrix from the list of ordinary differential equation and prints
    the output in the console.
    Args:
        dxdt : list of strings of mathematical expression
    """
    global DONE_PARSING
    DONE_PARSING = set()
    print()
    v_stoich, w_var = get_prop_stoich(dxdt)
    for tvar in v_stoich * w_var:
        print(tvar)
    print()

    for c_var in np.array(v_stoich):
        print([round(y, 4) for y in c_var])
    print()

    for c_var in np.array(w_var):
        print(c_var)


def grab_rxn_ksn(stch_var, xvar, w_var):
    """This function transform the stoichiometric matrix, species or
    components, and propensity vector into a list of reactions and list
    of rate constants.
    Args:
        stch_var : stoichiometric matrix or 2D Matrix of coefficient
        xvar     : list of components or species
        w_var    : propensity vector or 1D Matrix of fluxes
    Returns:
        rxn_var  : list of reaction similar to BioSANS reactions format
        ksn_var  : set of rate constant symbols if not numeric
    """
    rxn_var = []
    ksn_var = set()
    ind = 0

    for col in stch_var.T:
        r_big = ""
        p_big = ""
        for i, _ in enumerate(col):
            if col[i] != 0:
                if col[i] < 0:
                    r_big += str(abs(col[i])) + " " + xvar[i] + " " + "+ "
                else:
                    p_big += str(abs(col[i])) + " " + xvar[i] + " " + "+ "
        if r_big.strip() == "":
            r_big = "0 NONE"
        if p_big.strip() == "":
            p_big = "0 NONE"

        in_sp = []
        for s_var in w_var[ind].free_symbols:
            ss_var = str(s_var)
            if ss_var in xvar:
                in_sp.append(ss_var)
            else:
                ksn_var.add(ss_var)
        in_sp = ",".join(in_sp)
        rxn_var.append(r_big.strip("+ ") + " => " + p_big.strip("+ ") +
                       ", 1 ::::: lambda " + in_sp + " : " + str(w_var[ind]))
        ind = ind + 1

    return rxn_var, ksn_var


def transform_to_rxn(xvar, dxdt, x_ini=None, k_rc=None, items=None):
    """This function transform the list of components and list of ODE
    into BioSANS topology file format.
    Args:
        xvar     : list of components or species (variable in ODE)
                   example ["m","P"]
        dxdt     : list of differential equations (python string)
                   example ["km - m*rm","kp*m - P*rp"]
        x_ini    : dictionary of initial conecentration or value
                   example {"m" : 10, "p" : 0}
        k_rc     : dictionary of rate constant symbols abd values
                   example {"km" : 0.1, "rm" : 0.2}
        items    : list containing [canvas, scroll_x, scroll_y]
    Returns:
        text     : text area where the outputs are written
    """
    global DONE_PARSING
    # print(xvar, dxdt, x_ini, k_rc, items, sep="\n\n")
    if items:
        text = prepare_scroll_text(items)

        def ffrint(xvar):
            return text.insert(INSERT, xvar + "\n")
    else:
        def ffrint(xvar):
            return print(" ".join([str(y) for y in xvar]), end="")

    if x_ini is None:
        x_ini = {}
    if k_rc is None:
        k_rc = {}

    DONE_PARSING = set()
    v_stoich, w_var = get_prop_stoich(dxdt)
    stch_var = np.around(np.array(v_stoich).astype(float), 3)
    rxn_var, ksn_var = grab_rxn_ksn(stch_var, xvar, w_var)

    ffrint("Function_Definitions:")
    for k in ksn_var:
        if k.strip() not in k_rc:
            ffrint(k.strip() + " = type actual value")
        else:
            ffrint(k.strip() + " = " + k_rc[k.strip()])
    for c_var in xvar:
        if c_var.strip() not in x_ini:
            ffrint(c_var.strip() + "_ini = type actual value")
        else:
            ffrint(c_var.strip() + "_ini = " + x_ini[c_var.strip()])

    ffrint("")
    ffrint("#REACTIONS")
    for rx_i in rxn_var:
        ffrint(rx_i)

    ffrint("")
    ffrint("@CONCENTRATION")
    for c_var in xvar:
        ffrint(c_var + " , " + c_var.strip() + "_ini")
    ffrint("NONE, 1")
    return text


def odedxdt_to_topo(mfile, items):
    """This function reads BioSANS ODE file format and converts it to
    BioSANS topology file format. The ODE file contains ODE expression,
    initial concentration, and rate constants separated by tags. There
    are three tags that is currently supported : ODE_DECLARATIONS,
    INI_CONCENTRATIONS, and RATE_CONSTANTS. Provided below is an example
    of the content of a typical ODE file.

        ODE_DECLARATIONS:
        A = -ka*A*B/(1+C**2) + kf1/(1+B**2)
        B = -ka*A*B/(1+C**2)
        C = -kc*C + kf2
        D = ka*A*B/(1+C**2) - kf2

        INI_CONCENTRATIONS:
        A = 100
        B = 200
        C = 150
        D = 0

        RATE_CONSTANTS:
        ka = 0.02
        kf1 = 0.2
        kc = 0.03
        kf2 = 0.01

    Args:
        mfile    : file containing list of ODE expression
        items    : list containing [canvas, scroll_x, scroll_y]
    Returns:
        text     : text area where the outputs are written
    """
    ddvar = open(mfile, "r")
    print()
    xvar = []
    x_ini = {}
    k_rc = {}
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
                x_ini[row[0].strip()] = row[1].strip()
        elif last == "RATE_CONSTANTS" and xxvar.strip() \
                not in ["ODE_DECLARATIONS:", "INI_CONCENTRATIONS:"]:
            if xxvar.strip() != "":
                row = xxvar.split("=")
                k_rc[row[0].strip()] = row[1].strip()
        elif xxvar.strip() == "ODE_DECLARATIONS:":
            last = "ODE_DECLARATIONS"
        elif xxvar.strip() == "INI_CONCENTRATIONS:":
            last = "INI_CONCENTRATIONS"
        elif xxvar.strip() == "RATE_CONSTANTS:":
            last = "RATE_CONSTANTS"

    # print_stoich_prop(dxdt)
    ddvar.close()
    return transform_to_rxn(xvar, dxdt, x_ini, k_rc, items)

# odedxdt_to_topo("NewText.txt")
