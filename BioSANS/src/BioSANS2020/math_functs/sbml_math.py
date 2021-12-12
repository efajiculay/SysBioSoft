"""

                  This module in the sbml_math module

This contains function needed to interpret SBML files.

"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import sympy
import numpy as np

# def abs(xvar):
# return sympy.fabs(xvar)
# need to use trigonometric relationship to create non-built in numpy 
# trigonometric functions

NUMBER_TYPE = (float, int)


def acos(xvar):
    """returns sympy.acos(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arccos(xvar)
    return sympy.acos(xvar).evalf()


def arccos(xvar):
    """returns sympy.acos(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arccos(xvar)
    return sympy.acos(xvar).evalf()


def acosh(xvar):
    """returns sympy.acosh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arccosh(xvar)
    return sympy.acosh(xvar).evalf()


def arccosh(xvar):
    """returns sympy.acosh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arccosh(xvar)
    return sympy.acosh(xvar).evalf()


def acot(xvar):
    """returns sympy.acot(xvar).evalf()"""
    return sympy.acot(xvar).evalf()


def arccot(xvar):
    """returns sympy.acot(xvar).evalf()"""
    return sympy.acot(xvar).evalf()


def acoth(xvar):
    """returns sympy.acoth(xvar).evalf()"""
    return sympy.acoth(xvar).evalf()


def arccoth(xvar):
    """returns sympy.acoth(xvar).evalf()"""
    return sympy.acoth(xvar).evalf()


def acsc(xvar):
    """returns sympy.acsc(xvar).evalf()"""
    return sympy.acsc(xvar).evalf()


def arccsc(xvar):
    """returns sympy.acsc(xvar).evalf()"""
    return sympy.acsc(xvar).evalf()


def acsch(xvar):
    """returns sympy.acsch(xvar).evalf()"""
    return sympy.acsch(xvar).evalf()


def arccsch(xvar):
    """returns sympy.acsch(xvar).evalf()"""
    return sympy.acsch(xvar).evalf()


def arcsec(xvar):
    """returns sympy.asec(xvar).evalf()"""
    return sympy.asec(xvar).evalf()


def asech(xvar):
    """returns sympy.asech(xvar).evalf()"""
    return sympy.asech(xvar).evalf()


def arcsech(xvar):
    """returns sympy.asech(xvar).evalf()"""
    return sympy.asech(xvar).evalf()


def asin(xvar):
    """returns sympy.asin(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arcsin(xvar)
    return sympy.asin(xvar).evalf()


def asinh(xvar):
    """returns sympy.asinh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arcsinh(xvar)
    return sympy.asinh(xvar).evalf()


def arcsinh(xvar):
    """returns sympy.asinh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arcsinh(xvar)
    return sympy.asinh(xvar).evalf()


def arcsin(xvar):
    """returns sympy.asin(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arcsin(xvar)
    return sympy.asin(xvar).evalf()


def atan(xvar):
    """returns sympy.atan(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arctan(xvar)
    return sympy.atan(xvar).evalf()


def arctan(xvar):
    """returns sympy.atan(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arctan(xvar)
    return sympy.atan(xvar).evalf()


def atanh(xvar):
    """returns sympy.atanh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arctanh(xvar)
    return sympy.atanh(xvar).evalf()


def arctanh(xvar):
    """returns sympy.atanh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.arctanh(xvar)
    return sympy.atanh(xvar).evalf()


def ceil(xvar):
    """returns sympy.ceiling(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.ceil(xvar)
    return sympy.ceiling(xvar).evalf()


def ceiling(xvar):
    """returns sympy.ceiling(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.ceil(xvar)
    return sympy.ceiling(xvar).evalf()


def cos(xvar):
    """returns sympy.cos(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.cos(xvar)
    return sympy.cos(xvar).evalf()


def cosh(xvar):
    """returns sympy.cosh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.cosh(xvar)
    return sympy.cosh(xvar).evalf()


def cot(xvar):
    """returns sympy.cot(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return 1/np.tan(xvar)
    return sympy.cot(xvar).evalf()


def coth(xvar):
    """returns sympy.coth(xvar).evalf()"""
    return sympy.coth(xvar).evalf()


def csc(xvar):
    """returns sympy.csc(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return 1/np.sin(xvar)
    return sympy.csc(xvar).evalf()


def csch(xvar):
    """returns sympy.csch(xvar).evalf()"""
    return sympy.csch(xvar).evalf()


def factorial(xvar):
    """returns sympy.factorial(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.math.factorial(xvar)
    return sympy.factorial(xvar).evalf()


def exp(xvar):
    """returns sympy.exp(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.exp(xvar)
    return sympy.exp(xvar).evalf()


def floor(xvar):
    """returns sympy.floor(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.floor(xvar)
    return sympy.floor(xvar).evalf()


def ln(xvar):
    """returns sympy.ln(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.log(xvar)
    return sympy.ln(xvar).evalf()


def log(xvar):
    """returns sympy.log(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.log(xvar)
    return sympy.log(xvar).evalf()


def log10(xvar):
    """returns sympy.log(xvar, 10).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.log10(xvar)
    return sympy.log(xvar, 10).evalf()


def piecewise(*xvar):
    """returns the value before the first True value in xvar. If there
    is no True value, returns the last element in vxar"""
    if len(xvar) == 3:
        return xvar[0] if xvar[1] else xvar[2]

    # ans = []
    xlen = len(xvar) - 1
    for i in range(0, xlen, 2):
        if xvar[i + 1]:
            return xvar[i]
    return xvar[-1]


def pow(xvar, yvar):
    """sympy.Pow(xvar, yvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE) and isinstance(yvar, NUMBER_TYPE):
        return np.float_power(xvar, yvar)
    return sympy.Pow(xvar, yvar).evalf()


def power(xvar, yvar):
    """sympy.Pow(xvar, yvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE) and isinstance(yvar, NUMBER_TYPE):
        return np.float_power(xvar, yvar)
    return sympy.Pow(xvar, yvar).evalf()


def root(nvar, xvar):
    """sympy.root(xvar, nvar).evalf()"""
    return sympy.root(xvar, nvar).evalf()

def sec(xvar):
    """returns sympy.sec(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return 1/np.cos(xvar)
    return sympy.sec(xvar).evalf()


def sech(xvar):
    """sympy.sech(xvar).evalf()"""
    return sympy.sech(xvar).evalf()


def sqr(xvar):
    """sympy.sqrt(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.sqrt(xvar)
    return sympy.sqrt(xvar).evalf()


def sqrt(xvar):
    """sympy.sqrt(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.sqrt(xvar)
    return sympy.sqrt(xvar).evalf()


def sin(xvar):
    """sympy.sin(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.sin(xvar)
    return sympy.sin(xvar).evalf()


def sinh(xvar):
    """sympy.sinh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.sinh(xvar)
    return sympy.sinh(xvar).evalf()


def tan(xvar):
    """sympy.tan(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.tan(xvar)
    return sympy.tan(xvar).evalf()


def tanh(xvar):
    """sympy.tanh(xvar).evalf()"""
    if isinstance(xvar, NUMBER_TYPE):
        return np.tanh(xvar)
    return sympy.tanh(xvar).evalf()


def And(*xvar):
    """returns True if all elements of xvar is True else returns False"""
    for yvar in xvar:
        # if yvar == False:
        if not yvar:
            return False
    return True


def Not(xvar):
    """returns the reverse of the boolean value of xvar"""
    return not xvar


def Or(*xvar):
    """returns True if at least one value in xvar is True else False"""
    for yvar in xvar:
        # if yvar == True:
        if not yvar:
            return True
    return False


def xor(*xvar):
    """returns True if there is odd number of True else returns False"""
    odd = 0
    for yvar in xvar:
        # if yvar == True:
        if yvar:
            odd = odd + 1
    if odd % 2 == 0:
        return False

    return True


def eq(*xvar):
    """returns True if xvar[0] == xvar[1:] else returns False"""
    last = None
    for yvar in xvar:
        if last is None:
            last = yvar
        else:
            if yvar != last:
                return False
    return True


def geq(*xvar):
    """returns True if xvar[0] >= xvar[1:] else returns False"""
    last = None
    for yvar in xvar:
        if last is None:
            last = yvar
        else:
            if last < yvar:
                return False
    return True


def gt(*xvar):
    """returns True if xvar[0] > xvar[1:] else returns False"""
    last = None
    for yvar in xvar:
        if last is None:
            last = yvar
        else:
            if last <= yvar:
                return False
    return True


def leq(*xvar):
    """returns True if xvar[0] <= xvar[1:] else returns False"""
    last = None
    for yvar in xvar:
        if last is None:
            last = yvar
        else:
            if last > yvar:
                return False
    return True


def lt(*xvar):
    """returns True if xvar[0] < xvar[1:] else returns False"""
    last = None
    for yvar in xvar:
        if last is None:
            last = yvar
        else:
            if last >= yvar:
                return False
    return True


def neq(xvar, yvar):
    """returns True if xvar != yvar else returns False"""
    if xvar == yvar:
        return False
    return True


def plus(*xvar):
    """sum(xvar).evalf()"""
    return sum(xvar).evalf()


def times(*xvar):
    """returns the product of all elements in the list xvar"""
    pvar = 1
    for yvar in xvar:
        pvar = pvar * yvar
    return pvar


def minus(xvar, yvar):
    """returns xvar - yvar"""
    return xvar - yvar


def divide(xvar, yvar):
    """returns xvar / yvar"""
    return xvar / yvar


def multiply(*xvar):
    """returns the product of all elements in the list xvar"""
    return times(*xvar)


SBML_FUNCT_DICT = {
    "acos": acos, "arccos": arccos, "acosh": acosh, "arccosh": arccosh,
    "acot": acot, "arccot": arccot, "acoth": acoth, "arccoth": arccoth,
    "acsc": acsc, "arccsc": arccsc, "acsch": acsch, "arccsch": arccsch,
    "arcsec": arcsec, "asech": asech, "arcsech": arcsech, "asin": asin,
    "asinh": asinh, "arcsinh": arcsinh, "arcsin": arcsin, "atan": atan,
    "arctan": arctan, "atanh": atanh, "arctanh": arctanh, "ceil": ceil,
    "ceiling": ceiling, "cos": cos, "cosh": cosh, "cot": cot,
    "coth": coth, "csc": csc, "csch": csch, "factorial": factorial,
    "exp": exp, "floor": floor, "ln": ln, "log": log, "log10": log10,
    "piecewise": piecewise, "pow": pow, "power": power, "root": root,
    "sec": sec, "sech": sech, "sqr": sqr, "sqrt": sqrt, "sin": sin,
    "sinh": sinh, "tan": tan, "tanh": tanh, "And": And, "Not": Not,
    "Or": Or, "xor": xor, "eq": eq, "geq": geq, "gt": gt, "leq": leq,
    "lt": lt, "neq": neq, "plus": plus, "times": times, "minus": minus,
    "divide": divide, "multiply": multiply, "exponentiale":exp(1),
    "INF":np.inf, "NAN":np.nan, "inf":np.inf, "nan":np.nan
}
