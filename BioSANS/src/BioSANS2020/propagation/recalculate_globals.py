"""

             This module is the recalculate_globals module

The purpose of this module is to ensure that global variables holds True
each process for multitrajectory simulation.


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))


import random as rm
import re
from sympy import re as sympy_re
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT


reserve_events_words = {
    "t", "time", "status", "status2", "timer", "finish", "delay", "dtime"}


CONCP = None
RATESP = None
ORASP = None
SIP = None
ACTUALSP = None
CONCP2 = None


def eval2(to_eval):
    """Evaluate a string expression using the SBML_FUNCT_DICT dictionary

    Args:
        to_eval (str): expression

    Returns:
        multitype: evaluated result
    """
    return eval(to_eval, SBML_FUNCT_DICT)


def rate_of(xvar):
    """This is to supplement the SBML rateOf function. The return is the
    instantaneous rate at the current state of the system"""
    xvar = RATESP
    indx = SIP.index(xvar)
    if len(ORASP) < 2:
        return 0
    delt = ORASP[-1] - ORASP[-2]
    rate = (CONCP[-1][indx] - CONCP[-2][indx]) / delt
    return rate

SBML_FUNCT_DICT['rateOf'] = rate_of


def delay_part1(yvar, last):
    """Part of the delay function. See delay function for details
    """
    delt = ORASP[-1] - ORASP[-2]
    inc = int(yvar / delt)
    i = inc
    while abs(abs(last - ORASP[-i]) / yvar - 1.0) > 1.0e-2:
        if last - ORASP[-i] < yvar:
            i = i + 1
        else:
            i = i - 1
    return CONCP[-inc - 1][SIP.index(ACTUALSP)]


def delay_part2(last, delt, yvar):
    """Part of the delay function. See delay function for details
    """
    for yhvar in range(len(globals2.MODIFIED[ACTUALSP])):
        spvar, pfunc = globals2.MODIFIED[ACTUALSP][yhvar]
        spvar = [csp.strip() for csp in spvar]
        if "t" in spvar or "time" in spvar:
            try:
                return pfunc(last + delt - yvar)
            except:
                sprv = []
                for csp, _ in enumerate(spvar):
                    sprv.append(CONCP2[spvar[csp].strip()])
                return pfunc(*sprv)
    return None


def delay(xvar, yvar):
    """This is to supplement the SBML delay function"""
    if ORASP:
        last = ORASP[-1]
        # first = ORASP[0]
        if last - yvar >= 0:
            return delay_part1(yvar, last)
        try:
            # if True:
            if len(ORASP) > 1:
                delt = ORASP[-1] - ORASP[-2]
            else:
                delt = 0
            if ACTUALSP in globals2.MODIFIED:
                val = delay_part2(last, delt, yvar)
                if val:
                    return val
            else:
                return 0
        except BaseException:
            xvar = ACTUALSP
            return CONCP[0][SIP.index(xvar)]

    return CONCP[SIP.index(ACTUALSP)]


def none_to_list(xvar=None):
    """Transform a None variable into a list [] to avoid dangerous
    initial values."""
    if not xvar:
        return []
    return xvar


def apply_rules(conc, yconc, oras=None, spconc=None, slabels=None):
    """This function modifies the concentration base on the rules stated
    in BioSANS topology file.

    Args:
        conc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
        yconc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
        oras (float, optional): current time point. Defaults to None.
        spconc (list, optional): values of conc at the current time.
            Defaults to None.
        slabels (list, optional): names of components. Defaults to None.
    """    """"""
    global ORASP, CONCP, SIP, ACTUALSP, CONCP2, RATESP
    ORASP = none_to_list(oras)
    CONCP = none_to_list(spconc)
    CONCP2 = none_to_list(conc)
    SIP = none_to_list(slabels)

    tuples_pfunct = []
    for xvar in globals2.MODIFIED:
        for yvar in range(len(globals2.MODIFIED[xvar])):
            spvar, pfunc = globals2.MODIFIED[xvar][yvar]
            ACTUALSP = spvar[0].strip()
            RATESP = spvar[-1].strip()
            try:
                # if True:
                sprv = []
                for csp, _ in enumerate(spvar):
                    sprv.append(conc[spvar[csp].strip()])
                suby = pfunc(*sprv)
                try:
                    suby = sympy_re(suby.evalf())
                except BaseException:
                    pass
                if not isinstance(suby, tuple):
                    if suby is not None:
                        yconc[xvar] = suby
                else:
                    tuples_pfunct.append(
                        [suby[1], xvar, suby[0], spvar, pfunc, len(suby)])
            # """
            except:
                suby = pfunc()
                try:
                    suby = sympy_re(suby.evalf())
                except BaseException:
                    pass
                if not isinstance(suby, tuple):
                    if suby is not None:
                        yconc[xvar] = suby
                else:
                    tuples_pfunct.append(
                        [suby[1], xvar, suby[0], spvar, pfunc, len(suby)])
            # """

    current_key = None
    update_sp = []
    perst = True
    for val in tuples_pfunct:
        if val[5] == 3:
            perst = False

    while len(tuples_pfunct) > 0:
        tuples_pfunct.sort(key=lambda a: a[0])
        last = 0
        for ihvar, _ in enumerate(tuples_pfunct):
            if tuples_pfunct[last] == tuples_pfunct[ihvar]:
                pass
            else:
                bslice = tuples_pfunct[last:ihvar]
                rm.shuffle(bslice)
                tuples_pfunct[last:ihvar] = bslice
                last = ihvar

        row = tuples_pfunct.pop()
        if current_key is None:
            current_key = row[0]
        if row[2] is not None:
            xvar = row[1]
            yconc[xvar] = row[2]
            if current_key == row[0]:
                update_sp.append(xvar)
            if len(tuples_pfunct) == 0 or current_key != row[0]:
                if perst:
                    for x_o in update_sp:
                        if x_o.split("_")[0] not in reserve_events_words:
                            conc[x_o] = yconc[x_o]
                else:
                    for x_o in update_sp:
                        if row[5] == 3:
                            conc[x_o] = yconc[x_o]
                current_key = row[0]
                update_sp = [xvar]

                for tup in tuples_pfunct:
                    spvar, pfunc = tup[3:5]
                    sprv = []
                    for csp, _ in enumerate(spvar):
                        sprv.append(conc[spvar[csp].strip()])
                    suby = pfunc(*sprv)
                    tup[0] = suby[1]
                    tup[2] = suby[0]

    for xvar in globals2.MODIFIED:
        conc[xvar] = yconc[xvar]


def get_globals(rfile):
    """This function get the global variables from BioSANS topology file
     called rfile.

    Args:
        rfile (str): name of topology file
    """

    globals2.MODIFIED = {}
    globals2.PROP_MODIFIED = {}
    setattr(globals2, 'TCHECK', [])  # TCHECK

    try:
        with open(rfile, "r") as file:
            rows = []
            last = ""
            for row in file:
                if last == "Function_Definitions":
                    if row.strip() != "" and row[0] != "#":
                        exec(row.strip(), SBML_FUNCT_DICT)
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
                        if len(cvar) >= 3:
                            ccv = ",".join(cvar[2:])
                            try:
                                res = re.findall(r"\w{1,3}\(t,.+\)\s", ccv)
                                if res:
                                    try:
                                        rrv = float(str(res[0])
                                                    .split(",")[1]
                                                    .replace(")", ""))
                                    except:
                                        rrv = eval2(str(res[0])
                                                    .split(",")[1]
                                                    .replace(")", ""))
                                    globals2.TCHECK.append(rrv)
                            except BaseException:
                                pass
                            cc2 = ccv.split(":")[0].replace(
                                "lambda", "").split(",")

                            delay_pos = ccv.replace(
                                "delay_", "_yaled").find("delay")
                            if delay_pos >= 0:
                                dfirst = ""
                                for ibi in range(delay_pos + 6, len(ccv)):
                                    dvar = ccv[ibi]
                                    if dvar == ",":
                                        break
                                    dfirst = dfirst + dvar
                                dfirst = dfirst.strip()
                                cc2new = [
                                    dfirst] + [
                                        xvar for xvar in cc2
                                        if xvar.strip() != dfirst.strip()]
                                cc2 = cc2new
                                gcc = ccv.split(":")[1]
                                ccv = "lambda " + ",".join(cc2) + ":" + gcc

                            rate_pos = ccv.find("rateOf")
                            if rate_pos >= 0:
                                dfirst = ""
                                for ibi in range(rate_pos + 7, len(ccv)):
                                    dvar = ccv[ibi]
                                    if dvar == ")":
                                        break
                                    dfirst = dfirst + dvar
                                dfirst = dfirst.strip()
                                cc2new = [
                                    xvar for xvar
                                    in cc2 if xvar.strip() != dfirst.strip()
                                ] + [dfirst]
                                cc2 = cc2new
                                gcc = ccv.split(":")[1]
                                ccv = "lambda " + ",".join(cc2) + ":" + gcc

                            if cvar[0].strip() not in globals2.MODIFIED:
                                globals2.MODIFIED[cvar[0].strip()] = [
                                    [cc2, eval2(ccv)]]
                            else:
                                globals2.MODIFIED[cvar[0].strip()].append(
                                    [cc2, eval2(ccv)])
                elif row[0] == "#":
                    last = "#"
                    ggv = row.split(",")[1:]
                    try:
                        for xvar in ggv:
                            xxvar = xvar.split("=")
                            globals2.SETTINGS[xxvar[0].strip()
                                              ] = xxvar[1].strip()
                    except BaseException:
                        pass
                elif row[0] == "@":
                    last = "@"
                elif row.strip() == "Function_Definitions:":
                    last = "Function_Definitions"
            file.close()

        rxn = len(rows)
        for ihvar in range(rxn):
            col_row = rows[ihvar].split(":::::")
            row = col_row[0].strip().split(",")
            if len(col_row) > 1:
                krow = col_row[1].strip().split(":::")
                if len(krow) == 2:
                    cc2 = krow[0].split(":")[0].replace(
                        "lambda", "").split(",")
                    cc3 = krow[1].split(":")[0].replace(
                        "lambda", "").split(",")
                    globals2.PROP_MODIFIED["Prop_" + str(ihvar)] = \
                        [(cc2, eval2(krow[0])), (cc3, eval2(krow[1]))]
                else:
                    cc2 = krow[0].split(":")[0].replace(
                        "lambda", "").split(",")
                    globals2.PROP_MODIFIED["Prop_" + str(ihvar)] = \
                        [(cc2, eval2(krow[0]))]

        globals2.TCHECK.sort()
    except BaseException:
        pass
