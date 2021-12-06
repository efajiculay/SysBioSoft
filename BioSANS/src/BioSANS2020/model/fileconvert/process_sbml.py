# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from inspect import getargspec as inspect_getargspec
from sympy import Symbol as newSymbol
from sympy import solve, sympify
import libsbml as Mysbml

from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.math_functs.sbmlMath import acos, arccos, acosh, \
    arccosh, acot, arccot, acoth, arccoth, acsc, arccsc, acsch, \
    arccsch, arcsec, asech, arcsech, asin, asinh, arcsinh, arcsin, \
    atan, arctan, atanh, arctanh, ceil, ceiling, cos, cosh, cot, \
    coth, csc, csch, factorial, exp, floor, ln, log, log10, \
    piecewise, pow, power, root, sec, sech, sqr, sqrt, sin, sinh, \
    tan, tanh, And, Not, Or, xor, eq, geq, gt, leq, lt, neq, \
    plus, times, minus, divide, multiply


SBML_FUNCT_LIST = [
    acos, arccos, acosh, arccosh, acot, arccot, acoth, arccoth,
    acsc, arccsc, acsch, arccsch, arcsec, asech, arcsech, asin,
    asinh, arcsinh, arcsin, atan, arctan, atanh, arctanh, ceil,
    ceiling, cos, cosh, cot, coth, csc, csch, factorial, exp,
    floor, ln, log, log10, piecewise, pow, power, root, sec,
    sech, sqr, sqrt, sin, sinh, tan, tanh, And, Not, Or, xor,
    eq, geq, gt, leq, lt, neq, plus, times, minus, divide,
    multiply
]

OPERS_LIST = {"+", "-", "*", "/", "(", ")", ",", "=", ">", "<", ":"}
OPERS_LIST2 = {
    "abs", "acos", "arccos", "acosh", "arccosh", "acot", "arccot", "acoth",
    "arccoth", "acsc", "acsch", "arccsch", "arcsec", "true", "false", "True",
    "False", "asinh", "asech", "arcsech", "asin", "arcsin", "atan", "arctan",
    "atanh", "arctanh", "ceil", "ceiling", "cos", "cosh", "cot", "coth", "csc",
    "csch", "and", "not", "arcsinh", "factorial", "exp", "floor", "ln", "log",
    "log10", "piecewise", "pow", "power", "root", "sec", "sech", "sqr", "sqrt",
    "sin", "sinh", "tan", "None", "pi", "arccsc", "avogadro", "tanh", "And",
    "Not", "Or", "or", "xor", "eq", "geq", "gt", "leq", "lt", "neq", "plus",
    "times", "minus", "divide", "if", "else", "multiply", "lambda", "delay",
    "rateOf"
}


def get_exponent_sp(key, modk):
    exponent_in_formula = 1
    query = "pow(" + key + ","
    ind = modk.find(query)
    if ind >= 0:
        stvar = ind + len(query)
        icvar = 0
        exponent_in_formula = ""
        while modk[stvar + icvar] != ")":
            exponent_in_formula = exponent_in_formula + modk[stvar + icvar]
            icvar = icvar + 1
    exponent_in_formula = float(exponent_in_formula)

    query = key + "**"
    ind = modk.find(query)
    if ind >= 0:
        stvar = ind + len(query)
        icvar = 0
        exponent_in_formula = ""
        while modk[stvar + icvar] != ")":
            exponent_in_formula = exponent_in_formula + modk[stvar + icvar]
            icvar = icvar + 1

    query = "/" + key
    ind = modk.find(query)
    if ind >= 0:
        exponent_in_formula = -1

    exponent_in_formula = float(exponent_in_formula)
    return exponent_in_formula


def replace_crucial_funct(trep):
    return Add2spacessep(trep) \
        .replace("and", "And") \
        .replace("not", "Not") \
        .replace("or", "Or") \
        .replace("factOrial", "factorial") \
        .replace("floOr", "floor") \
        .replace("xOr", "xor") \
        .replace("true", "True") \
        .replace("false", "False") \
        .replace(" ", "") \
        .replace("+-", "+ -")


def extract_species(modk):
    # global OPERS_LIST, OPERS_LIST2
    here = " " + str(modk).replace("time-", "emit-") \
        .replace("time+", "emit+").replace("e+", " ").replace("E+", " ") \
        .replace("e-", " ").replace("E-", " ").replace("emit-", "time-") \
        .replace("emit+", "time+") + " "
    for xvar in OPERS_LIST:
        here = here.replace(xvar, "  ")
    here = " " + here + " "
    for xvar in OPERS_LIST2:
        here = here.replace(" " + xvar + " ", " ")
    here = here.split()
    sp_comp = set()
    for xvar in here:
        try:
            float(xvar)
        except BaseException:
            try:
                float(eval(xvar))
            except BaseException:
                sp_comp.add(xvar)
    return ",".join(sp_comp)


def extract_par_num(modk):
    # global OPERS_LIST, OPERS_LIST2
    here = " " + str(modk) + " "
    for xvar in OPERS_LIST:
        here = here.replace(xvar, "  ")

    here = " " + here + " "
    for xvar in OPERS_LIST2:
        here = here.replace(" " + xvar + " ", " ")
    here = here.split()
    sp_comp = set()
    for xvar in here:
        try:
            float(xvar)
            sp_comp.add(xvar)
        except BaseException:
            pass

    return sp_comp


def extract_var_func(ssv):
    # v = ""
    last_comma = 0
    oper = {"+", "-", "*", "/", "(", ")"}
    pvar = 0
    for xvar in ssv + ")":
        if xvar == ",":
            last_comma = pvar
        elif xvar not in oper:
            pass
        else:
            return [ssv[0:last_comma], ssv[last_comma + 1:]]
        pvar = pvar + 1
    return ["", ""]


def extract_function(pforms, rbig_params, compartments, functions,
                     functions_str):
    # global OPERS_LIST, OPERS_LIST2
    spcorm = pforms.replace("(", " ").replace(")", " ").replace(
        ",", " ").replace("*", " ").split()
    modk = ""
    ind2 = 0
    for sp_comp in spcorm:
        if sp_comp in rbig_params:
            pass
        elif sp_comp in compartments:
            pass
        elif sp_comp in functions:
            npar = len(inspect_getargspec(functions[sp_comp])[0])
            ind = spcorm.index(sp_comp)
            ind2 = pforms.index(sp_comp)
            ssv = []
            for yvar in spcorm[ind + 1:]:
                try:
                    float(yvar)
                except BaseException:
                    if yvar not in OPERS_LIST and yvar not in OPERS_LIST2 \
                            and yvar not in functions:
                        ssv.append(newSymbol(yvar))
            if len(ssv) == npar:
                try:
                    modk = str(functions[sp_comp](*ssv))
                except BaseException:
                    modk = functions_str[sp_comp]
                    modk = FunctRedefineVar(modk, ssv)
            else:
                modk = "lambda " \
                    + ",".join([str(s) for s in ssv]) + " : " + pforms[ind2:]
                try:
                    modk = str(eval(modk)(*ssv))
                except BaseException:
                    modk = str(eval(modk[0:-1])(*ssv))
            break
        else:
            pass
    if len(modk) > 2:
        modk = "1*(" + modk + ")"
    modk = pforms[:ind2] + modk
    return modk


def parSubstitution(modk, parameters, index=0):
    # global OPERS_LIST
    modk = str(modk)
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, " " + yvar + " ")
    modk = " " + modk + " "
    for yvar in parameters:
        if parameters[yvar][0] is not None:
            modk = modk.replace(" " + yvar + " ", str(parameters[yvar][index]))
    return modk.replace(" ", "")


def Add2spacessep(modk):
    modk = str(modk)
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, "  " + yvar + "  ")
    return modk


def SpSubstitution(modk, parameters):
    # global OPERS_LIST
    modk = str(modk)
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, " " + yvar + " ")
    modk = " " + modk + " "
    for yvar in parameters:
        if parameters[yvar] is not None:
            modk = modk.replace(" " + yvar + " ", str(parameters[yvar]))
    return modk.replace(" ", "")


def varSubstitution(modk, rbig_params, parameters, compartments, RateRules):
    # global OPERS_LIST
    modk = str(modk)
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, "  " + yvar + "  ")
    modk = " " + modk + " "

    for yvar in rbig_params:
        if yvar not in RateRules:
            modk = modk.replace(" " + yvar + " ", "("
                                + str(rbig_params[yvar]) + ")")
    for yvar in parameters:
        if yvar not in RateRules:
            modk = modk.replace(" " + yvar + " ",
                                "(" + str(parameters[yvar][0]) + ")")
    for yvar in compartments:
        modk = modk.replace(" " + yvar + " ", "("
                            + str(compartments[yvar][0]) + ")")
    return modk.replace(" ", "")


def FunctRedefineVar(modk, ssv):
    spset = [b.strip() for b in modk.split(
        ":")[0].split("lambda")[1].split(",")]
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, " " + yvar + " ")
    modk = " " + modk + " "
    for s1 in range(len(ssv)):
        modk = modk.replace(" " + spset[s1] + " ", str(ssv[s1]))
    return modk.split(":")[1]


def process_sbml(file, molar=False, variables=None):
    # global OPERS_LIST, OPERS_LIST2

    fftopofile = open(file + ".topo", "w")

    SbmlUnits = {
        0: "ampere", 1: "avogadro", 2: "becquerel", 3: "candela", 5: "coulomb",
        6: "dimensionless", 7: "farad", 8: "gram", 9: "gray", 10: "henry",
        11: "hertz", 13: "joule", 14: "katal", 15: "kelvin", 18: "litre",
        19: "lumen", 20: "lux", 22: "metre", 23: "mole", 24: "newton",
        25: "ohm", 26: "pascal", 27: "radian", 28: "second", 29: "siemens",
        30: "sievert", 31: "steradian", 32: "tesla", 34: "watt", 35: "weber"
    }

    reader = Mysbml.SBMLReader()
    document = reader.readSBML(file)
    model = document.getModel()
    # errors = document.getNumErrors()

    time_var = None
    GCFactor = model.getConversionFactor()

    units = {}
    if model.getNumUnitDefinitions() > 0:
        for xvar in model.getListOfUnitDefinitions():
            c = xvar.getListOfUnits()
            units[xvar.getId()] = []
            for yvar in c:
                F = (yvar.getMultiplier() * 10 ** yvar.getScale()) \
                    ** yvar.getExponent()
                units[xvar.getId()].append(
                    [F, SbmlUnits[yvar.getKind()], yvar.getExponent()])

    compartments = {}
    constant_comp = {}
    nonConstant_comp = {}
    orig_size = {}
    for xvar in model.getListOfCompartments():
        if xvar.isSetSize():
            orig_size[xvar.getId()] = xvar.getSize()
            if not molar:
                compartments[xvar.getId()] = [xvar.getSize(), xvar.getUnits()]
            else:
                compartments[xvar.getId()] = [1, xvar.getUnits()]
        else:
            orig_size[xvar.getId()] = 1
            compartments[xvar.getId()] = [1, xvar.getUnits()]
        # if xvar.getConstant() == True:
        if xvar.getConstant():
            constant_comp[xvar.getId()] = compartments[xvar.getId()]
        else:
            nonConstant_comp[xvar.getId()] = compartments[xvar.getId()]

    reactions = {}
    for xvar in model.getListOfReactions():
        reactions[xvar.getId()] = xvar

    species = {}
    species_comp = {}
    constant_species = {}
    HasOnlySUnits = {}
    sp_wCFactor = {}
    for xvar in model.getListOfSpecies():
        species[xvar.getId()] = xvar
        species_comp[xvar.getId()] = xvar.getCompartment()
        HasOnlySUnits[xvar.getId()] = xvar.getHasOnlySubstanceUnits()
        # if xvar.getConstant() == True:
        if xvar.getConstant():
            constant_species[xvar.getId()] = True
        if xvar.getConversionFactor():
            sp_wCFactor[xvar.getId()] = xvar.getConversionFactor()
    # print(sp_wCFactor)

    SpInitialConc = {}
    for xvar in species:
        sp_comp = species[xvar]
        if sp_comp.isSetInitialConcentration():
            val = sp_comp.getInitialConcentration()
            SpInitialConc[xvar] = val
        elif sp_comp.isSetInitialAmount():
            val = sp_comp.getInitialAmount()
            SpInitialConc[xvar] = val
    # print(SpInitialConc)

    parameters = {}
    nonConstPar = set()
    constant_par = {}
    for xvar in model.getListOfParameters():
        if xvar.isSetValue():
            ssv = xvar.getValue()
        else:
            ssv = None
        parameters[xvar.getId()] = [ssv, xvar.getUnits()]
        # if xvar.getConstant() == False:
        if not xvar.getConstant():
            nonConstPar.add(xvar.getId())
        else:
            constant_par[xvar.getId()] = parameters[xvar.getId()]

    AssignRules = {}
    RateRules = {}
    algebr_rules = []
    for xvar in model.getListOfRules():
        if xvar.getMath().isAvogadro():
            ssv = str(6.02214179e+23)
        elif Mysbml.formulaToL3String(xvar.getMath()).find("avogadro") >= 0:
            ssv = replace_crucial_funct(Mysbml.formulaToL3String(
                xvar.getMath()).replace("avogadro", "(6.02214179e+23)"))
        else:
            ssv = replace_crucial_funct(Mysbml.formulaToString(xvar.getMath()))
        if xvar.isAssignment():
            ssv = parSubstitution(ssv, constant_par)
            AssignRules[xvar.getVariable()] = ssv
        elif xvar.isRate():
            RateRules[xvar.getVariable()] = ssv
        elif xvar.isAlgebraic():
            algebr_rules.append(ssv)

    functions = {}
    functions_str = {}
    fftopofile.write("Function_Definitions:\n")
    for xvar in model.getListOfFunctionDefinitions():
        ssv = replace_crucial_funct(Mysbml.formulaToString(xvar.getMath()))
        ssv = ssv.replace("lambda(", "")[:-1]
        ssv = extract_var_func(ssv)
        ssv = " : ".join(ssv)
        functions[xvar.getId()] = eval("lambda " + ssv)
        functions_str[xvar.getId()] = "lambda " + ssv
        globals2.execFunctions.append(xvar.getId() + " = lambda " + ssv)
        exec(xvar.getId() + " = lambda " + ssv, globals())
        fftopofile.write(xvar.getId() + " = lambda " + ssv + "\n")
        OPERS_LIST2.add(xvar.getId())

    InitialAssign = {}
    for xvar in model.getListOfInitialAssignments():
        sp_comp = xvar.getId()
        if xvar.getMath().isAvogadro():
            ssv = str(6.02214179e+23)
        elif Mysbml.formulaToL3String(xvar.getMath()).find("avogadro") >= 0:
            ssv = Mysbml.formulaToL3String(xvar.getMath()).replace(
                "avogadro", "(6.02214179e+23)")
            ssv = replace_crucial_funct(ssv.replace("lambda(", ""))
        else:
            ssv = replace_crucial_funct(Mysbml.formulaToString(
                xvar.getMath()).replace("lambda(", ""))
        if sp_comp in species:
            ssv = str(compartments[species_comp[sp_comp]][0]) \
                + "*(" + ssv + ")"
        try:
            ssv = eval(parSubstitution(ssv, parameters))
        except BaseException:
            ssv = parSubstitution(ssv, parameters)
        InitialAssign[sp_comp] = ssv

        if sp_comp in constant_par:
            constant_par[sp_comp][0] = ssv
            parameters[sp_comp][0] = ssv
        elif sp_comp in nonConstPar:
            parameters[sp_comp][0] = ssv
        elif sp_comp in compartments:
            orig_size[sp_comp] = ssv
            if not molar:
                compartments[sp_comp][0] = ssv
            else:
                compartments[sp_comp][0] = 1
        elif sp_comp in constant_comp:
            constant_comp[sp_comp][0] = ssv
        elif sp_comp in nonConstant_comp:
            nonConstant_comp[sp_comp][0] = ssv
    # print(InitialAssign)

    Events = {}
    EventAssign = {}
    Delays = {}
    IniValTrig = {}
    Priorities = []
    for xvar in model.getListOfEvents():
        uVFTT = xvar.getUseValuesFromTriggerTime()
        IniTV = xvar.getTrigger().getInitialValue()
        PersT = xvar.getTrigger().getPersistent()

        Prior = ""
        if xvar.isSetPriority():
            ssv = replace_crucial_funct(
                Mysbml.formulaToString(xvar.getPriority().getMath()))
            # ssv = SpSubstitution(ssv,SpInitialConc)
            Prior = ssv  # eval(ssv)
            Priorities.append(Prior)

        ssv = replace_crucial_funct(parSubstitution(
            Mysbml.formulaToString(xvar.getTrigger().getMath()), constant_par))
        try:
            for sp_comp in species:
                ssv = Add2spacessep(ssv) \
                    .replace(" " + sp_comp + " ", sp_comp + "/"
                             + str(compartments[species_comp[sp_comp]][0])) \
                    .replace(" ", "")
        except BaseException:
            pass

        mods = extract_species(ssv)
        for ttvar in mods.split(","):
            if ttvar not in species and ttvar not in parameters and ttvar \
                not in compartments and ttvar.strip() != "" \
                    and ttvar.find("delay") == -1:
                time_var = ttvar
                break

        Events[xvar.getId()] = ssv
        if not xvar.isSetDelay():
            for yvar in xvar.getListOfEventAssignments():
                IniValTrig[yvar.getVariable()] = IniTV
                try:
                    ssv = str(eval(parSubstitution(
                        Mysbml.formulaToString(yvar.getMath()), constant_par)))
                except BaseException:
                    ssv = parSubstitution(
                        Mysbml.formulaToString(yvar.getMath()), constant_par)
                    ssv = parSubstitution(ssv, constant_comp)
                ssv = replace_crucial_funct(ssv)
                mods = extract_species(ssv)
                for ttvar in mods.split(","):
                    if ttvar not in species and ttvar not in parameters \
                        and ttvar not in compartments and ttvar.strip() != "" \
                            and ttvar.find("delay") == -1:
                        time_var = ttvar
                        break

                sp_comp = yvar.getVariable()
                keyS = "status_" + sp_comp

                if sp_comp in species_comp:
                    ssv = str(compartments[species_comp[sp_comp]]
                              [0]) + "*(" + ssv + ")"

                if Prior == "":
                    estatus = " : True if " + Events[xvar.getId()] \
                        + " else False"
                    if PersT:
                        express = " : " + ssv + " if " + \
                            Events[xvar.getId()] + " and not " + keyS + \
                            " else " + yvar.getVariable()
                    else:
                        express = " : " + ssv + " if " + \
                            Events[xvar.getId()] + " and " + keyS + \
                            " else " + yvar.getVariable()
                else:
                    keyS = keyS + "_" + str(Priorities.index(Prior))
                    if PersT:
                        estatus = " : (True if " + \
                            Events[xvar.getId()] + " else False" + \
                            "," + Prior + ")"
                        express = " : (" + ssv + " if " \
                            + Events[xvar.getId()] \
                            + " and not " + keyS + " else " \
                            + yvar.getVariable() + "," + Prior + ")"
                    else:
                        estatus = " : (True if " + \
                            Events[xvar.getId()] + " else False" + \
                            "," + Prior + ",1)"
                        express = " : (" + ssv + " if " \
                            + Events[xvar.getId()] \
                            + " and " + keyS + " else " + yvar.getVariable() \
                            + "," + Prior + ",1)"

                estatus = "lambda " + extract_species(estatus) + estatus
                express = "lambda " + extract_species(express) + express
                if keyS in EventAssign:
                    EventAssign[keyS].append(estatus)
                else:
                    EventAssign[keyS] = [estatus]
                if yvar.getVariable() in EventAssign:
                    EventAssign[yvar.getVariable()].append(express)
                else:
                    EventAssign[yvar.getVariable()] = [express]
        else:
            Delays[xvar.getId()] = replace_crucial_funct(parSubstitution(
                Mysbml.formulaToString(xvar.getDelay().getMath()),
                constant_par))
            tdepDelay = False
            mods = extract_species(Delays[xvar.getId()])
            for ttvar in mods.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    tdepDelay = True
                    break

            for yvar in xvar.getListOfEventAssignments():
                IniValTrig[yvar.getVariable()] = IniTV

                sp_comp = yvar.getVariable()
                keyT = "timer_" + sp_comp
                keyS = "status_" + sp_comp
                keyS2 = "status2_" + sp_comp
                keyF = "finish_" + sp_comp
                keyC = "dtime_" + sp_comp
                keyD = "delay_" + sp_comp
                TheMath = replace_crucial_funct(parSubstitution(
                    Mysbml.formulaToString(yvar.getMath()), constant_par))

                if sp_comp in species_comp:
                    TheMath = str(compartments[species_comp[sp_comp]][0]) \
                        + "*(" + TheMath + ")"

                mods = extract_species(TheMath)
                for ttvar in mods.split(","):
                    if ttvar not in species and ttvar not in parameters \
                            and ttvar not in compartments \
                            and ttvar.strip() != "" \
                            and ttvar.find("delay") == -1:
                        time_var = ttvar
                        tdepDelay = True
                        break

                # if tdepDelay == True:
                if tdepDelay:
                    estatus2 = " : True if " + \
                        Events[xvar.getId()] + " else False"
                    estatus2 = "lambda " + extract_species(estatus2) + estatus2
                    delay_val = " : " + \
                        Delays[xvar.getId()] + " if " + Events[xvar.getId()] \
                        + " and " + keyF + " and not " + keyS2 \
                        + " else None"
                    delay_val = "lambda " + extract_species(delay_val) \
                        + delay_val
                    DelayCon = keyT + " >= " + keyD
                else:
                    DelayCon = keyT + " >= " + Delays[xvar.getId()]

                finishFr = " : 1 if " + DelayCon + " else 0 if " + \
                    Events[xvar.getId()] + " else None"
                finishFr = "lambda " + extract_species(finishFr) + finishFr

                estatus = " : True if " + DelayCon + " else False"
                estatus = "lambda " + extract_species(estatus) + estatus

                express1 = " : " + keyT + " if " + \
                    Events[xvar.getId()] + " or not " + keyF + " else 0"
                express1 = "lambda " + extract_species(express1) + express1

                CurrentT = " : " + TheMath + " if " + \
                    Events[xvar.getId()] + " and " + keyF + " else None"
                CurrentT = "lambda " + extract_species(CurrentT) + CurrentT

                if PersT:
                    keySCond = "not " + keyS
                    pass
                else:
                    estatus = " : True if " + Events[xvar.getId()] \
                        + " else False"
                    estatus = "lambda " + extract_species(estatus) + estatus
                    keySCond = keyS

                if uVFTT:
                    express2 = " : " + keyC + " if " + DelayCon + " and " \
                        + keySCond + " else " + sp_comp
                else:
                    express2 = " : " + TheMath + " if " + \
                        DelayCon + " and " + keySCond + " else " + sp_comp

                express2 = "lambda " + extract_species(express2) + express2

                if yvar.getVariable() in EventAssign:
                    EventAssign[keyT].append(express1)
                    EventAssign[sp_comp].append(express2)
                    EventAssign[keyS].append(estatus)
                    EventAssign[keyF].append(finishFr)
                    EventAssign[keyC].append(CurrentT)
                    # if tdepDelay == True:
                    if tdepDelay:
                        EventAssign[keyS2].append(estatus2)
                        EventAssign[keyD].append(delay_val)
                else:
                    EventAssign[keyT] = [express1]
                    EventAssign[sp_comp] = [express2]
                    EventAssign[keyS] = [estatus]
                    EventAssign[keyF] = [finishFr]
                    EventAssign[keyC] = [CurrentT]
                    # if tdepDelay == True:
                    if tdepDelay:
                        EventAssign[keyD] = [delay_val]
                        EventAssign[keyS2] = [estatus2]

    for xvar in AssignRules:
        ssv = AssignRules[xvar]
        ssv = parSubstitution(ssv, parameters)
        if xvar in compartments:
            compartments[xvar][0] = eval(ssv)
        elif xvar in parameters:
            try:
                parameters[xvar][0] = eval(ssv)
            except BaseException:
                pass

    big_none_added = False
    fftopofile.write("\n")
    if len(orig_size) == 1:
        for cx in orig_size:
            if cx in constant_comp:
                fftopofile.write("#Reactions Volume = " +
                                 str(orig_size[cx]) + "\n")
            else:
                fftopofile.write("#Reactions\n")
    else:
        fftopofile.write("#Reactions\n")

    rbig_params = {}
    use_species = set()
    fast_rxn_exp = {}
    fast_irrev = {}
    use_species_inrxn = set()
    StoichVar = {}
    StoichPar = {}
    CFinReact = []
    CFinProdu = []
    CFkeys = []

    for xvar in reactions:
        rbig_params = {}
        SFactor = 1
        react = ""
        # row = []
        if reactions[xvar].isSetKineticLaw():
            for kx in reactions[xvar].getKineticLaw().getListOfParameters():
                rbig_params[kx.getId()] = kx.getValue()

        withCF = False
        for rx in reactions[xvar].getListOfReactants():
            key = rx.getSpecies()
            if key in sp_wCFactor:
                CFkeys.append(sp_wCFactor[key])
                sp_wCFactor[key] = parameters[sp_wCFactor[key]][0]
                withCF = True
                CFinReact.append(key)
            elif GCFactor:
                CFkeys.append(GCFactor)
                sp_wCFactor[key] = parameters[GCFactor][0]
                withCF = True
                CFinReact.append(key)
            use_species.add(key)
            use_species_inrxn.add(key)
            if rx.isSetStoichiometryMath():
                ddvar = rx.getStoichiometryMath().getMath()
                ddvar = replace_crucial_funct(Mysbml.formulaToString(ddvar))
                ddvar = varSubstitution(
                    ddvar, rbig_params, parameters, constant_comp, RateRules)
                ddvar = eval(ddvar)
                react = react + str(ddvar) + " " + key + " + "
            elif rx.isSetId():
                idrx = rx.getId()
                if idrx in InitialAssign:
                    ddvar = varSubstitution(InitialAssign[idrx], rbig_params,
                                            parameters, constant_comp,
                                            RateRules)
                    ddvar = eval(ddvar)
                    react = react + str(ddvar) + " " + key + " + "
                elif idrx in RateRules:
                    ddvar = varSubstitution(RateRules[idrx], rbig_params,
                                            parameters, constant_comp,
                                            RateRules)
                    ddvar = eval(ddvar)
                    kxid = Mysbml.formulaToString(
                        reactions[xvar].getKineticLaw().getMath())
                    if kxid in nonConstPar:
                        StoichVar[kxid] = ddvar
                    else:
                        StoichPar[idrx] = ddvar
                        SFactor = idrx
                    react = react + str(1) + " " + key + " + "
            else:
                react = react + str(rx.getStoichiometry()) + " " + key + " + "

        if len(react) == 0:
            big_none_added = True
            react = "0 NONE => "
        else:
            react = react[0:-2] + " => "

        produ = ""
        SPactor = 1
        for rx in reactions[xvar].getListOfProducts():
            key = rx.getSpecies()
            if key in sp_wCFactor:
                CFkeys.append(sp_wCFactor[key])
                sp_wCFactor[key] = parameters[sp_wCFactor[key]][0]
                withCF = True
                CFinProdu.append(key)
            elif GCFactor:
                CFkeys.append(GCFactor)
                sp_wCFactor[key] = parameters[GCFactor][0]
                withCF = True
                CFinProdu.append(key)
            use_species.add(key)
            use_species_inrxn.add(key)
            if rx.isSetStoichiometryMath():
                ddvar = rx.getStoichiometryMath().getMath()
                ddvar = replace_crucial_funct(Mysbml.formulaToString(ddvar))
                ddvar = varSubstitution(
                    ddvar, rbig_params, parameters, constant_comp, RateRules)
                ddvar = eval(ddvar)
                produ = produ + str(ddvar) + " " + key + " + "
            elif rx.isSetId():
                idrx = rx.getId()
                if idrx in InitialAssign:
                    ddvar = varSubstitution(InitialAssign[idrx], rbig_params,
                                            parameters, constant_comp,
                                            RateRules)
                    ddvar = eval(ddvar)
                    produ = produ + str(ddvar) + " " + key + " + "
                elif idrx in RateRules:
                    ddvar = varSubstitution(RateRules[idrx], rbig_params,
                                            parameters, constant_comp,
                                            RateRules)
                    ddvar = eval(ddvar)
                    kxid = Mysbml.formulaToString(
                        reactions[xvar].getKineticLaw().getMath())
                    if kxid in nonConstPar:
                        StoichVar[kxid] = ddvar
                    else:
                        StoichPar[idrx] = ddvar
                        SPactor = idrx
                    produ = produ + str(1) + " " + key + " + "
                elif idrx in AssignRules:
                    ddvar = varSubstitution(AssignRules[idrx], rbig_params,
                                            parameters, constant_comp,
                                            RateRules)
                    try:
                        ddvar = eval(ddvar)
                        kxid = Mysbml.formulaToString(
                            reactions[xvar].getKineticLaw().getMath())
                        SPactor = idrx
                        if kxid in nonConstPar:
                            StoichVar[kxid] = ddvar
                        produ = produ + str(1) + " " + key + " + "
                    except BaseException:
                        kxid = Mysbml.formulaToString(
                            reactions[xvar].getKineticLaw().getMath())
                        if kxid in nonConstPar:
                            StoichVar[kxid] = ddvar
                        else:
                            StoichPar[idrx] = ddvar
                            SPactor = idrx
                        produ = produ + str(1) + " " + key + " + "

            else:
                produ = produ + str(rx.getStoichiometry()) + " " + key + " + "

        if len(produ) == 0:
            big_none_added = True
            produ = "0 NONE"
            react = react + produ
            react = react + ",1   :::::   "
        # elif withCF == True:
        elif withCF:
            react = react + "0 NONE" + ",1   :::::   "
            produ = "0 NONE => " + produ[0:-2] + ",1   :::::   "
        else:
            react = react + produ
            react = react[0:-2] + ",1   :::::   "

        kforts = reactions[xvar].getKineticLaw()
        pforms = replace_crucial_funct(
            Mysbml.formulaToString(kforts.getMath()))

        try:
            modk = extract_function(
                pforms, rbig_params, constant_comp, functions, functions_str)
        except BaseException:
            modk = ""

        if len(modk) == 0:
            modk = pforms
        else:
            pass

        modk = varSubstitution(
            modk, rbig_params, constant_par, constant_comp, RateRules)
        # if reactions[xvar].getReversible() == True:
        if reactions[xvar].getReversible():
            if modk.find("+ -") > -1:
                modk = modk.split("+ -")
            else:
                if modk[0] == "-":
                    modk = modk[1:].replace("-", "+-").split("+")
                    modk[0] = "-" + modk[0]
                else:
                    modk = modk.replace("-", "+-").split("+")

            factor = 1
            mods = extract_species(modk[0])
            for key in mods.split(","):
                if key in species:
                    if HasOnlySUnits[key]:
                        pass
                    else:
                        exponent_in_formula = get_exponent_sp(key, modk[0])
                        if species_comp[key] in constant_comp:
                            factor = factor * \
                                (1 / compartments[species_comp[key]]
                                 [0])**exponent_in_formula
                        else:
                            factor = factor * \
                                (1 / newSymbol(species_comp[key])
                                 )**exponent_in_formula

            if SFactor != 1:
                factor = "((" + str(factor) + ")*(" + str(SFactor) + "))"

            pactor = 1
            mods = extract_species(modk[-1])
            for key in mods.split(","):
                if key in species:
                    if HasOnlySUnits[key]:
                        pass
                    else:
                        exponent_in_formula = get_exponent_sp(key, modk[-1])
                        if species_comp[key] in constant_comp:
                            pactor = pactor * \
                                (1 / compartments[species_comp[key]]
                                 [0])**exponent_in_formula
                        else:
                            pactor = pactor * \
                                (1 / newSymbol(species_comp[key])
                                 )**exponent_in_formula

            if SPactor != 1:
                pactor = "((" + str(pactor) + ")*(" + str(SPactor) + "))"

            if modk[0] == modk[-1]:
                if SFactor != 1:
                    modk = "(" + modk[0] + "*(" + str(factor) + ")"
                elif SPactor != 1:
                    modk = "(" + modk[0] + "*(" + str(pactor) + ")"
                else:
                    modk = modk[0]
            # elif reactions[xvar].getFast() == True:
            elif reactions[xvar].getFast():
                if modk[0][0] == "-" or modk[0][1] == "-":
                    modk = "1000*(" + modk[0] + "*(" + str(pactor) + \
                        ")+(" + str(factor) + ")*" + modk[-1] + ")"
                else:
                    modk = "1000*(" + modk[0] + "*(" + str(factor) + \
                        ")+(" + str(pactor) + ")*" + modk[-1] + ")"
                sphere = extract_species(modk).split(",")
                for ssv in sphere:
                    fast_rxn_exp[ssv] = modk
            else:
                if modk[0][0] == "-" or modk[0][1] == "-":
                    modk = "(" + modk[0] + "*(" + str(pactor) + \
                        ")+(" + str(factor) + ")*" + modk[-1] + ")"
                else:
                    modk = "(" + modk[0] + "*(" + str(factor) + \
                        ")+(" + str(pactor) + ")*" + modk[-1] + ")"

            if modk.count("(") > modk.count(")"):
                modk = modk + ")"

            mods = extract_species(modk)
            modk = "lambda " + mods + " : " + modk
            for ttvar in mods.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    break
        else:
            mods = extract_species(modk)
            # if reactions[xvar].getFast() == True:
            if reactions[xvar].getFast():
                sphere = mods.split(",")
                for ssv in sphere:
                    fast_irrev[ssv] = modk

                for ssv in reactions[xvar].getListOfProducts():
                    s = ssv.getSpecies()
                    vvar = 0
                    if species[s].isSetInitialConcentration():
                        vvar = species[s].getInitialConcentration()
                    elif species[s].isSetInitialAmount():
                        vvar = species[s].getInitialAmount()
                    modk2 = modk
                    for sp_comp in mods.split(","):
                        vv2 = 0
                        if species[sp_comp].isSetInitialConcentration():
                            vv2 = species[sp_comp].getInitialConcentration()
                        elif species[sp_comp].isSetInitialAmount():
                            vv2 = species[sp_comp].getInitialAmount()
                        modk2 = modk2.replace(sp_comp, str(vv2))
                    fast_irrev[ssv.getSpecies()] = modk2 + "+" + str(vvar) + \
                        "-" + str(newSymbol(s))

            factor = 1
            for key in mods.split(","):
                if key in species:
                    if HasOnlySUnits[key]:
                        pass
                    else:
                        exponent_in_formula = get_exponent_sp(key, modk)
                        if species_comp[key] in constant_comp:
                            factor = factor * \
                                (1 / compartments[species_comp[key]]
                                 [0])**exponent_in_formula
                        else:
                            factor = factor * \
                                (1 / newSymbol(species_comp[key])
                                 )**exponent_in_formula

            if SFactor != 1:
                factor = "((" + str(factor) + ")*(" + str(SFactor) + "))"
            if SPactor != 1:
                factor = "((" + str(factor) + ")*(" + str(SPactor) + "))"

            modk = "lambda " + mods + " :" + str(factor) + "*" + modk
            for ttvar in mods.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    break
        # if withCF == True:
        if withCF:
            moda = modk.split(":")
            for key in CFinReact:
                val = str(sp_wCFactor[key])
                moda[1] = "(" + moda[1] + ")*(" + val + ")"
            react = react + moda[0] + ":" + moda[1]
            fftopofile.write(react + "\n")
            moda = modk.split(":")
            for key in CFinProdu:
                val = str(sp_wCFactor[key])
                moda[1] = "(" + moda[1] + ")*(" + val + ")"
            produ = produ + moda[0] + ":" + moda[1]
            fftopofile.write(produ + "\n")
            big_none_added = True
        else:
            react = react + modk
            fftopofile.write(react + "\n")
        # print(fast_rxn_exp)
        # print(fast_irrev)

    for key in CFkeys:
        fftopofile.write(key + " => 0 NONE, 0" + "\n")

    done_comp = set()
    for key in sp_wCFactor:
        comp = species_comp[key]
        if comp not in done_comp and comp in constant_comp:
            done_comp.add(comp)
            fftopofile.write(comp + " => 0 NONE, 0" + "\n")

    for yvar in range(len(algebr_rules)):
        ssv = varSubstitution(
            algebr_rules[yvar], {}, constant_par, constant_comp, RateRules)
        spss = extract_species(ssv).split(",")
        try:
            ssv = eval("lambda " + ",".join(spss) + " : " + ssv)
            sf = [newSymbol(xvar) for xvar in spss]
            ssv = ssv(*sf)
        except BaseException:
            pass
        algebr_rules[yvar] = "Eq(" + str(ssv) + ",0)"

    modifierssp = {}
    for xvar in reactions:
        for rx in reactions[xvar].getListOfModifiers():
            use_species.add(rx.getSpecies())
            modifierssp[rx.getSpecies()] = newSymbol(rx.getSpecies())
            fftopofile.write("0 NONE => " + rx.getSpecies() + ", 0" + "\n")
            big_none_added = True

    algebr_rules = sympify(algebr_rules)

    for xvar in EventAssign:
        if xvar not in species and xvar not in parameters and xvar \
                not in compartments and xvar.find("status") == -1 \
                and xvar.find("finish") == -1 and xvar.find("dtime") == -1 \
                and xvar.find("delay") == -1:
            use_species.add(xvar)
            fftopofile.write("0 NONE => " + xvar + ", 1\n")
            big_none_added = True

    for xvar in RateRules:
        if xvar in parameters or xvar in compartments or xvar in species:
            use_species.add(xvar)
            modk = extract_function(
                RateRules[xvar], {}, compartments, functions, functions_str)

            mods = extract_species(RateRules[xvar])
            for ttvar in mods.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    break

            if modk == "":
                modk = RateRules[xvar]

            reversible = False
            if modk.find("+ -") > -1:
                modk = modk.split("+ -")
                for iv in range(1, len(modk)):
                    modk[iv] = "-" + modk[iv]
                reversible = True
            elif modk.find("+-") > -1:
                modk = modk.split("+-")
                for iv in range(1, len(modk)):
                    modk[iv] = "-" + modk[iv]
                reversible = True
            else:
                pass

            if reversible:
                Pactors = []
                for iv in range(len(modk)):
                    pactor = 1
                    mods = extract_species(modk[iv])
                    for key in mods.split(","):
                        if key in species:
                            if HasOnlySUnits[key]:
                                pass
                            else:
                                exponent_in_formula = get_exponent_sp(
                                    key, modk[iv])
                                if species_comp[key] in constant_comp:
                                    pactor = pactor * \
                                        (1 / compartments[species_comp[key]]
                                         [0])**exponent_in_formula
                                else:
                                    pactor = pactor * \
                                        (1 / newSymbol(species_comp[key])
                                         )**exponent_in_formula
                    Pactors.append(pactor)

            else:
                factor = 1
                mods = extract_species(modk)
                for key in mods.split(","):
                    if key in species:
                        if HasOnlySUnits[key]:
                            pass
                        else:
                            exponent_in_formula = get_exponent_sp(key, modk)
                            if species_comp[key] in constant_comp:
                                factor = factor * \
                                    (1 / compartments[species_comp[key]]
                                     [0])**exponent_in_formula
                            else:
                                factor = factor * \
                                    (1 / newSymbol(species_comp[key])
                                     )**exponent_in_formula

            try:
                float(modk)
                if xvar in species:
                    if HasOnlySUnits[xvar]:
                        pass
                    else:
                        comp = compartments[species_comp[xvar]][0]
                        modk = str(
                            eval(
                                modk +
                                "*" +
                                str(comp) +
                                "*" +
                                str(factor)))
                fftopofile.write("0 NONE" + " => " + xvar + ", " + modk + "\n")
            except BaseException:
                if not reversible:
                    modk = modk.replace(xvar, "$" + xvar)
                    modk = varSubstitution(
                        modk, {}, constant_par, constant_comp, RateRules)
                    modk = "lambda " + extract_species(modk) + " : " + modk
                    modk = modk.replace("$" + xvar, xvar)
                    modk = modk + "*" + str(factor)
                else:
                    modks = "("
                    for iv in range(len(modk) - 1):
                        modks = modks + modk[iv] + \
                            "*(" + str(Pactors[iv]) + ")" + "+"
                    modks = modks + \
                        modk[len(modk) - 1] + \
                        "*(" + str(Pactors[len(modk) - 1]) + ")" + ")"
                    modk = modks
                    modk = modk.replace(xvar, "$" + xvar)
                    modk = varSubstitution(
                        modk, {}, constant_par, constant_comp, RateRules)
                    modk = "lambda " + extract_species(modk) + " : " + modk
                    modk = modk.replace("$" + xvar, xvar)

                if xvar in species:
                    if HasOnlySUnits[xvar]:
                        pass
                    else:
                        comp = compartments[species_comp[xvar]][0]
                        modk = modk + "*" + str(comp)
                    fftopofile.write("0 NONE" + " => " + xvar +
                                     ",1   :::::   " + modk + "\n")
                else:
                    fftopofile.write("0 NONE" + " => " + xvar +
                                     ",1   :::::   " + modk + "\n")
        else:
            pass
            # modk = RateRules[xvar]
            # fftopofile.write("0 NONE"+" => "+ xvar +", "+modk+"\n")

        big_none_added = True

    for xvar in AssignRules:
        if xvar in species or xvar in parameters or xvar in compartments:
            use_species.add(xvar)
            fftopofile.write("0 NONE" + " => " + xvar + ",0" + "\n")
        big_none_added = True
        ssv = varSubstitution(
            AssignRules[xvar], rbig_params, parameters, compartments,
            RateRules)
        sss = extract_species(ssv)
        for ttvar in sss.split(","):
            if ttvar not in species and ttvar not in parameters and ttvar \
                not in compartments and ttvar.strip() != "" \
                    and ttvar.find("delay") == -1:
                time_var = ttvar
                break

    timevar = None
    for xvar in InitialAssign:
        ssv = varSubstitution(
            InitialAssign[xvar], rbig_params, parameters, compartments,
            RateRules)
        sss = extract_species(ssv)
        for ttvar in sss.split(","):
            if ttvar not in species and ttvar not in parameters and ttvar \
                not in compartments and ttvar.strip() != "" \
                    and ttvar.find("delay") == -1:
                time_var = ttvar
                timevar = 0
                break

    for xvar in nonConstant_comp:
        if xvar not in RateRules and xvar not in AssignRules and xvar \
                not in algebr_rules:
            fftopofile.write("0 NONE" + " => " + xvar + ",0" + "\n")
        big_none_added = True

    for xvar in nonConstPar:
        if xvar not in RateRules and xvar not in AssignRules:
            if xvar in StoichVar:
                ssv = StoichVar[xvar]
                try:
                    float(ssv)
                    fftopofile.write(
                        "0 NONE" + " => " + xvar + "," + str(ssv) + "\n")
                except BaseException:
                    ssv = varSubstitution(
                        ssv, rbig_params, constant_par, compartments,
                        RateRules)
                    fftopofile.write("0 NONE" + " => " + xvar
                                     + ",1 ::::: lambda "
                                     + extract_species(ssv) + " : " + ssv
                                     + "\n")
            else:
                fftopofile.write("0 NONE" + " => " + xvar + ",0" + "\n")
        big_none_added = True

    for xvar in StoichPar:
        ssv = StoichPar[xvar]
        if xvar not in RateRules and xvar not in AssignRules:
            try:
                float(ssv)
                fftopofile.write("0 NONE" + " => " + xvar + ","
                                 + str(ssv) + "\n")
            except BaseException:
                ssv = varSubstitution(
                    ssv, rbig_params, constant_par, compartments, RateRules)
                fftopofile.write(
                    "0 NONE" + " => " + xvar + ",1 ::::: lambda "
                    + extract_species(ssv) + " : " + ssv + "\n")

    if len(reactions) == 0:
        for xvar in species:
            sp_comp = species[xvar]
            if xvar in AssignRules:
                pass
            elif xvar in RateRules:
                pass
            else:
                use_species.add(xvar)
                modifierssp[xvar] = newSymbol(xvar)
                fftopofile.write(xvar + " => 0 NONE, 0.0\n")

        for xvar in parameters:
            if xvar in InitialAssign and xvar not in RateRules:
                fftopofile.write(xvar + " => 0 NONE, 0.0\n")

        if len(nonConstPar) == 0:
            for xvar in parameters:
                fftopofile.write(xvar + " => 0 NONE, 0.0\n")
        else:
            for xvar in constant_par:
                fftopofile.write(xvar + " => 0 NONE, 0.0\n")

        big_none_added = True

    for xvar in species:
        if xvar not in use_species:
            fftopofile.write(xvar + " => 0 NONE, 0.0\n")
            big_none_added = True

    if time_var is not None:
        if time_var in StoichPar:
            fftopofile.write("0 NONE => " + time_var + ", 0\n")
        elif timevar is None:
            fftopofile.write("0 NONE => " + time_var + ", 1\n")
        elif timevar == 0:
            fftopofile.write("0 NONE => " + time_var + ", 0\n")

    for key in constant_comp:
        if variables:
            if key in variables:
                fftopofile.write(key + " => 0 NONE, 0\n")
                big_none_added = True

    fftopofile.write("\n@Concentrations\n")

    for xvar in species:
        sp_comp = species[xvar]
        # if xvar in fast_rxn_exp and sp_comp.getBoundaryCondition() == False:
        if xvar in fast_rxn_exp and not sp_comp.getBoundaryCondition():
            ssv = extract_species(fast_rxn_exp[xvar])
            spsum = ""
            sbsum = ""
            val = 0
            bal = 0
            syms = []
            noBoundary = True
            for s in ssv.split(","):
                syms.append(newSymbol(s))
                if species[s].isSetInitialConcentration():
                    vvar = species[s].getInitialConcentration()
                elif species[s].isSetInitialAmount():
                    vvar = species[s].getInitialAmount()

                # if species[s].getBoundaryCondition() == False:
                if not species[s].getBoundaryCondition():
                    val = val + vvar
                    spsum = spsum + s + "+"
                else:
                    noBoundary = False
                    sbsum = s
                    bal = bal + vvar

            if noBoundary:
                spsum = spsum[0:-1]
                eqst = sympify(["Eq(" + fast_rxn_exp[xvar] + ",0)",
                                "Eq(" + spsum + "," + str(val) + ")"])
            else:
                eqst = sympify(["Eq(" + fast_rxn_exp[xvar] + ",0)",
                                "Eq(" + sbsum + "," + str(bal) + ")"])
            sol = solve(eqst, syms)
            if sol:
                val1 = sol[newSymbol(xvar)]
            else:
                val1 = 0
        # elif xvar in fast_irrev and sp_comp.getBoundaryCondition() == False:
        elif xvar in fast_irrev and not sp_comp.getBoundaryCondition():
            eqst = sympify("Eq(" + fast_irrev[xvar] + ",0)")
            ssz = extract_species(fast_irrev[xvar])
            sol = solve(eqst, ssz)
            if sol:
                val1 = sol[0]
            else:
                val1 = 0
        elif sp_comp.isSetInitialConcentration():
            val1 = sp_comp.getInitialConcentration()
        elif sp_comp.isSetInitialAmount():
            val1 = sp_comp.getInitialAmount()
        else:
            val1 = ""

        if xvar in InitialAssign or xvar in SpInitialConc and xvar \
                not in fast_rxn_exp and xvar not in fast_irrev:
            try:
                ssv = varSubstitution(
                    InitialAssign[xvar],
                    rbig_params,
                    parameters,
                    compartments,
                    RateRules)
                for ini in SpInitialConc:
                    ssv = ssv.replace(ini, str(SpInitialConc[ini]))
            except BaseException:
                ssv = str(SpInitialConc[xvar])
            try:
                val1 = eval(ssv)
            except BaseException:
                val1 = 0
                val2 = ssv
        if xvar in AssignRules:
            ssv = varSubstitution(
                AssignRules[xvar], rbig_params, parameters, compartments,
                RateRules)
            sss = extract_species(ssv)
            for ttvar in sss.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    break
            try:
                val2 = eval(ssv)
                if float(val2):
                    if val1 == "":
                        val1 = val2
            except BaseException:
                val2 = ssv
                if val1 == "":
                    val1 = 0

        if xvar not in EventAssign:
            cmod = ", "
            if xvar in modifierssp or xvar not in use_species and xvar \
                    not in constant_species:
                ssphere = newSymbol(xvar)
                ddvar = solve(algebr_rules, ssphere)
                if ddvar:
                    cmod = cmod + "lambda " + \
                        extract_species(str(ddvar[ssphere])) \
                        + " : " + str(ddvar[ssphere])
            if cmod == ", ":
                cmod = ""
            # if sp_comp.getBoundaryCondition() == False and xvar \
            if not sp_comp.getBoundaryCondition() and xvar \
                    not in AssignRules:
                if species_comp[xvar] in nonConstant_comp and molar and xvar \
                        not in use_species_inrxn:
                    fftopofile.write(
                        sp_comp.getId() + " , " + str(val1) + cmod
                        + ", lambda " + species_comp[xvar] + " :"
                        + str(val1) + "/" + species_comp[xvar] + "\n")
                else:
                    # if str(val1) in val2:
                    # fftopofile.write(
                    #     sp_comp.getId() + " , " + str(val1) + cmod
                    #     + ", lambda " + sss + ":" + val2 + "\n");
                    # else:
                    fftopofile.write(
                        sp_comp.getId() + " , " + str(val1) + cmod + "\n")
            elif xvar in AssignRules:
                new_funct = extract_function(
                    str(val2), rbig_params, compartments, functions,
                    functions_str)
                if new_funct == "":
                    new_funct = val2
                fftopofile.write(sp_comp.getId() + " ," + str(val1) +
                                 ", lambda " +
                                 extract_species(new_funct) + ": " +
                                 str(new_funct) + str(cmod) + "\n")
            # elif sp_comp.getBoundaryCondition() == True \
            #     and xvar in RateRules:
            elif sp_comp.getBoundaryCondition() and xvar in RateRules:
                fftopofile.write(sp_comp.getId() + " , " + str(val1)
                                 + cmod + "\n")
            # elif sp_comp.getBoundaryCondition() == True and xvar \
            elif sp_comp.getBoundaryCondition() and xvar \
                not in RateRules and xvar not in AssignRules \
                    and cmod.strip() != "" and xvar not in EventAssign:
                fftopofile.write(sp_comp.getId() + " , " + str(val1)
                                 + cmod + "\n")
            else:
                fftopofile.write(sp_comp.getId() + " , " + str(val1) +
                                 ", lambda " + sp_comp.getId() + ": "
                                 + str(val1) + cmod + "\n")
        else:
            for ig in range(len(EventAssign[xvar])):
                cmod = ", "
                cmod = cmod + EventAssign[xvar][ig]
                # if sp_comp.getBoundaryCondition() == False:
                if not sp_comp.getBoundaryCondition():
                    fftopofile.write(
                        sp_comp.getId() + " , " + str(val1) + cmod + "\n")
                else:
                    ssv = cmod[1:].split(":")
                    lam_pa = ssv[0].replace("lambda", "").strip()
                    cond = ssv[1].split("if")[1].split("and")[0].strip()
                    con_pa = ssv[1].split("if")[0].strip()
                    fftopofile.write(sp_comp.getId() + " , " + str(val1)
                                     + ", lambda " + lam_pa + ": None if not "
                                     + cond + " else " + con_pa + "\n")

    if time_var is not None and time_var not in StoichPar:
        fftopofile.write(time_var + ", 0, lambda " + time_var +
                         " : round(" + time_var + ",10)\n")
        big_none_added = True

    for xvar in RateRules:
        if xvar in parameters and xvar not in EventAssign:
            fftopofile.write(xvar + "," + str(parameters[xvar][0]) + "\n")
            big_none_added = True
        elif xvar in compartments:
            fftopofile.write(xvar + " , " + str(compartments[xvar][0]) + "\n")
            big_none_added = True
        else:
            pass
            # fftopofile.write(xvar+" , 1\n")
            # big_none_added = True

    writtenEvent = {}
    for xvar in EventAssign:
        if xvar not in species:
            for ih in range(len(EventAssign[xvar])):
                if xvar in parameters:
                    ssv = str(parameters[xvar][0])
                elif xvar in nonConstant_comp:
                    ssv = str(nonConstant_comp[xvar][0])
                else:
                    if xvar.find("finish") > -1 or xvar.find("delay") > -1:
                        ssv = "1"
                    elif xvar.find("status") > -1:
                        # if IniValTrig[xvar.split("_")[1]] == True:
                        if IniValTrig[xvar.split("_")[1]]:
                            ssv = "1"
                        # elif IniValTrig[xvar.split("_")[1]] == False:
                        elif not IniValTrig[xvar.split("_")[1]]:
                            ssv = "0"
                    else:
                        ssv = "0"
                Eev = xvar + "," + ssv + "," + EventAssign[xvar][ih] + "\n"
                if Eev not in writtenEvent:
                    writtenEvent[Eev] = True
                    fftopofile.write(Eev)

    for xvar in AssignRules:
        if xvar not in species:
            try:
                float(AssignRules[xvar])
                fftopofile.write(xvar + ", " + AssignRules[xvar] + "\n")
            except BaseException:
                ssv = varSubstitution(AssignRules[xvar], rbig_params,
                                      constant_par, compartments, RateRules)
                if xvar in StoichPar:
                    fftopofile.write(xvar + ", 1, lambda " +
                                     extract_species(ssv) + " : " + ssv + "\n")
                else:
                    fftopofile.write(xvar + ", 0, lambda " +
                                     extract_species(ssv) + " : " + ssv + "\n")

    for xvar in InitialAssign:
        if xvar not in species and xvar not in RateRules:
            try:
                float(AssignRules[xvar])
                fftopofile.write(xvar + ", " + InitialAssign[xvar] + "\n")
            except BaseException:
                ssv = varSubstitution(InitialAssign[xvar], rbig_params,
                                      parameters, compartments, RateRules)
                if not molar:
                    fftopofile.write(xvar + ", 0, lambda " +
                                     extract_species(ssv) + " : " + ssv + "\n")

    for xvar in nonConstPar:
        if xvar not in RateRules and xvar not in AssignRules:
            # and xvar not in InitialAssign:
            ssphere = newSymbol(xvar)
            ddvar = solve(algebr_rules, ssphere)
            ssv = str(parameters[xvar][0] if parameters[xvar][0]
                      is not None else 0)
            try:
                fftopofile.write(xvar + "," + ssv + "," + "lambda " +
                                 extract_species(str(ddvar[ssphere]))
                                 + " : " + str(ddvar[ssphere]) + "\n")
            except BaseException:
                try:
                    float(ssv)
                    fftopofile.write(xvar + "," + ssv + "\n")
                except BaseException:
                    pass
        big_none_added = True

    for xvar in nonConstant_comp:
        if xvar in InitialAssign and xvar not in RateRules and xvar \
                not in AssignRules:
            fftopofile.write(
                xvar + "," + str(eval(str(nonConstant_comp[xvar][0])))
                + "\n" "\n")
            big_none_added = True
        elif xvar not in RateRules and xvar not in AssignRules \
                and xvar not in EventAssign:
            ssphere = newSymbol(xvar)
            ddvar = solve(algebr_rules, ssphere)
            fftopofile.write(xvar + ","
                             + str(eval(str(nonConstant_comp[xvar][0])))
                             + "," + "lambda " +
                             extract_species(str(ddvar[ssphere]))
                             + " : " + str(ddvar[ssphere]) + "\n" "\n")
            big_none_added = True

    if len(reactions) == 0:
        if len(nonConstPar) == 0:
            for xvar in parameters:
                try:
                    fftopofile.write(xvar + "," + str(parameters[xvar][0])
                                     + ", lambda :"
                                     + str(eval(parameters[xvar][0])) + "\n")
                except BaseException:
                    fftopofile.write(xvar + "," + "0" + ", lambda :" +
                                     str(parameters[xvar][0]) + "\n")
        else:
            for xvar in constant_par:
                try:
                    fftopofile.write(xvar + "," + str(parameters[xvar][0])
                                     + ", lambda :"
                                     + str(eval(parameters[xvar][0])) + "\n")
                except BaseException:
                    fftopofile.write(xvar + "," + "0" + ", lambda :" +
                                     str(parameters[xvar][0]) + "\n")

    for key in CFkeys:
        fftopofile.write(key + ", " + str(parameters[key][0]) + "\n")

    for key in StoichPar:
        if key not in AssignRules and key not in InitialAssign \
                and key not in RateRules:
            ssv = StoichPar[key]
            fftopofile.write(key + ", " + str(1) + "\n")

    done_comp = set()
    for key in sp_wCFactor:
        comp = species_comp[key]
        if comp not in done_comp and comp in constant_comp:
            done_comp.add(comp)
            fftopofile.write(comp + ", " + str(constant_comp[comp][0]) + "\n")

    for key in constant_comp:
        if variables:
            if key in variables:
                fftopofile.write(
                    key + ", " + str(constant_comp[key][0]) + "\n")

    if big_none_added:
        fftopofile.write("NONE , 1.0\n")

    fftopofile.close()
    return file + ".topo"
