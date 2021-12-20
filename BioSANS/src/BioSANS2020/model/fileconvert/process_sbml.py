"""

                        The process_sbml module

This module contains functions that facilitates the conversion of sbml
file format into BioSANS topology format. Currently, this module highly
relies on string manipulations, sympify, eval, and exec python functions

The following are the list of function in this module

1. get_exponent_sp
2. replace_crucial_funct
3. extract_species
4. extract_par_num
5. extract_var_func
6. extract_function
7. par_substitution
8. add2_spaces_sep
9. sp_substitution
10. var_substitution
11. funct_redefine_var
12. get_sbml_units
13. get_compartment_details
14. get_species_details
15. get_initial_conc
16. get_param_DETAILS
17. get_rule_details
18. process_sbml

"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from inspect import getargspec as inspect_getargspec
from sympy import Symbol as newSymbol
from sympy import solve, sympify
import libsbml as Mysbml

from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT


# "abs" : abs,
# "true" : true, "false" : false, "True" : True, "False" : False,
# "and" : and, "not" : not, "None" : None, "pi" : pi,
# "avogadro" : avogadro,
# "or" : or, "if" : if, "else" : else, "lambda" : lambda,
# "delay" : delay, "rateOf" : rateOf

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


def eval_exp(xvar):
    return eval(xvar, globals(), SBML_FUNCT_DICT)


def get_exponent_sp(key, modk):
    """This function extract the exponent of component or species key
    from a given propensity expression modk.
    Args:
        key : species/component string
        modk : propensity expression converted to string
               example "2*0.5*A*B"
    Returns:
        exponent_in_formula : the power of species key
    Example:
        modk = (1.0)*(1.0)*S1
        key = S1
        exponent_in_formula = 1
    """
    exponent_in_formula = 1
    # Checking if pow is present in the string
    query = "pow(" + key + ","
    ind = modk.find(query)
    if ind >= 0:
        # if found, grab the exponent after "," but before ")"
        stvar = ind + len(query)
        icvar = 0
        exponent_in_formula = ""
        while modk[stvar + icvar] != ")":
            exponent_in_formula = exponent_in_formula + modk[stvar + icvar]
            icvar = icvar + 1
    exponent_in_formula = float(exponent_in_formula)

    # Checking if ** is present in the string
    query = key + "**"
    ind = modk.find(query)
    if ind >= 0:
        # if found, grab the exponent after "**" but before ")"
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
    """this function converts some problematic operator/logical/function
    and etc. to its equivalent in SBML.
    Args:
        trep : result of SBML formulaToString(xvar.getMath())
               xvar is an object in SBML getListOfFunctionDefinitions()
    """
    return add2_spaces_sep(trep) \
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
    """This function extracts the components or species from a given ex-
    pression and returns a comma concatenated string of them.
    Args:
        modk : propensity expression converted to string
               example "2*0.5*A*B"
    Returns:
        ",".join(sp_comp) : a string of comma concatenated components
                            example "A,B"
    """
    # global OPERS_LIST, OPERS_LIST2
    # manipulation of string to avoid conflicts with scientific notation
    here = " " + str(modk).replace("time-", "emit-") \
        .replace("time+", "emit+").replace("e+", " ").replace("E+", " ") \
        .replace("e-", " ").replace("E-", " ").replace("emit-", "time-") \
        .replace("emit+", "time+") + " "
    for xvar in OPERS_LIST:
        here = here.replace(xvar, "  ")
    here = " " + here + " "
    for xvar in OPERS_LIST2:
        here = here.replace(" " + xvar + " ", " ")

    # grabbing the individual components
    here = here.split()
    sp_comp = set()
    for xvar in here:
        try:
            float(xvar)
        except:
            try:
                float(eval_exp(xvar))
            except:
                sp_comp.add(xvar)
    return ",".join(sp_comp)


def extract_par_num(modk):
    """This function is currently not used but it can extract the
    numeric values from a given expression
    Example:
        print(extract_par_num("(5.0)*(6.0)*S1+S2**3")) will give
        {'5.0', '6.0', '3'}
    """
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
        except:
            pass

    return sp_comp


def extract_var_func(ssv):
    """This function extract the variables from an expression and gives
    a 2D list containing a comma concatenated string of variable and a
    string of the expression.
    Args:
        ssv : comma concatenated string of variables and expression.
              example :  "x,y,x*y"
    Returns:
        ["comma concatenated variables","expression"] :
              example : ['x,y', 'x*y']
    """
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
    """This function converts a user defined function or function call
    in the SBML file to actual expression.
    Args:
        pforms : from formulaToString(xvar.getMath()) where xvar can be
                an object from getListOfRules() in libsbml.
                example : compartment*multiply(k1,S1)
        rbig_params : dictionary of parameter name : parameter value
                      from kinetic rate law expression in sbml
        compartments : dictionary of compartments name : [size, units]
        functions : dictionary of function name : function definition
        functions_str : the string equivalent of functions above
                        example : {'multiply': 'lambda x,y : x*y'}
    Returns:
        modk : expanded form of pforms where every function definition
               is evaluated
               example : compartment*1*(S1*k1)
    """
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
                except:
                    if yvar not in OPERS_LIST and yvar not in OPERS_LIST2 \
                            and yvar not in functions:
                        ssv.append(newSymbol(yvar))
            if len(ssv) == npar:
                try:
                    modk = str(functions[sp_comp](*ssv))
                except:
                    modk = functions_str[sp_comp]
                    modk = funct_redefine_var(modk, ssv)
            else:
                modk = "lambda " \
                    + ",".join([str(ssi) for ssi in ssv]) + " : " \
                    + pforms[ind2:]
                try:
                    modk = str(eval_exp(modk)(*ssv))
                except:
                    modk = str(eval_exp(modk[0:-1])(*ssv))
            break
        else:
            pass
    if len(modk) > 2:
        modk = "1*(" + modk + ")"
    modk = pforms[:ind2] + modk
    return modk


def par_substitution(modk, parameters, index=0):
    """This function substitute parameters to modk expression by first
    padding the operators to ensure parameters can be distinguished from
    the expression.
    Args:
        modk : math expression string
        parameters : dictionary of parameters name : [value, unit]
        index : 0 or 1. 0 will give value and 1 will give unit
    Returns:
        modk : expression where parameters are now numeric
    """
    # global OPERS_LIST
    modk = str(modk)
    for yvar in OPERS_LIST:
        # adding space to every operator
        modk = modk.replace(yvar, " " + yvar + " ")
    # padding the expression
    modk = " " + modk + " "
    for yvar in parameters:
        # substitution of numeric values to the parameter
        if parameters[yvar][0] is not None:
            modk = modk.replace(" " + yvar + " ", str(parameters[yvar][index]))
    # return without the padding
    return modk.replace(" ", "")


def add2_spaces_sep(modk):
    """This function adds paddings on both sides of the operators and
    return the expression modk with the operators padded.
    Args:
        modk : string math expression
    Returns:
        modk : string math expression with the operators padded
    """
    modk = str(modk)
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, "  " + yvar + "  ")
    return modk


def sp_substitution(modk, parameters):
    """ This is a redundant function to par_substitution. This function
    substitute parameters to modk expression by first padding the
    operators to ensure parameters can be distinguished fromthe
    expression.
    Args:
        modk : math expression string
        parameters : dictionary of parameters name : [value, unit]
                     Intended for constant parameters
    Returns:
        modk : expression where parameters are now numeric
    """
    # global OPERS_LIST
    modk = str(modk)
    for yvar in OPERS_LIST:
        # padding the operator
        modk = modk.replace(yvar, " " + yvar + " ")
    # padding the entire expression
    modk = " " + modk + " "
    for yvar in parameters:
        if parameters[yvar] is not None:
            # substitution of parameters
            modk = modk.replace(" " + yvar + " ", str(parameters[yvar]))
    return modk.replace(" ", "")


def var_substitution(modk, rbig_params, parameters, compartments, rate_rules):
    """This function substitute variables to modk expression by first
    padding the operators to ensure variables can be distinguished from
    the expression.
    Args:
        modk : math expression string
        rbig_params : dictionary of parameter name : parameter value
                      from kinetic rate law expression in sbml
        parameters : dictionary of parameters name : [value, unit]
        compartments : dictionary of compartments name : [size, units]
        rate_rules : dictionary of variable : variable modifier function
    Returns:
        modk : expression where variables are substituted
    """
    # global OPERS_LIST
    modk = str(modk)
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, "  " + yvar + "  ")
    modk = " " + modk + " "

    for yvar in rbig_params:
        if yvar not in rate_rules:
            modk = modk.replace(" " + yvar + " ", "("
                                + str(rbig_params[yvar]) + ")")
    for yvar in parameters:
        if yvar not in rate_rules:
            modk = modk.replace(" " + yvar + " ",
                                "(" + str(parameters[yvar][0]) + ")")
    for yvar in compartments:
        modk = modk.replace(" " + yvar + " ", "("
                            + str(compartments[yvar][0]) + ")")
    return modk.replace(" ", "")


def funct_redefine_var(modk, ssv):
    """This function transform the variables in a lambda function to the
    actual variables needed.
    Args:
        modk : string math lambda expression with abstract variables
        ssv : sympy symbols
    Returns:
        string definition expression with actual varaibles
    """
    # spset is a string list of components or species
    spset = [b.strip() for b in modk.split(
        ":")[0].split("lambda")[1].split(",")]
    # padding operators in modk
    for yvar in OPERS_LIST:
        modk = modk.replace(yvar, " " + yvar + " ")
    modk = " " + modk + " "
    for ss1, _ in enumerate(ssv):
        # substitution of sympy variables which is defined from the file
        modk = modk.replace(" " + spset[ss1] + " ", str(ssv[ss1]))
    return modk.split(":")[1]


def get_sbml_units(model, sbml_units):
    """Returns units dictionary from sbml file. The model is defined by
    the following python syntax;
        model = libsbml.SBMLReader().readSBML(sbml_file).getModel()
    Args:
        model : model definition like shown above
        sbml_units : a dictionary of sbml unit equivalence
     Returns :
        units_sbml : dictionary of sbml units associated to the file
    """
    units_sbml = {}
    if model.getNumUnitDefinitions() > 0:
        for xvar in model.getListOfUnitDefinitions():
            cvar = xvar.getListOfUnits()
            units_sbml[xvar.getId()] = []
            for yvar in cvar:
                fvar = (yvar.getMultiplier() * 10 ** yvar.getScale()) \
                    ** yvar.getExponent()
                units_sbml[xvar.getId()].append(
                    [fvar, sbml_units[yvar.getKind()], yvar.getExponent()])
    return units_sbml


def get_compartments_details(model, molar):
    """Returns compartment dictionary from sbml file. The model is
    defined by the following python syntax;
        model = libsbml.SBMLReader().readSBML(sbml_file).getModel()
    Args:
        model : model definition like shown above
        molar : True if conc. is in molar else False
    Returns:
        compartments : dictionary of all compartments with the following
                       example format { 'c1' : [size, units]}
        constant_comp : dictionary of constant compartments
                        example format { 'c1' : [size, units]}
        non_constant_comp : dictionary of non constant compartments
                            example format { 'c1' : [size, units]}
        orig_size : dictionary of original size of compartment
                    example format { 'c1' : size}
    """
    compartments = {}
    constant_comp = {}
    non_constant_comp = {}
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
            non_constant_comp[xvar.getId()] = compartments[xvar.getId()]
    return compartments, constant_comp, non_constant_comp, orig_size


def get_species_details(model):
    """Returns species dictionary from sbml file. The model is
    defined by the following python syntax;
        model = libsbml.SBMLReader().readSBML(sbml_file).getModel()
    Args:
        model : model definition like shown above
    Returns:
        species : dictionary of all species with the following
                       example format { id : species_label }
        species_comp : dictionary of species compartments
                       example format { species_label : compartment }
        constant_species : dictionary of constant species
                           example format { species_label : True }
        has_only_sunits : dictionary of species with only substance as
                          example format { id : True }
        sp_wcfactor : dictionary of convertion factor
                      example format { id : conversion factor }
    """
    species = {}
    species_comp = {}
    constant_species = {}
    has_only_sunits = {}
    sp_wcfactor = {}
    for xvar in model.getListOfSpecies():
        species[xvar.getId()] = xvar
        species_comp[xvar.getId()] = xvar.getCompartment()
        has_only_sunits[xvar.getId()] = xvar.getHasOnlySubstanceUnits()
        # if xvar.getConstant() == True:
        if xvar.getConstant():
            constant_species[xvar.getId()] = True
        if xvar.getConversionFactor():
            sp_wcfactor[xvar.getId()] = xvar.getConversionFactor()
    # print(sp_wcfactor)
    return species, species_comp, constant_species, has_only_sunits, \
        sp_wcfactor


def get_initial_conc(species):
    """Returns a dictionary of initial concentration from species"""
    sp_initial_conc = {}
    for xvar in species:
        sp_comp = species[xvar]
        if sp_comp.isSetInitialConcentration():
            val = sp_comp.getInitialConcentration()
            sp_initial_conc[xvar] = val
        elif sp_comp.isSetInitialAmount():
            val = sp_comp.getInitialAmount()
            sp_initial_conc[xvar] = val
    # print(sp_initial_conc)
    return sp_initial_conc


def get_param_details(model):
    """Returns a dictionary of parameter name : [value, unit]"""
    parameters = {}
    non_const_par = set()
    constant_par = {}
    for xvar in model.getListOfParameters():
        if xvar.isSetValue():
            ssv = xvar.getValue()
        else:
            ssv = None
        parameters[xvar.getId()] = [ssv, xvar.getUnits()]
        # if xvar.getConstant() == False:
        if not xvar.getConstant():
            non_const_par.add(xvar.getId())
        else:
            constant_par[xvar.getId()] = parameters[xvar.getId()]
    return parameters, non_const_par, constant_par


def get_rule_details(model, constant_par):
    """This function returns the assignment rules, rate rules, and
    algebraic rules from sbml model given the disctionary of constant
    parameters. The constants are already substituted in the return"""
    assign_rules = {}
    rate_rules = {}
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
            ssv = par_substitution(ssv, constant_par)
            assign_rules[xvar.getVariable()] = ssv
        elif xvar.isRate():
            rate_rules[xvar.getVariable()] = ssv
        elif xvar.isAlgebraic():
            algebr_rules.append(ssv)
    return assign_rules, rate_rules, algebr_rules


def process_sbml(file, molar=False, variables=None):
    """This function read the sbml file and extracts the necessary most
    of the details to construct the corresponding topology file of the
    system described by the SBML tags.
    Args:
        file : this is the sbml file
        molar : True or False, tells if the systems unit molar or not
        variables : None or constant compartment to show in plot
    """
    # global OPERS_LIST, OPERS_LIST2

    fftopofile = open(file + ".topo", "w")

    sbml_units = {
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
    gcfactor = model.getConversionFactor()

    units_sbml = get_sbml_units(model, sbml_units)
    compartments, constant_comp, non_constant_comp, orig_size \
        = get_compartments_details(model, molar)
    reactions = {xvar.getId(): xvar for xvar in model.getListOfReactions()}
    species, species_comp, constant_species, has_only_sunits, sp_wcfactor \
        = get_species_details(model)
    sp_initial_conc = get_initial_conc(species)
    parameters, non_const_par, constant_par = get_param_details(model)
    assign_rules, rate_rules, algebr_rules \
        = get_rule_details(model, constant_par)

    functions = {}
    functions_str = {}
    fftopofile.write("Function_Definitions:\n")
    for xvar in model.getListOfFunctionDefinitions():
        ssv = replace_crucial_funct(Mysbml.formulaToString(xvar.getMath()))
        ssv = ssv.replace("lambda(", "")[:-1]
        ssv = extract_var_func(ssv)
        ssv = " : ".join(ssv)
        functions[xvar.getId()] = eval_exp("lambda " + ssv)
        functions_str[xvar.getId()] = "lambda " + ssv
        globals2.EXEC_FUNCTIONS.append(xvar.getId() + " = lambda " + ssv)
        exec(xvar.getId() + " = lambda " + ssv, globals(), SBML_FUNCT_DICT)
        fftopofile.write(xvar.getId() + " = lambda " + ssv + "\n")
        OPERS_LIST2.add(xvar.getId())

    initial_assign = {}
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
            ssv = eval_exp(par_substitution(ssv, parameters))
        except:
            ssv = par_substitution(ssv, parameters)
        initial_assign[sp_comp] = ssv

        if sp_comp in constant_par:
            constant_par[sp_comp][0] = ssv
            parameters[sp_comp][0] = ssv
        elif sp_comp in non_const_par:
            parameters[sp_comp][0] = ssv
        elif sp_comp in compartments:
            orig_size[sp_comp] = ssv
            if not molar:
                compartments[sp_comp][0] = ssv
            else:
                compartments[sp_comp][0] = 1
        elif sp_comp in constant_comp:
            constant_comp[sp_comp][0] = ssv
        elif sp_comp in non_constant_comp:
            non_constant_comp[sp_comp][0] = ssv
    # print(initial_assign)

    events_var = {}
    event_assign = {}
    delays_var = {}
    ini_val_trig = {}
    priorities_var = []
    for xvar in model.getListOfEvents():
        uv_fft = xvar.getUseValuesFromTriggerTime()
        ini_tv = xvar.getTrigger().getInitialValue()
        persi_t = xvar.getTrigger().getPersistent()

        prior_var = ""
        if xvar.isSetPriority():
            ssv = replace_crucial_funct(
                Mysbml.formulaToString(xvar.getPriority().getMath()))
            # ssv = sp_substitution(ssv,sp_initial_conc)
            prior_var = ssv  # eval_exp(ssv)
            priorities_var.append(prior_var)

        ssv = replace_crucial_funct(par_substitution(
            Mysbml.formulaToString(xvar.getTrigger().getMath()), constant_par))
        try:
            for sp_comp in species:
                ssv = add2_spaces_sep(ssv) \
                    .replace(" " + sp_comp + " ", sp_comp + "/"
                             + str(compartments[species_comp[sp_comp]][0])) \
                    .replace(" ", "")
        except:
            pass

        mods = extract_species(ssv)
        for ttvar in mods.split(","):
            if ttvar not in species and ttvar not in parameters and ttvar \
                not in compartments and ttvar.strip() != "" \
                    and ttvar.find("delay") == -1:
                time_var = ttvar
                break

        events_var[xvar.getId()] = ssv
        if not xvar.isSetDelay():
            for yvar in xvar.getListOfEventAssignments():
                ini_val_trig[yvar.getVariable()] = ini_tv
                try:
                    ssv = str(eval_exp(par_substitution(
                        Mysbml.formulaToString(yvar.getMath()), constant_par)))
                except:
                    ssv = par_substitution(
                        Mysbml.formulaToString(yvar.getMath()), constant_par)
                    ssv = par_substitution(ssv, constant_comp)
                ssv = replace_crucial_funct(ssv)
                mods = extract_species(ssv)
                for ttvar in mods.split(","):
                    if ttvar not in species and ttvar not in parameters \
                        and ttvar not in compartments and ttvar.strip() != "" \
                            and ttvar.find("delay") == -1:
                        time_var = ttvar
                        break

                sp_comp = yvar.getVariable()
                key_s = "status_" + sp_comp

                if sp_comp in species_comp:
                    ssv = str(compartments[species_comp[sp_comp]]
                              [0]) + "*(" + ssv + ")"

                if prior_var == "":
                    estatus = " : True if " + events_var[xvar.getId()] \
                        + " else False"
                    if persi_t:
                        express = " : " + ssv + " if " + \
                            events_var[xvar.getId()] + " and not " + key_s + \
                            " else " + yvar.getVariable()
                    else:
                        express = " : " + ssv + " if " + \
                            events_var[xvar.getId()] + " and " + key_s + \
                            " else " + yvar.getVariable()
                else:
                    key_s = key_s + "_" + str(priorities_var.index(prior_var))
                    if persi_t:
                        estatus = " : (True if " + \
                            events_var[xvar.getId()] + " else False" + \
                            "," + prior_var + ")"
                        express = " : (" + ssv + " if " \
                            + events_var[xvar.getId()] \
                            + " and not " + key_s + " else " \
                            + yvar.getVariable() + "," + prior_var + ")"
                    else:
                        estatus = " : (True if " + \
                            events_var[xvar.getId()] + " else False" + \
                            "," + prior_var + ",1)"
                        express = " : (" + ssv + " if " \
                            + events_var[xvar.getId()] \
                            + " and " + key_s + " else " + yvar.getVariable() \
                            + "," + prior_var + ",1)"

                estatus = "lambda " + extract_species(estatus) + estatus
                express = "lambda " + extract_species(express) + express
                if key_s in event_assign:
                    event_assign[key_s].append(estatus)
                else:
                    event_assign[key_s] = [estatus]
                if yvar.getVariable() in event_assign:
                    event_assign[yvar.getVariable()].append(express)
                else:
                    event_assign[yvar.getVariable()] = [express]
        else:
            delays_var[xvar.getId()] = replace_crucial_funct(par_substitution(
                Mysbml.formulaToString(xvar.getDelay().getMath()),
                constant_par))
            tdep_delay = False
            mods = extract_species(delays_var[xvar.getId()])
            for ttvar in mods.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    tdep_delay = True
                    break

            for yvar in xvar.getListOfEventAssignments():
                ini_val_trig[yvar.getVariable()] = ini_tv

                sp_comp = yvar.getVariable()
                key_t = "timer_" + sp_comp
                key_s = "status_" + sp_comp
                key_s2 = "status2_" + sp_comp
                key_f = "finish_" + sp_comp
                key_c = "dtime_" + sp_comp
                key_d = "delay_" + sp_comp
                the_math = replace_crucial_funct(par_substitution(
                    Mysbml.formulaToString(yvar.getMath()), constant_par))

                if sp_comp in species_comp:
                    the_math = str(compartments[species_comp[sp_comp]][0]) \
                        + "*(" + the_math + ")"

                mods = extract_species(the_math)
                for ttvar in mods.split(","):
                    if ttvar not in species and ttvar not in parameters \
                            and ttvar not in compartments \
                            and ttvar.strip() != "" \
                            and ttvar.find("delay") == -1:
                        time_var = ttvar
                        tdep_delay = True
                        break

                # if tdep_delay == True:
                if tdep_delay:
                    estatus2 = " : True if " + \
                        events_var[xvar.getId()] + " else False"
                    estatus2 = "lambda " + extract_species(estatus2) + estatus2
                    delay_val = " : " + \
                        delays_var[xvar.getId()] + " if " \
                        + events_var[xvar.getId()] \
                        + " and " + key_f + " and not " + key_s2 \
                        + " else None"
                    delay_val = "lambda " + extract_species(delay_val) \
                        + delay_val
                    delay_con = key_t + " >= " + key_d
                else:
                    delay_con = key_t + " >= " + delays_var[xvar.getId()]

                finish_fr = " : 1 if " + delay_con + " else 0 if " + \
                    events_var[xvar.getId()] + " else None"
                finish_fr = "lambda " + extract_species(finish_fr) + finish_fr

                estatus = " : True if " + delay_con + " else False"
                estatus = "lambda " + extract_species(estatus) + estatus

                express1 = " : " + key_t + " if " + \
                    events_var[xvar.getId()] + " or not " + key_f + " else 0"
                express1 = "lambda " + extract_species(express1) + express1

                current_t = " : " + the_math + " if " + \
                    events_var[xvar.getId()] + " and " + key_f + " else None"
                current_t = "lambda " + extract_species(current_t) + current_t

                if persi_t:
                    key_scond = "not " + key_s
                    # pass
                else:
                    estatus = " : True if " + events_var[xvar.getId()] \
                        + " else False"
                    estatus = "lambda " + extract_species(estatus) + estatus
                    key_scond = key_s

                if uv_fft:
                    express2 = " : " + key_c + " if " + delay_con + " and " \
                        + key_scond + " else " + sp_comp
                else:
                    express2 = " : " + the_math + " if " + \
                        delay_con + " and " + key_scond + " else " + sp_comp

                express2 = "lambda " + extract_species(express2) + express2

                if yvar.getVariable() in event_assign:
                    event_assign[key_t].append(express1)
                    event_assign[sp_comp].append(express2)
                    event_assign[key_s].append(estatus)
                    event_assign[key_f].append(finish_fr)
                    event_assign[key_c].append(current_t)
                    # if tdep_delay == True:
                    if tdep_delay:
                        event_assign[key_s2].append(estatus2)
                        event_assign[key_d].append(delay_val)
                else:
                    event_assign[key_t] = [express1]
                    event_assign[sp_comp] = [express2]
                    event_assign[key_s] = [estatus]
                    event_assign[key_f] = [finish_fr]
                    event_assign[key_c] = [current_t]
                    # if tdep_delay == True:
                    if tdep_delay:
                        event_assign[key_d] = [delay_val]
                        event_assign[key_s2] = [estatus2]

    for xvar in assign_rules:
        ssv = assign_rules[xvar]
        ssv = par_substitution(ssv, parameters)
        if xvar in compartments:
            compartments[xvar][0] = eval_exp(ssv)
        elif xvar in parameters:
            try:
                parameters[xvar][0] = eval_exp(ssv)
            except:
                pass

    big_none_added = False
    fftopofile.write("\n")
    if len(orig_size) == 1:
        for ccx in orig_size:
            if ccx in constant_comp:
                fftopofile.write("#Reactions Volume = " +
                                 str(orig_size[ccx]) + "\n")
            else:
                fftopofile.write("#Reactions\n")
    else:
        fftopofile.write("#Reactions\n")

    rbig_params = {}
    use_species = set()
    fast_rxn_exp = {}
    fast_irrev = {}
    use_species_inrxn = set()
    stoich_var = {}
    stoich_par = {}
    cf_in_react = []
    cf_in_produ = []
    cf_keys = []

    for xvar in reactions:
        rbig_params = {}
        sfactor_var = 1
        react = ""
        # row = []
        if reactions[xvar].isSetKineticLaw():
            for kkx in reactions[xvar].getKineticLaw().getListOfParameters():
                rbig_params[kkx.getId()] = kkx.getValue()

        with_cf = False
        for rrx in reactions[xvar].getListOfReactants():
            key = rrx.getSpecies()
            if key in sp_wcfactor:
                cf_keys.append(sp_wcfactor[key])
                sp_wcfactor[key] = parameters[sp_wcfactor[key]][0]
                with_cf = True
                cf_in_react.append(key)
            elif gcfactor:
                cf_keys.append(gcfactor)
                sp_wcfactor[key] = parameters[gcfactor][0]
                with_cf = True
                cf_in_react.append(key)
            use_species.add(key)
            use_species_inrxn.add(key)
            if rrx.isSetStoichiometryMath():
                ddvar = rrx.getStoichiometryMath().getMath()
                ddvar = replace_crucial_funct(Mysbml.formulaToString(ddvar))
                ddvar = var_substitution(
                    ddvar, rbig_params, parameters, constant_comp, rate_rules)
                ddvar = eval_exp(ddvar)
                react = react + str(ddvar) + " " + key + " + "
            elif rrx.isSetId():
                idrx = rrx.getId()
                if idrx in initial_assign:
                    ddvar = var_substitution(initial_assign[idrx], rbig_params,
                                             parameters, constant_comp,
                                             rate_rules)
                    ddvar = eval_exp(ddvar)
                    react = react + str(ddvar) + " " + key + " + "
                elif idrx in rate_rules:
                    ddvar = var_substitution(rate_rules[idrx], rbig_params,
                                             parameters, constant_comp,
                                             rate_rules)
                    ddvar = eval_exp(ddvar)
                    kxid = Mysbml.formulaToString(
                        reactions[xvar].getKineticLaw().getMath())
                    if kxid in non_const_par:
                        stoich_var[kxid] = ddvar
                    else:
                        stoich_par[idrx] = ddvar
                        sfactor_var = idrx
                    react = react + str(1) + " " + key + " + "
            else:
                react = react + str(rrx.getStoichiometry()) + " " + key + " + "

        if len(react) == 0:
            big_none_added = True
            react = "0 NONE => "
        else:
            react = react[0:-2] + " => "

        produ = ""
        spactor_var = 1
        for rrx in reactions[xvar].getListOfProducts():
            key = rrx.getSpecies()
            if key in sp_wcfactor:
                cf_keys.append(sp_wcfactor[key])
                sp_wcfactor[key] = parameters[sp_wcfactor[key]][0]
                with_cf = True
                cf_in_produ.append(key)
            elif gcfactor:
                cf_keys.append(gcfactor)
                sp_wcfactor[key] = parameters[gcfactor][0]
                with_cf = True
                cf_in_produ.append(key)
            use_species.add(key)
            use_species_inrxn.add(key)
            if rrx.isSetStoichiometryMath():
                ddvar = rrx.getStoichiometryMath().getMath()
                ddvar = replace_crucial_funct(Mysbml.formulaToString(ddvar))
                ddvar = var_substitution(
                    ddvar, rbig_params, parameters, constant_comp, rate_rules)
                ddvar = eval_exp(ddvar)
                produ = produ + str(ddvar) + " " + key + " + "
            elif rrx.isSetId():
                idrx = rrx.getId()
                if idrx in initial_assign:
                    ddvar = var_substitution(initial_assign[idrx], rbig_params,
                                             parameters, constant_comp,
                                             rate_rules)
                    ddvar = eval_exp(ddvar)
                    produ = produ + str(ddvar) + " " + key + " + "
                elif idrx in rate_rules:
                    ddvar = var_substitution(rate_rules[idrx], rbig_params,
                                             parameters, constant_comp,
                                             rate_rules)
                    ddvar = eval_exp(ddvar)
                    kxid = Mysbml.formulaToString(
                        reactions[xvar].getKineticLaw().getMath())
                    if kxid in non_const_par:
                        stoich_var[kxid] = ddvar
                    else:
                        stoich_par[idrx] = ddvar
                        spactor_var = idrx
                    produ = produ + str(1) + " " + key + " + "
                elif idrx in assign_rules:
                    ddvar = var_substitution(assign_rules[idrx], rbig_params,
                                             parameters, constant_comp,
                                             rate_rules)
                    try:
                        ddvar = eval_exp(ddvar)
                        kxid = Mysbml.formulaToString(
                            reactions[xvar].getKineticLaw().getMath())
                        spactor_var = idrx
                        if kxid in non_const_par:
                            stoich_var[kxid] = ddvar
                        produ = produ + str(1) + " " + key + " + "
                    except:
                        kxid = Mysbml.formulaToString(
                            reactions[xvar].getKineticLaw().getMath())
                        if kxid in non_const_par:
                            stoich_var[kxid] = ddvar
                        else:
                            stoich_par[idrx] = ddvar
                            spactor_var = idrx
                        produ = produ + str(1) + " " + key + " + "

            else:
                produ = produ + str(rrx.getStoichiometry()) + " " + key + " + "

        if len(produ) == 0:
            big_none_added = True
            produ = "0 NONE"
            react = react + produ
            react = react + ",1   :::::   "
        # elif with_cf == True:
        elif with_cf:
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
        except:
            modk = ""

        if len(modk) == 0:
            modk = pforms
        else:
            pass

        modk = var_substitution(
            modk, rbig_params, constant_par, constant_comp, rate_rules)
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
                    if has_only_sunits[key]:
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

            if sfactor_var != 1:
                factor = "((" + str(factor) + ")*(" + str(sfactor_var) + "))"

            pactor = 1
            mods = extract_species(modk[-1])
            for key in mods.split(","):
                if key in species:
                    if has_only_sunits[key]:
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

            if spactor_var != 1:
                pactor = "((" + str(pactor) + ")*(" + str(spactor_var) + "))"

            if modk[0] == modk[-1]:
                if sfactor_var != 1:
                    modk = "(" + modk[0] + "*(" + str(factor) + ")"
                elif spactor_var != 1:
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
                    ssi = ssv.getSpecies()
                    vvar = 0
                    if species[ssi].isSetInitialConcentration():
                        vvar = species[ssi].getInitialConcentration()
                    elif species[ssi].isSetInitialAmount():
                        vvar = species[ssi].getInitialAmount()
                    modk2 = modk
                    for sp_comp in mods.split(","):
                        vv2 = 0
                        if species[sp_comp].isSetInitialConcentration():
                            vv2 = species[sp_comp].getInitialConcentration()
                        elif species[sp_comp].isSetInitialAmount():
                            vv2 = species[sp_comp].getInitialAmount()
                        modk2 = modk2.replace(sp_comp, str(vv2))
                    fast_irrev[ssv.getSpecies()] = modk2 + "+" + str(vvar) + \
                        "-" + str(newSymbol(ssi))

            factor = 1
            for key in mods.split(","):
                if key in species:
                    if has_only_sunits[key]:
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

            if sfactor_var != 1:
                factor = "((" + str(factor) + ")*(" + str(sfactor_var) + "))"
            if spactor_var != 1:
                factor = "((" + str(factor) + ")*(" + str(spactor_var) + "))"

            modk = "lambda " + mods + " :" + str(factor) + "*" + modk
            for ttvar in mods.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    break
        # if with_cf == True:
        if with_cf:
            moda = modk.split(":")
            for key in cf_in_react:
                val = str(sp_wcfactor[key])
                moda[1] = "(" + moda[1] + ")*(" + val + ")"
            react = react + moda[0] + ":" + moda[1]
            fftopofile.write(react + "\n")
            moda = modk.split(":")
            for key in cf_in_produ:
                val = str(sp_wcfactor[key])
                moda[1] = "(" + moda[1] + ")*(" + val + ")"
            produ = produ + moda[0] + ":" + moda[1]
            fftopofile.write(produ + "\n")
            big_none_added = True
        else:
            react = react + modk
            fftopofile.write(react + "\n")
        # print(fast_rxn_exp)
        # print(fast_irrev)

    for key in cf_keys:
        fftopofile.write(key + " => 0 NONE, 0" + "\n")

    done_comp = set()
    for key in sp_wcfactor:
        comp = species_comp[key]
        if comp not in done_comp and comp in constant_comp:
            done_comp.add(comp)
            fftopofile.write(comp + " => 0 NONE, 0" + "\n")

    for yvar in range(len(algebr_rules)):
        ssv = var_substitution(
            algebr_rules[yvar], {}, constant_par, constant_comp, rate_rules)
        spss = extract_species(ssv).split(",")
        try:
            ssv = eval_exp("lambda " + ",".join(spss) + " : " + ssv)
            ssi = [newSymbol(xvar) for xvar in spss]
            ssv = ssv(*ssi)
        except:
            pass
        algebr_rules[yvar] = "Eq(" + str(ssv) + ",0)"

    modifierssp = {}
    for xvar in reactions:
        for rrx in reactions[xvar].getListOfModifiers():
            use_species.add(rrx.getSpecies())
            modifierssp[rrx.getSpecies()] = newSymbol(rrx.getSpecies())
            fftopofile.write("0 NONE => " + rrx.getSpecies() + ", 0" + "\n")
            big_none_added = True

    algebr_rules = sympify(algebr_rules)

    for xvar in event_assign:
        if xvar not in species and xvar not in parameters and xvar \
                not in compartments and xvar.find("status") == -1 \
                and xvar.find("finish") == -1 and xvar.find("dtime") == -1 \
                and xvar.find("delay") == -1:
            use_species.add(xvar)
            fftopofile.write("0 NONE => " + xvar + ", 1\n")
            big_none_added = True

    for xvar in rate_rules:
        if xvar in parameters or xvar in compartments or xvar in species:
            use_species.add(xvar)
            modk = extract_function(
                rate_rules[xvar], {}, compartments, functions, functions_str)

            mods = extract_species(rate_rules[xvar])
            for ttvar in mods.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    break

            if modk == "":
                modk = rate_rules[xvar]

            reversible = False
            if modk.find("+ -") > -1:
                modk = modk.split("+ -")
                for ivi in range(1, len(modk)):
                    modk[ivi] = "-" + modk[ivi]
                reversible = True
            elif modk.find("+-") > -1:
                modk = modk.split("+-")
                for ivi in range(1, len(modk)):
                    modk[ivi] = "-" + modk[ivi]
                reversible = True
            else:
                pass

            if reversible:
                pactors_var = []
                for ivi, _ in enumerate(modk):
                    pactor = 1
                    mods = extract_species(modk[ivi])
                    for key in mods.split(","):
                        if key in species:
                            if has_only_sunits[key]:
                                pass
                            else:
                                exponent_in_formula = get_exponent_sp(
                                    key, modk[ivi])
                                if species_comp[key] in constant_comp:
                                    pactor = pactor * \
                                        (1 / compartments[species_comp[key]]
                                         [0])**exponent_in_formula
                                else:
                                    pactor = pactor * \
                                        (1 / newSymbol(species_comp[key])
                                         )**exponent_in_formula
                    pactors_var.append(pactor)

            else:
                factor = 1
                mods = extract_species(modk)
                for key in mods.split(","):
                    if key in species:
                        if has_only_sunits[key]:
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
                    if has_only_sunits[xvar]:
                        pass
                    else:
                        comp = compartments[species_comp[xvar]][0]
                        modk = str(
                            eval_exp(
                                modk +
                                "*" +
                                str(comp) +
                                "*" +
                                str(factor)))
                fftopofile.write("0 NONE" + " => " + xvar + ", " + modk + "\n")
            except:
                if not reversible:
                    modk = modk.replace(xvar, "$" + xvar)
                    modk = var_substitution(
                        modk, {}, constant_par, constant_comp, rate_rules)
                    modk = "lambda " + extract_species(modk) + " : " + modk
                    modk = modk.replace("$" + xvar, xvar)
                    modk = modk + "*" + str(factor)
                else:
                    modks = "("
                    for ivi in range(len(modk) - 1):
                        modks = modks + modk[ivi] + \
                            "*(" + str(pactors_var[ivi]) + ")" + "+"
                    modks = modks + \
                        modk[len(modk) - 1] + \
                        "*(" + str(pactors_var[len(modk) - 1]) + ")" + ")"
                    modk = modks
                    modk = modk.replace(xvar, "$" + xvar)
                    modk = var_substitution(
                        modk, {}, constant_par, constant_comp, rate_rules)
                    modk = "lambda " + extract_species(modk) + " : " + modk
                    modk = modk.replace("$" + xvar, xvar)

                if xvar in species:
                    if has_only_sunits[xvar]:
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
            # modk = rate_rules[xvar]
            # fftopofile.write("0 NONE"+" => "+ xvar +", "+modk+"\n")

        big_none_added = True

    for xvar in assign_rules:
        if xvar in species or xvar in parameters or xvar in compartments:
            use_species.add(xvar)
            fftopofile.write("0 NONE" + " => " + xvar + ",0" + "\n")
        big_none_added = True
        ssv = var_substitution(
            assign_rules[xvar], rbig_params, parameters, compartments,
            rate_rules)
        sss = extract_species(ssv)
        for ttvar in sss.split(","):
            if ttvar not in species and ttvar not in parameters and ttvar \
                not in compartments and ttvar.strip() != "" \
                    and ttvar.find("delay") == -1:
                time_var = ttvar
                break

    timevar = None
    for xvar in initial_assign:
        ssv = var_substitution(
            initial_assign[xvar], rbig_params, parameters, compartments,
            rate_rules)
        sss = extract_species(ssv)
        for ttvar in sss.split(","):
            if ttvar not in species and ttvar not in parameters and ttvar \
                not in compartments and ttvar.strip() != "" \
                    and ttvar.find("delay") == -1:
                time_var = ttvar
                timevar = 0
                break

    for xvar in non_constant_comp:
        if xvar not in rate_rules and xvar not in assign_rules and xvar \
                not in algebr_rules:
            fftopofile.write("0 NONE" + " => " + xvar + ",0" + "\n")
        big_none_added = True

    for xvar in non_const_par:
        if xvar not in rate_rules and xvar not in assign_rules:
            if xvar in stoich_var:
                ssv = stoich_var[xvar]
                try:
                    float(ssv)
                    fftopofile.write(
                        "0 NONE" + " => " + xvar + "," + str(ssv) + "\n")
                except:
                    ssv = var_substitution(
                        ssv, rbig_params, constant_par, compartments,
                        rate_rules)
                    fftopofile.write("0 NONE" + " => " + xvar
                                     + ",1 ::::: lambda "
                                     + extract_species(ssv) + " : " + ssv
                                     + "\n")
            else:
                fftopofile.write("0 NONE" + " => " + xvar + ",0" + "\n")
        big_none_added = True

    for xvar in stoich_par:
        ssv = stoich_par[xvar]
        if xvar not in rate_rules and xvar not in assign_rules:
            try:
                float(ssv)
                fftopofile.write("0 NONE" + " => " + xvar + ","
                                 + str(ssv) + "\n")
            except:
                ssv = var_substitution(
                    ssv, rbig_params, constant_par, compartments, rate_rules)
                fftopofile.write(
                    "0 NONE" + " => " + xvar + ",1 ::::: lambda "
                    + extract_species(ssv) + " : " + ssv + "\n")

    if len(reactions) == 0:
        for xvar in species:
            sp_comp = species[xvar]
            if xvar in assign_rules:
                pass
            elif xvar in rate_rules:
                pass
            else:
                use_species.add(xvar)
                modifierssp[xvar] = newSymbol(xvar)
                fftopofile.write(xvar + " => 0 NONE, 0.0\n")

        for xvar in parameters:
            if xvar in initial_assign and xvar not in rate_rules:
                fftopofile.write(xvar + " => 0 NONE, 0.0\n")

        if len(non_const_par) == 0:
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
        if time_var in stoich_par:
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
            no_boundary = True
            for ssi in ssv.split(","):
                syms.append(newSymbol(ssi))
                if species[ssi].isSetInitialConcentration():
                    vvar = species[ssi].getInitialConcentration()
                elif species[ssi].isSetInitialAmount():
                    vvar = species[ssi].getInitialAmount()

                # if species[ssi].getBoundaryCondition() == False:
                if not species[ssi].getBoundaryCondition():
                    val = val + vvar
                    spsum = spsum + ssi + "+"
                else:
                    no_boundary = False
                    sbsum = ssi
                    bal = bal + vvar

            if no_boundary:
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

        if xvar in initial_assign or xvar in sp_initial_conc and xvar \
                not in fast_rxn_exp and xvar not in fast_irrev:
            try:
                ssv = var_substitution(
                    initial_assign[xvar],
                    rbig_params,
                    parameters,
                    compartments,
                    rate_rules)
                for ini in sp_initial_conc:
                    ssv = ssv.replace(ini, str(sp_initial_conc[ini]))
            except:
                ssv = str(sp_initial_conc[xvar])
            try:
                val1 = eval_exp(ssv)
            except:
                val1 = 0
                val2 = ssv
        if xvar in assign_rules:
            ssv = var_substitution(
                assign_rules[xvar], rbig_params, parameters, compartments,
                rate_rules)
            sss = extract_species(ssv)
            for ttvar in sss.split(","):
                if ttvar not in species and ttvar not in parameters and ttvar \
                    not in compartments and ttvar.strip() != "" \
                        and ttvar.find("delay") == -1:
                    time_var = ttvar
                    break
            try:
                val2 = eval_exp(ssv)
                if float(val2):
                    if val1 == "":
                        val1 = val2
            except:
                val2 = ssv
                if val1 == "":
                    val1 = 0

        if xvar not in event_assign:
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
                    not in assign_rules:
                if species_comp[xvar] in non_constant_comp and molar and xvar \
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
            elif xvar in assign_rules:
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
            #     and xvar in rate_rules:
            elif sp_comp.getBoundaryCondition() and xvar in rate_rules:
                fftopofile.write(sp_comp.getId() + " , " + str(val1)
                                 + cmod + "\n")
            # elif sp_comp.getBoundaryCondition() == True and xvar \
            elif sp_comp.getBoundaryCondition() and xvar \
                not in rate_rules and xvar not in assign_rules \
                    and cmod.strip() != "" and xvar not in event_assign:
                fftopofile.write(sp_comp.getId() + " , " + str(val1)
                                 + cmod + "\n")
            else:
                fftopofile.write(sp_comp.getId() + " , " + str(val1) +
                                 ", lambda " + sp_comp.getId() + ": "
                                 + str(val1) + cmod + "\n")
        else:
            for igi in range(len(event_assign[xvar])):
                cmod = ", "
                cmod = cmod + event_assign[xvar][igi]
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

    if time_var is not None and time_var not in stoich_par:
        fftopofile.write(time_var + ", 0, lambda " + time_var +
                         " : round(" + time_var + ",10)\n")
        big_none_added = True

    for xvar in rate_rules:
        if xvar in parameters and xvar not in event_assign:
            fftopofile.write(xvar + "," + str(parameters[xvar][0]) + "\n")
            big_none_added = True
        elif xvar in compartments:
            fftopofile.write(xvar + " , " + str(compartments[xvar][0]) + "\n")
            big_none_added = True
        else:
            pass
            # fftopofile.write(xvar+" , 1\n")
            # big_none_added = True

    written_event = {}
    for xvar in event_assign:
        if xvar not in species:
            for ihi in range(len(event_assign[xvar])):
                if xvar in parameters:
                    ssv = str(parameters[xvar][0])
                elif xvar in non_constant_comp:
                    ssv = str(non_constant_comp[xvar][0])
                else:
                    if xvar.find("finish") > -1 or xvar.find("delay") > -1:
                        ssv = "1"
                    elif xvar.find("status") > -1:
                        # if ini_val_trig[xvar.split("_")[1]] == True:
                        if ini_val_trig[xvar.split("_")[1]]:
                            ssv = "1"
                        # elif ini_val_trig[xvar.split("_")[1]] == False:
                        elif not ini_val_trig[xvar.split("_")[1]]:
                            ssv = "0"
                    else:
                        ssv = "0"
                eev_var = xvar + "," + ssv + "," + event_assign[xvar][ihi] \
                    + "\n"
                if eev_var not in written_event:
                    written_event[eev_var] = True
                    fftopofile.write(eev_var)

    for xvar in assign_rules:
        if xvar not in species:
            try:
                float(assign_rules[xvar])
                fftopofile.write(xvar + ", " + assign_rules[xvar] + "\n")
            except:
                ssv = var_substitution(assign_rules[xvar], rbig_params,
                                       constant_par, compartments, rate_rules)
                if xvar in stoich_par:
                    fftopofile.write(xvar + ", 1, lambda " +
                                     extract_species(ssv) + " : " + ssv + "\n")
                else:
                    fftopofile.write(xvar + ", 0, lambda " +
                                     extract_species(ssv) + " : " + ssv + "\n")

    for xvar in initial_assign:
        if xvar not in species and xvar not in rate_rules:
            try:
                float(assign_rules[xvar])
                fftopofile.write(xvar + ", " + initial_assign[xvar] + "\n")
            except:
                ssv = var_substitution(initial_assign[xvar], rbig_params,
                                       parameters, compartments, rate_rules)
                if not molar:
                    fftopofile.write(xvar + ", 0, lambda " +
                                     extract_species(ssv) + " : " + ssv + "\n")

    for xvar in non_const_par:
        if xvar not in rate_rules and xvar not in assign_rules:
            # and xvar not in initial_assign:
            ssphere = newSymbol(xvar)
            ddvar = solve(algebr_rules, ssphere)
            ssv = str(parameters[xvar][0] if parameters[xvar][0]
                      is not None else 0)
            try:
                fftopofile.write(xvar + "," + ssv + "," + "lambda " +
                                 extract_species(str(ddvar[ssphere]))
                                 + " : " + str(ddvar[ssphere]) + "\n")
            except:
                try:
                    float(ssv)
                    fftopofile.write(xvar + "," + ssv + "\n")
                except:
                    pass
        big_none_added = True

    for xvar in non_constant_comp:
        if xvar in initial_assign and xvar not in rate_rules and xvar \
                not in assign_rules:
            fftopofile.write(
                xvar + "," + str(eval_exp(str(non_constant_comp[xvar][0])))
                + "\n" "\n")
            big_none_added = True
        elif xvar not in rate_rules and xvar not in assign_rules \
                and xvar not in event_assign:
            ssphere = newSymbol(xvar)
            ddvar = solve(algebr_rules, ssphere)
            fftopofile.write(xvar + ","
                             + str(eval_exp(str(non_constant_comp[xvar][0])))
                             + "," + "lambda " +
                             extract_species(str(ddvar[ssphere]))
                             + " : " + str(ddvar[ssphere]) + "\n" "\n")
            big_none_added = True

    if len(reactions) == 0:
        if len(non_const_par) == 0:
            for xvar in parameters:
                try:
                    fftopofile.write(xvar + "," + str(parameters[xvar][0])
                                     + ", lambda :"
                                     + str(eval_exp(parameters[xvar][0]))
                                     + "\n")
                except:
                    fftopofile.write(xvar + "," + "0" + ", lambda :" +
                                     str(parameters[xvar][0]) + "\n")
        else:
            for xvar in constant_par:
                try:
                    fftopofile.write(xvar + "," + str(parameters[xvar][0])
                                     + ", lambda :"
                                     + str(eval_exp(parameters[xvar][0]))
                                     + "\n")
                except:
                    fftopofile.write(xvar + "," + "0" + ", lambda :" +
                                     str(parameters[xvar][0]) + "\n")

    for key in cf_keys:
        fftopofile.write(key + ", " + str(parameters[key][0]) + "\n")

    for key in stoich_par:
        if key not in assign_rules and key not in initial_assign \
                and key not in rate_rules:
            ssv = stoich_par[key]
            fftopofile.write(key + ", " + str(1) + "\n")

    done_comp = set()
    for key in sp_wcfactor:
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
