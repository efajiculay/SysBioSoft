"""

             This module is the convtopotosbml module

The sole purpose of this module is to faciliate conversion of topology
files into SBML files. Currently, this module still have a lot of
problems and the converted file may not run properly.

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from sympy import Symbol, Matrix
from libsbml import SBMLDocument, UNIT_KIND_ITEM, UNIT_KIND_SECOND, \
    UNIT_KIND_MOLE, UNIT_KIND_LITER, parseL3Formula, writeSBMLToFile

# from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.propagation.propensity import propensity_vec, \
    propensity_vec_molar
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.propagation.recalculate_globals import get_globals


def topo_to_sbml(sp_comp, ks_vals, conc, r_react, r_prods, vbig, v_volms,
                 items=None, molar=False, rfile=None):
    """This function helps in the conversion of topology file to SBML.
    Consider a file containing the following reaction topology;

        #REACTIONS
        A <=> B, 0.5, 0.3
        3 B => 2 C, 0.04

        @CONCENTRATION
        A, 100
        B, 90
        C, 80

    Args:
        sp_comp : dictionary of components where the key is the
                  component and the value is a set of integers that tell
                  the rank or order of that reaction as it appears in
                  the file. For the above topology
                  sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}
        ks_vals : dictioary of parameters where the key is the line in
                  the file and the value is the list of parameter found
                  i.e. for the above topology;
                  ks_vals = {0: [0.5, 0.3], 1: [0.04]}
        conc : dictionary of initial concentration i.e.
               conc = {'A': 100.0, 'B': 90.0, 'C': 80}
        r_react : dictionary where the key is the line in the file and
                  the value is a dictionary of reactants where the key
                  is the component and the value is its stoichiometric
                  coefficient in that line i.e. for the above topology;
                  r_react = {0: {'A': 1}, 1: {'B': 3}}
        r_prods : dictionary where the key is the line in the file and
                  the value is a dictionary of products where the key
                  is the component and the value is its stoichiometric
                  coefficient in that line i.e. for the above topology;
                  r_prods = {0: {'B': 1}, 1: {'C': 2}}
        vbig : not used but pertains to stoichiometric matrix
        v_volms : volume of compartment
        items : three dimensional array of [canvas, scroll_x, scroll_y]
        molar : True if conc. is in molar otherwise False
        rfile : path of the topology file + name of topology file
    """
    get_globals(rfile)
    cs_comp = {}
    for xvar in sp_comp:
        cs_comp[xvar] = Symbol(xvar, real=True, negative=False)

    kcs_vals = []
    for i, _ in enumerate(ks_vals):
        row = []
        if len(ks_vals[i]) == 1:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
        else:
            key = 'kf' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
            key = 'kb' + str(i + 1)
            row.append(Symbol(key, real=True, positive=True))
        kcs_vals.append(row)

    document = SBMLDocument(3, 1)
    model = document.createModel()
    model.setTimeUnits("second")

    comp = model.createCompartment()
    comp.setName('BasicCompartment')
    comp.setId('Basic')
    comp.setVolume(v_volms)
    comp.setUnits('litre')
    comp.setConstant(True)
    comp.setSpatialDimensions(3)

    substance = model.createUnitDefinition()
    substance.setId('substance')
    unit = substance.createUnit()
    unit.setKind(UNIT_KIND_ITEM)

    if not molar:
        fprop = Matrix(propensity_vec(kcs_vals, cs_comp, r_react, r_prods))
        units = []
        item_per_second = model.createUnitDefinition()
        item_per_second.setId('molecules_per_second')
        units.append(item_per_second.createUnit())
        units[-1].setKind(UNIT_KIND_ITEM)
        units[-1].setExponent(1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        units.append(item_per_second.createUnit())
        units[-1].setKind(UNIT_KIND_SECOND)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)

        per_second = model.createUnitDefinition()
        per_second.setId('per_second')
        units.append(per_second.createUnit())
        units[-1].setKind(UNIT_KIND_SECOND)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)

        per_item_second = model.createUnitDefinition()
        per_item_second.setId('per_molecules_second')
        units.append(per_item_second.createUnit())
        units[-1].setKind(UNIT_KIND_ITEM)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        units.append(per_item_second.createUnit())
        units[-1].setKind(UNIT_KIND_SECOND)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)

        molecules = model.createUnitDefinition()
        molecules.setId('molecules')
        units.append(molecules.createUnit())
        units[-1].setKind(UNIT_KIND_ITEM)
        units[-1].setExponent(1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)

        model.setSubstanceUnits("molecules")
    else:
        fprop = Matrix(
            propensity_vec_molar(
                kcs_vals,
                cs_comp,
                r_react,
                r_prods))
        units = []
        item_per_second = model.createUnitDefinition()
        item_per_second.setId('molar_per_second')
        units.append(item_per_second.createUnit())
        units[-1].setKind(UNIT_KIND_MOLE)
        units[-1].setExponent(1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        units.append(item_per_second.createUnit())
        units[-1].setKind(UNIT_KIND_LITER)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        units.append(item_per_second.createUnit())
        units[-1].setKind(UNIT_KIND_SECOND)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)

        per_second = model.createUnitDefinition()
        per_second.setId('per_second')
        units.append(per_second.createUnit())
        units[-1].setKind(UNIT_KIND_SECOND)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)

        per_item_second = model.createUnitDefinition()
        per_item_second.setId('per_molar_second')
        units.append(per_item_second.createUnit())
        units[-1].setKind(UNIT_KIND_LITER)
        units[-1].setExponent(1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        units.append(per_item_second.createUnit())
        units[-1].setKind(UNIT_KIND_MOLE)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        units.append(per_item_second.createUnit())
        units[-1].setKind(UNIT_KIND_SECOND)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)

        molarity = model.createUnitDefinition()
        molarity.setId('moles_per_liter')
        units.append(molarity.createUnit())
        units[-1].setKind(UNIT_KIND_MOLE)
        units[-1].setExponent(1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        units.append(molarity.createUnit())
        units[-1].setKind(UNIT_KIND_LITER)
        units[-1].setExponent(-1)
        units[-1].setScale(0)
        units[-1].setMultiplier(1)
        model.setSubstanceUnits("moles_per_liter")

    spcs = []
    for sp_i in sp_comp:
        spcs.append(model.createSpecies())
        spcs[-1].setName(sp_i)
        spcs[-1].setId(sp_i)
        spcs[-1].setCompartment('Basic')
        spcs[-1].setInitialAmount(float(conc[sp_i]))
        # if molar == False:
        if not molar:
            spcs[-1].setSubstanceUnits('molecules')
        # elif molar == True:
        elif molar:
            spcs[-1].setSubstanceUnits('moles_per_liter')
        else:
            spcs[-1].setHasOnlySubstanceUnits(True)

    rcks_var = []
    for i, _ in enumerate(kcs_vals):
        if len(kcs_vals[i]) == 1:
            rcks_var.append(model.createParameter())
            rcks_var[-1].setId(str(kcs_vals[i][0]))
            rcks_var[-1].setConstant(True)
            rcks_var[-1].setValue(ks_vals[i][0])
            dvar = sum([r_react[i][x] for x in r_react[i]])
            # if molar == False:
            if not molar:
                if dvar == 0:
                    rcks_var[-1].setUnits('molecules_per_second')
                elif dvar == 1:
                    rcks_var[-1].setUnits('per_second')
                elif dvar == 2:
                    rcks_var[-1].setUnits('per_molecules_second')
            # elif molar == True:
            elif molar:
                if dvar == 0:
                    rcks_var[-1].setUnits('molar_per_second')
                elif dvar == 1:
                    rcks_var[-1].setUnits('per_second')
                elif dvar == 2:
                    rcks_var[-1].setUnits('per_molar_second')
            else:
                pass
            rcks_var[-1].setConstant(True)
        else:
            rcks_var.append(model.createParameter())
            rcks_var[-1].setId(str(kcs_vals[i][0]))
            rcks_var[-1].setConstant(True)
            rcks_var[-1].setValue(ks_vals[i][0])
            dvar = sum([r_react[i][x] for x in r_react[i]])
            # if molar == False:
            if not molar:
                if dvar == 0:
                    rcks_var[-1].setUnits('molecules_per_second')
                elif dvar == 1:
                    rcks_var[-1].setUnits('per_second')
                elif dvar == 2:
                    rcks_var[-1].setUnits('per_molecules_second')
            # elif molar == True:
            elif molar:
                if dvar == 0:
                    rcks_var[-1].setUnits('molar_per_second')
                elif dvar == 1:
                    rcks_var[-1].setUnits('per_second')
                elif dvar == 2:
                    rcks_var[-1].setUnits('per_molar_second')
            else:
                pass
            rcks_var[-1].setConstant(True)

            rcks_var.append(model.createParameter())
            rcks_var[-1].setId(str(kcs_vals[i][1]))
            rcks_var[-1].setConstant(True)
            rcks_var[-1].setValue(ks_vals[i][1])
            dvar = sum([r_react[i][x] for x in r_react[i]])
            # if molar == False:
            if not molar:
                if dvar == 0:
                    rcks_var[-1].setUnits('molecules_per_second')
                elif dvar == 1:
                    rcks_var[-1].setUnits('per_second')
                elif dvar == 2:
                    rcks_var[-1].setUnits('per_molecules_second')
            # elif molar == True:
            elif molar:
                if dvar == 0:
                    rcks_var[-1].setUnits('molar_per_second')
                elif dvar == 1:
                    rcks_var[-1].setUnits('per_second')
                elif dvar == 2:
                    rcks_var[-1].setUnits('per_molar_second')
            else:
                pass
            rcks_var[-1].setConstant(True)

    rxns = []
    for rvar in r_react:
        rxns.append(model.createReaction())
        rxns[-1].setId("Rf" + str(int(rvar) + 1))
        rxns[-1].setReversible(False)
        rbig_var = []
        for avar in r_react[rvar]:
            if int(r_react[rvar][avar]) > 0 and avar[0] != "-":
                rbig_var.append(rxns[-1].createReactant())
                rbig_var[-1].setSpecies(avar)
                rbig_var[-1].setStoichiometry(r_react[rvar][avar])

        pbig_var = []
        for bvar in r_prods[rvar]:
            if int(r_prods[rvar][bvar]) > 0 and bvar[0] != "-":
                pbig_var.append(rxns[-1].createProduct())
                pbig_var[-1].setSpecies(bvar)
                pbig_var[-1].setStoichiometry(r_prods[rvar][bvar])
        if len(ks_vals[rvar]) == 2:
            rxns.append(model.createReaction())
            rxns[-1].setId("Rb" + str(int(rvar) + 1))
            rxns[-1].setReversible(False)
            rbig_var = []
            for avar in r_react[rvar]:
                if int(r_react[rvar][avar]) > 0 and avar[0] != "-":
                    rbig_var.append(rxns[-1].createProduct())
                    rbig_var[-1].setSpecies(avar)
                    rbig_var[-1].setStoichiometry(r_react[rvar][avar])

            pbig_var = []
            for bvar in r_prods[rvar]:
                if int(r_prods[rvar][bvar]) > 0 and bvar[0] != "-":
                    pbig_var.append(rxns[-1].createReactant())
                    pbig_var[-1].setSpecies(bvar)
                    pbig_var[-1].setStoichiometry(r_prods[rvar][bvar])

    kinetic_law = []
    for rvar, _ in enumerate(rxns):
        kinetic_law.append(rxns[rvar].createKineticLaw())
        kinetic_law[-1].setMath(parseL3Formula(str(fprop[rvar]
                                                   ).replace("**", "^")))

    document.setModel(model)
    writeSBMLToFile(document, globals2.TO_CONVERT + ".xml")
    print(vbig, items)
    return [0, 0]
