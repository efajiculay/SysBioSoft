#import sys
#import os
#sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.gui_functs.scrollable_text import *
from BioSANS2020.propagation.propensity import *
from sympy import *
from libsbml import *
from BioSANS2020.myglobal import mglobals as globals2

	
def topo_to_sbml(Sp,Ks,conc,Rr,Rp,V,Vm,items=None,molar=False):
	
	Cs = {}
	for x in Sp:
		Cs[x] = Symbol(x, real = True, negative=False)
	
	KCs = []
	for i in range(len(Ks)):
		row = []
		if len(Ks[i])==1:
			key = 'kf'+str(i+1)
			row.append(Symbol(key, real = True, positive=True))
		else:
			key = 'kf'+str(i+1)
			row.append(Symbol(key, real = True, positive=True))
			key = 'kb'+str(i+1)
			row.append(Symbol(key, real = True, positive=True))		
		KCs.append(row)  
			
	document = SBMLDocument(3, 1)
	model = document.createModel() 
	model.setTimeUnits("second")
		
	comp =  model.createCompartment()
	comp.setName('BasicCompartment')
	comp.setId('Basic')
	comp.setVolume(Vm)
	comp.setUnits('litre')
	comp.setConstant(True)
	comp.setSpatialDimensions(3)
	
	substance = model.createUnitDefinition()
	substance.setId('substance')
	unit = substance.createUnit()
	unit.setKind(UNIT_KIND_ITEM)
	
	if not molar:
		f = Matrix(propensity_vec(KCs,Cs,Rr,Rp))
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
		f = Matrix(propensity_vec_molar(KCs,Cs,Rr,Rp))
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
	for s in Sp:
		spcs.append(model.createSpecies())
		spcs[-1].setName(s)
		spcs[-1].setId(s)
		spcs[-1].setCompartment('Basic')
		spcs[-1].setInitialAmount(float(conc[s]))
		if molar == False:
			spcs[-1].setSubstanceUnits('molecules')
		elif molar == True:
			spcs[-1].setSubstanceUnits('moles_per_liter')
		else:
			spcs[-1].setHasOnlySubstanceUnits(True)
		
	rCks = []
	for i in range(len(KCs)):
		if len(KCs[i])==1:
			rCks.append(model.createParameter())
			rCks[-1].setId(str(KCs[i][0]))
			rCks[-1].setConstant(True)
			rCks[-1].setValue(Ks[i][0])
			d = sum([Rr[i][x] for x in Rr[i]])
			if molar == False:
				if d == 0:
					rCks[-1].setUnits('molecules_per_second')
				elif d == 1:
					rCks[-1].setUnits('per_second')
				elif d == 2:
					rCks[-1].setUnits('per_molecules_second')
			elif molar == True:
				if d == 0:
					rCks[-1].setUnits('molar_per_second')
				elif d == 1:
					rCks[-1].setUnits('per_second')
				elif d == 2:
					rCks[-1].setUnits('per_molar_second')
			else:
				pass
			rCks[-1].setConstant(True)
		else:
			rCks.append(model.createParameter())
			rCks[-1].setId(str(KCs[i][0]))
			rCks[-1].setConstant(True)
			rCks[-1].setValue(Ks[i][0])
			d = sum([Rr[i][x] for x in Rr[i]])
			if molar == False:
				if d == 0:
					rCks[-1].setUnits('molecules_per_second')
				elif d == 1:
					rCks[-1].setUnits('per_second')
				elif d == 2:
					rCks[-1].setUnits('per_molecules_second')
			elif molar == True:
				if d == 0:
					rCks[-1].setUnits('molar_per_second')
				elif d == 1:
					rCks[-1].setUnits('per_second')
				elif d == 2:
					rCks[-1].setUnits('per_molar_second')	
			else:
				pass
			rCks[-1].setConstant(True)
			
			rCks.append(model.createParameter())
			rCks[-1].setId(str(KCs[i][1]))
			rCks[-1].setConstant(True)
			rCks[-1].setValue(Ks[i][1])
			d = sum([Rr[i][x] for x in Rr[i]])
			if molar == False:
				if d == 0:
					rCks[-1].setUnits('molecules_per_second')
				elif d == 1:
					rCks[-1].setUnits('per_second')
				elif d == 2:
					rCks[-1].setUnits('per_molecules_second')
			elif molar == True:
				if d == 0:
					rCks[-1].setUnits('molar_per_second')
				elif d == 1:
					rCks[-1].setUnits('per_second')
				elif d == 2:
					rCks[-1].setUnits('per_molar_second')	
			else:
				pass
			rCks[-1].setConstant(True)
			
	rxns = []
	for r in Rr:
		rxns.append(model.createReaction())
		rxns[-1].setId("Rf"+str(int(r)+1))
		rxns[-1].setReversible(False)
		R = []
		for a in Rr[r]:
			if int(Rr[r][a])>0 and a[0]!="-":
				R.append(rxns[-1].createReactant())
				R[-1].setSpecies(a)
				R[-1].setStoichiometry(Rr[r][a])
		
		P = []
		for b in Rp[r]:
			if int(Rp[r][b])>0 and a[0]!="-":
				P.append(rxns[-1].createProduct())
				P[-1].setSpecies(b)
				P[-1].setStoichiometry(Rp[r][b])
		if len(Ks[r])==2:
			rxns.append(model.createReaction())
			rxns[-1].setId("Rb"+str(int(r)+1))
			rxns[-1].setReversible(False)
			R = []
			for a in Rr[r]:
				if int(Rr[r][a])>0 and a[0]!="-":
					R.append(rxns[-1].createProduct())
					R[-1].setSpecies(a)
					R[-1].setStoichiometry(Rr[r][a])
			
			P = []
			for b in Rp[r]:
				if int(Rp[r][b])>0 and a[0]!="-":
					P.append(rxns[-1].createReactant())
					P[-1].setSpecies(b)
					P[-1].setStoichiometry(Rp[r][b])		
	
	kinetic_law = []
	for r in range(len(rxns)):
		kinetic_law.append(rxns[r].createKineticLaw())
		kinetic_law[-1].setMath(parseL3Formula(str(f[r]).replace("**","^")))
			
	document.setModel(model)
	writeSBMLToFile(document,globals2.toConvert+".xml")
		
	return [0,0]