import libsbml as Mysbml
from sympy import Symbol as newSymbol
from sympy import solve, sympify
from sympy import sympify
import mglobals as globals2
from sbmlMath import *
import inspect

OPERS_list = {"+","-","*","/","(",")",",","=",">","<",":"}
OPERS_list2 = {
	"abs","acos","arccos","acosh","arccosh","acot","arccot","acoth","arccoth","acsc","acsch","arccsch","arcsec","true","false","True","False",
	"asech","arcsech","asin","arcsin","atan","arctan","atanh","arctanh","ceil","ceiling","cos","cosh","cot","coth","csc","csch","and","not",
	"factorial","exp","floor","ln","log","log10","piecewise","pow","power","root","sec","sech","sqr","sqrt","sin","sinh","tan","None",
	"tanh","And","Not","Or","or","xor","eq","geq","gt","leq","lt","neq","plus","times","minus","divide","if","else","multiply","lambda"
}

def replace_crucial_funct(trep):
	return 	Add2spacesSep(trep) \
	.replace("and","And") \
	.replace("not","Not") \
	.replace("or","Or") \
	.replace("factOrial","factorial") \
	.replace("floOr","floor") \
	.replace("xOr","xor") \
	.replace("true","True") \
	.replace("false","False") \
	.replace(" ","") \
	.replace("+-","+ -")

def extractSpecies(modk):
	global OPERS_list, OPERS_list2
	here = " "+str(modk).replace("time-","emit-").replace("e+"," ").replace("E+"," ").replace("e-"," ").replace("E-"," ").replace("emit-","time-")+" "
	for x in OPERS_list:
		here = here.replace(x,"  ")
	here = " "+here+" "
	for x in OPERS_list2:
		here = here.replace(" "+x+" "," ")
	here = here.split()	
	sp = set()
	for x in here:
		try :
			float(x)
		except:
			try:
				float(eval(x))
			except:
				sp.add(x)
	return ",".join(sp)
	
def extractParNum(modk):
	global OPERS_list, OPERS_list2
	here = " "+str(modk)+" "
	for x in OPERS_list:
		here = here.replace(x,"  ")
	
	here = " "+here+" "
	for x in OPERS_list2:
		here = here.replace(" "+x+" "," ")
	here = here.split()	
	sp = set()
	for x in here:
		try :
			float(x)
			sp.add(x)
		except:
			pass
				
	return sp

def extractVarFunc(ss):
	v = ""
	lastComma = 0
	oper = {"+","-","*","/","(",")"}
	p = 0
	for x in ss+")":
		if x == ",":
			lastComma = p
		elif x not in oper:
			pass
		else:
			return [ss[0:lastComma],ss[lastComma+1:]]
		p = p + 1
	return ["",""]
	
def extractFunction(pforms,Rparams,compartments,functions,functions_str):
	global OPERS_list, OPERS_list2
	spcorm = pforms.replace("("," ").replace(")"," ").replace(","," ").replace("*"," ").split()	
	modk = ""
	ind2 = 0
	for sp in spcorm:
		if sp in Rparams:
			pass
		elif sp in compartments:
			pass
		elif sp in functions:
			npar = len(inspect.getargspec(functions[sp])[0])
			ind = spcorm.index(sp)
			ind2 = pforms.index(sp)
			ss = []
			for y in spcorm[ind+1:]:
				try:
					float(y)
				except:
					if y not in OPERS_list and y not in OPERS_list2 and y not in functions:
						ss.append(newSymbol(y))
			if len(ss) == npar:
				try:
					modk = str(functions[sp](*ss))
				except:
					modk = functions_str[sp]
					modk = FunctRedefineVar(modk,ss)
			else:
				modk = "lambda "+",".join([str(s) for s in ss])+" : "+pforms[ind2:]
				try:
					modk = str(eval(modk)(*ss))
				except:
					modk = str(eval(modk[0:-1])(*ss))
			break
		else:
			pass
	if len(modk)>2:
		modk = "1*("+modk+")"
	modk = pforms[:ind2]+modk
	return modk
	
def parSubstitution(modk,parameters,index=0):
	global OPERS_list
	modk = str(modk)
	for y in OPERS_list:
		modk = modk.replace(y," "+y+" ")
	modk = " "+modk+" "	
	for y in parameters:
		if parameters[y][0] is not None:
			modk = modk.replace(" "+y+" ",str(parameters[y][index]))
	return modk.replace(" ","")
	
def Add2spacesSep(modk):
	modk = str(modk)
	for y in OPERS_list:
		modk = modk.replace(y,"  "+y+"  ")
	return modk
	
def SpSubstitution(modk,parameters):
	global OPERS_list
	modk = str(modk)
	for y in OPERS_list:
		modk = modk.replace(y," "+y+" ")
	modk = " "+modk+" "	
	for y in parameters:
		if parameters[y] is not None:
			modk = modk.replace(" "+y+" ",str(parameters[y]))
	return modk.replace(" ","")
	
def varSubstitution(modk,Rparams,parameters,compartments,RateRules):
	global OPERS_list
	modk = str(modk)
	for y in OPERS_list:
		modk = modk.replace(y,"  "+y+"  ")
	modk = " "+modk+" "
	
	for y in Rparams:
		if y not in RateRules:
			modk = modk.replace(" "+y+" ","("+str(Rparams[y])+")")
	for y in parameters:
		if y not in RateRules:
			modk = modk.replace(" "+y+" ","("+str(parameters[y][0])+")")
	for y in compartments:
		modk = modk.replace(" "+y+" ","("+str(compartments[y][0])+")")
	return modk.replace(" ","")
	
def FunctRedefineVar(modk,ss):
	spset = [ b.strip() for b in modk.split(":")[0].split("lambda")[1].split(",") ]
	for y in OPERS_list:
		modk = modk.replace(y," "+y+" ")
	modk = " "+modk+" "
	for s1 in range(len(ss)):
		modk = modk.replace(" "+spset[s1]+" ",str(ss[s1]))
	return modk.split(":")[1]
	
def process_sbml(file,molar=False):
	global OPERS_list, OPERS_list2
	
	fftopofile = open(file+".topo","w")
	
	SbmlUnits = {
	0:"ampere",1:"avogadro",2:"becquerel",3:"candela",5:"coulomb",6:"dimensionless",7:"farad",8:"gram",9:"gray",10:"henry",11:"hertz",13:"joule",14:"katal",15:"kelvin",18:"litre",19:"lumen",20:"lux",22:"metre",23:"mole",24:"newton",25:"ohm",26:"pascal",27:"radian",28:"second",29:"siemens",30:"sievert",31:"steradian",32:"tesla",34:"watt",35:"weber"
	}
	
	reader = Mysbml.SBMLReader()
	document = reader.readSBML(file)
	model = document.getModel()
	errors = document.getNumErrors();
	
	time_var = None
		
	units = {}
	for x in model.getListOfUnitDefinitions():
		c = x.getListOfUnits()
		units[x.getId()] = []
		for y in c:
			F = ( y.getMultiplier() * 10 ** y.getScale() ) ** y.getExponent()
			units[x.getId()].append([F,SbmlUnits[y.getKind()], y.getExponent()])
	
	compartments = {}
	constant_comp = {}
	nonConstant_comp = {}
	orig_size = {}
	for x in model.getListOfCompartments():
		if x.isSetSize():
			orig_size[x.getId()] = x.getSize()
			if not molar:
				compartments[x.getId()] = [ x.getSize(), x.getUnits() ]
			else:
				compartments[x.getId()] = [ 1, x.getUnits() ]
		else:
			orig_size[x.getId()] = 1
			compartments[x.getId()] = [ 1, x.getUnits() ]
		if x.getConstant() == True:
			constant_comp[x.getId()] = compartments[x.getId()]
		else:
			nonConstant_comp[x.getId()] = compartments[x.getId()]
				
	reactions = {}
	for x in model.getListOfReactions():
		reactions[x.getId()] = x
				
	species = {}
	species_comp = {}
	constant_species = {}
	HasOnlySUnits = {}
	for x in model.getListOfSpecies():
		species[x.getId()] = x
		species_comp[x.getId()] = x.getCompartment()
		HasOnlySUnits[x.getId()] = x.getHasOnlySubstanceUnits()
		if x.getConstant()==True:
			constant_species[x.getId()] = True
			
	SpInitialConc = {}
	for x in species:
		sp = species[x]
		if sp.isSetInitialConcentration():
			val = sp.getInitialConcentration();
			SpInitialConc[x] = val
		elif sp.isSetInitialAmount():
			val = sp.getInitialAmount()	
			SpInitialConc[x] = val
				
	parameters = {}
	nonConstPar = set()
	constant_par = {}
	for x in model.getListOfParameters():
		if x.isSetValue():
			ss = x.getValue()
		else:
			ss = None
		parameters[x.getId()] = [ss, x.getUnits()]
		if x.getConstant()==False:
			nonConstPar.add(x.getId())
		else:
			constant_par[x.getId()] = parameters[x.getId()]
			
	AssignRules = {}
	RateRules = {}
	AlgebrRules = []
	for x in model.getListOfRules():
		ss = replace_crucial_funct(Mysbml.formulaToString(x.getMath()))
		if x.isAssignment():
			ss = parSubstitution(ss,constant_par)
			AssignRules[x.getVariable()] = ss
		elif x.isRate():
			RateRules[x.getVariable()] = ss
		elif x.isAlgebraic():
			AlgebrRules.append(ss)
					
	functions = {}
	functions_str = {}
	fftopofile.write("Function_Definitions:\n")
	for x in model.getListOfFunctionDefinitions():
		ss = replace_crucial_funct(Mysbml.formulaToString(x.getMath()))
		ss = ss.replace("lambda(","")[:-1]
		ss = extractVarFunc(ss)
		ss = " : ".join(ss)
		functions[x.getId()] = eval("lambda "+ss)
		functions_str[x.getId()] = "lambda "+ss
		globals2.execFunctions.append(x.getId()+" = lambda "+ss)
		exec(x.getId()+" = lambda "+ss,globals())
		fftopofile.write(x.getId()+" = lambda "+ss+"\n")
		OPERS_list2.add(x.getId())
	
	InitialAssign = {}
	for x in model.getListOfInitialAssignments():
		sp = x.getId()
		ss = replace_crucial_funct(Mysbml.formulaToString(x.getMath()).replace("lambda(",""))
		if sp in species:
			ss = str(compartments[species_comp[sp]][0])+"*("+ss+")"
		try: 
			ss = eval(parSubstitution(ss,parameters))
		except:
			ss = parSubstitution(ss,parameters)
		InitialAssign[sp] = ss
		
		if sp in constant_par:
			constant_par[sp][0] = ss
			parameters[sp][0] = ss	
		elif sp in nonConstPar:
			parameters[sp][0] = ss		
		elif sp in compartments:
			orig_size[sp] = ss
			if not molar: 
				compartments[sp][0] = ss
			else:
				compartments[sp][0] = 1
		elif sp in constant_comp:
			constant_comp[sp][0] = ss
		elif sp in nonConstant_comp:
			nonConstant_comp[sp][0] = ss
					
	Events = {}
	EventAssign = {}
	Delays = {}
	IniValTrig = {}
	Priorities = []
	for x in model.getListOfEvents():
		uVFTT = x.getUseValuesFromTriggerTime()
		IniTV = x.getTrigger().getInitialValue()
		PersT = x.getTrigger().getPersistent()
		
		Prior = ""
		if x.isSetPriority():
			ss = replace_crucial_funct(Mysbml.formulaToString(x.getPriority().getMath()))
			#ss = SpSubstitution(ss,SpInitialConc)
			Prior = ss #eval(ss)
			Priorities.append(Prior)
		
		ss = replace_crucial_funct(parSubstitution(Mysbml.formulaToString(x.getTrigger().getMath()),constant_par))
		try:
			for sp in species:
				ss = Add2spacesSep(ss).replace(" "+sp+" ",sp+"/"+str(compartments[species_comp[sp]][0])).replace(" ","")
		except:
			pass
			
		mods = extractSpecies(ss)
		for tt in mods.split(","):
			if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
				time_var  = tt
				break			
					
		Events[x.getId()] = ss
		if not x.isSetDelay():
			for y in x.getListOfEventAssignments():
				IniValTrig[y.getVariable()] = IniTV 
				try:
					ss = str(eval(parSubstitution(Mysbml.formulaToString(y.getMath()),constant_par)))
				except:
					ss = parSubstitution(Mysbml.formulaToString(y.getMath()),constant_par)
					ss = parSubstitution(ss,constant_comp)
				ss = replace_crucial_funct(ss)
				mods = extractSpecies(ss)
				for tt in mods.split(","):
					if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
						time_var  = tt
						break
						
				sp = y.getVariable()
				keyS = "status_"+sp
				
				if sp in species_comp:
					ss = str(compartments[species_comp[sp]][0])+"*("+ss+")"
					
				if Prior == "":
					estatus = " : True if "+Events[x.getId()]+" else False"
					if PersT:
						express = " : "+ss+" if "+Events[x.getId()]+" and not "+keyS +" else "+y.getVariable()
					else:
						express = " : "+ss+" if "+Events[x.getId()]+" and "+keyS +" else "+y.getVariable()
				else:
					keyS = keyS +"_"+str(Priorities.index(Prior))
					if PersT:
						estatus = " : (True if "+Events[x.getId()]+" else False"+","+Prior+")"
						express = " : ("+ss+" if "+Events[x.getId()]+" and not "+keyS +" else "+y.getVariable()+","+Prior+")"
					else:
						estatus = " : (True if "+Events[x.getId()]+" else False"+","+Prior+",1)"
						express = " : ("+ss+" if "+Events[x.getId()]+" and "+keyS +" else "+y.getVariable()+","+Prior+",1)"
						
				estatus = "lambda "+extractSpecies(estatus)+estatus			
				express = "lambda "+extractSpecies(express)+express				
				if keyS in EventAssign:
					EventAssign[keyS].append(estatus)
				else:
					EventAssign[keyS] = [estatus]
				if y.getVariable() in EventAssign:
					EventAssign[y.getVariable()].append(express)
				else:
					EventAssign[y.getVariable()] = [express]
		else:
			Delays[x.getId()] = replace_crucial_funct(parSubstitution(Mysbml.formulaToString(x.getDelay().getMath()),constant_par))
			tdepDelay = False
			mods = extractSpecies(Delays[x.getId()])
			for tt in mods.split(","):
				if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
					time_var  = tt
					tdepDelay = True
					break			
			
			for y in x.getListOfEventAssignments():
				IniValTrig[y.getVariable()] = IniTV 
			
				sp = y.getVariable()
				keyT = "timer_"+sp
				keyS = "status_"+sp
				keyS2 = "status2_"+sp
				keyF = "finish_"+sp
				keyC = "dtime_"+sp
				keyD = "delay_"+sp
				TheMath = replace_crucial_funct(parSubstitution(Mysbml.formulaToString(y.getMath()),constant_par))
			
				if sp in species_comp:
					TheMath = str(compartments[species_comp[sp]][0])+"*("+TheMath+")"
					
				mods = extractSpecies(TheMath)
				for tt in mods.split(","):
					if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
						time_var  = tt
						tdepDelay = True
						break					
				
				if tdepDelay == True:	
					estatus2 = " : True if "+Events[x.getId()]+" else False" 
					estatus2 = "lambda "+extractSpecies(estatus2)+estatus2		
					DelayVal = " : "+Delays[x.getId()]+" if "+Events[x.getId()]+" and "+keyF+" and not "+keyS2+" else None"
					DelayVal = "lambda "+extractSpecies(DelayVal)+DelayVal
					DelayCon = keyT+" >= "+keyD
				else:
					DelayCon = keyT+" >= "+Delays[x.getId()]
				
				finishFr = " : 1 if "+DelayCon+" else 0 if "+Events[x.getId()]+" else None"
				finishFr = "lambda "+extractSpecies(finishFr)+finishFr

				estatus = " : True if "+DelayCon+" else False" 
				estatus = "lambda "+extractSpecies(estatus)+estatus	
				
				express1 = " : "+keyT+" if "+Events[x.getId()]+" or not "+keyF+" else 0"
				express1 = "lambda "+extractSpecies(express1)+express1
				
				CurrentT = " : "+TheMath+" if "+Events[x.getId()]+" and "+keyF+" else None"
				CurrentT = "lambda "+extractSpecies(CurrentT)+CurrentT			
				
				if PersT:
					keySCond = "not "+keyS
					pass
				else:
					estatus = " : True if "+Events[x.getId()]+" else False" 
					estatus = "lambda "+extractSpecies(estatus)+estatus
					keySCond = keyS
					
				if uVFTT:
					express2 = " : "+keyC+" if "+DelayCon+" and "+keySCond+" else "+sp
				else:
					express2 = " : "+TheMath+" if "+DelayCon+" and "+keySCond+" else "+sp
					
				express2 = "lambda "+extractSpecies(express2)+express2			
							
				if y.getVariable() in EventAssign:
					EventAssign[keyT].append(express1)
					EventAssign[sp].append(express2)
					EventAssign[keyS].append(estatus)
					EventAssign[keyF].append(finishFr)
					EventAssign[keyC].append(CurrentT)
					if tdepDelay == True:
						EventAssign[keyS2].append(estatus2)
						EventAssign[keyD].append(DelayVal)
				else:
					EventAssign[keyT] = [express1]
					EventAssign[sp] = [express2]
					EventAssign[keyS] = [estatus]
					EventAssign[keyF] = [finishFr]
					EventAssign[keyC] = [CurrentT]
					if tdepDelay == True:	
						EventAssign[keyD] = [DelayVal]
						EventAssign[keyS2] = [estatus2]
						
			
	for x in AssignRules:
		ss = AssignRules[x]
		ss = parSubstitution(ss,parameters)
		if x in compartments:
			compartments[x][0] = eval(ss)
		elif x in parameters:
			try:
				parameters[x][0] = eval(ss)	
			except:
				pass
	
	NONE_added = False
	fftopofile.write("\n")
	if len(orig_size) == 1:
		for cx in orig_size:
			if cx in constant_comp:
				fftopofile.write("#Reactions Volume = "+str(orig_size[cx])+"\n")
			else:
				fftopofile.write("#Reactions\n")
	else:
		fftopofile.write("#Reactions\n")
	Rparams = {}
	UseSpecies = set()
	fastRxnExp = {}
	for x in reactions:
		Rparams = {}
		Factor = 1
		react = ""
		row = []
		if (reactions[x].isSetKineticLaw()):
			for kx in reactions[x].getKineticLaw().getListOfParameters():
				Rparams[kx.getId()] = kx.getValue()
				
		for rx in reactions[x].getListOfReactants():
			key = rx.getSpecies()
			UseSpecies.add(key)
			if rx.isSetStoichiometryMath():
				dd = rx.getStoichiometryMath().getMath()
				dd = replace_crucial_funct(Mysbml.formulaToString(dd))
				dd = varSubstitution(dd,Rparams,parameters,constant_comp,RateRules)
				dd = eval(dd)
				react = react + str(dd) + " "  + key + " + "
			else:			
				react = react + str(rx.getStoichiometry()) + " "  + key + " + "
			if HasOnlySUnits[key]:
				pass
			else:
				if species_comp[key] in constant_comp:
					Factor = Factor*(1/compartments[species_comp[key]][0])**rx.getStoichiometry()
				else:
					Factor = Factor*(1/newSymbol(species_comp[key]))**rx.getStoichiometry()

		if len(react) == 0:
			NONE_added = True
			react = "0 NONE => "
		else:	
			react = react[0:-2]+" => "
		
		produ = ""
		Pactor = 1
		for rx in reactions[x].getListOfProducts():
			key = rx.getSpecies()
			UseSpecies.add(key)
			if rx.isSetStoichiometryMath():
				dd = rx.getStoichiometryMath().getMath()
				dd = replace_crucial_funct(Mysbml.formulaToString(dd))
				dd = varSubstitution(dd,Rparams,parameters,constant_comp,RateRules)
				dd = eval(dd)
				produ = produ + str(dd) + " "  + key + " + "
			else:
				produ = produ + str(rx.getStoichiometry()) + " "  + key + " + "
			if HasOnlySUnits[key]:
				pass
			else:
				if species_comp[key] in constant_comp:
					Pactor = Pactor*(1/compartments[species_comp[key]][0])**rx.getStoichiometry()
				else:
					Pactor = Pactor*(1/newSymbol(species_comp[key]))**rx.getStoichiometry()

		if len(produ) == 0:
			NONE_added = True
			produ = "0 NONE"
			react = react + produ
			react = react + ",1   :::::   "		
		else:
			react = react + produ
			react = react[0:-2] + ",1   :::::   "		
				
		kforts = reactions[x].getKineticLaw()
		pforms = replace_crucial_funct(Mysbml.formulaToString(kforts.getMath()))
		
		try:
			modk = extractFunction(pforms,Rparams,constant_comp,functions,functions_str)
		except:
			modk = ""
		
		if len(modk) == 0:
			modk = pforms
		else:
			pass			
			
		modk = varSubstitution(modk,Rparams,constant_par,constant_comp,RateRules)	
		
		if reactions[x].getReversible() == True:
			if modk.find("+ -")>-1:
				modk = modk.split("+ -")
			else:
				if modk[0] == "-":
					modk = modk[1:].replace("-","+-").split("+")
					modk[0] = "-"+modk[0]
				else:
					modk = modk.replace("-","+-").split("+")
			
			if modk[0] == modk[-1]:	
				modk = "("+modk[0]+"*("+str(Factor)+")"
			elif reactions[x].getFast() == True:
				if modk[0][0] == "-" or modk[0][1] == "-":
					modk = "1000*("+modk[0]+"*("+str(Pactor)+")+("+str(Factor)+")*"+modk[-1]+")"
				else:
					modk = "1000*("+modk[0]+"*("+str(Factor)+")+("+str(Pactor)+")*"+modk[-1]+")"
				sphere = extractSpecies(modk).split(",")
				for ss in sphere:
					fastRxnExp[ss] = modk
			else:
				if modk[0][0] == "-" or modk[0][1] == "-":
					modk = "("+modk[0]+"*("+str(Pactor)+")+("+str(Factor)+")*"+modk[-1]+")"
				else:
					modk = "("+modk[0]+"*("+str(Factor)+")+("+str(Pactor)+")*"+modk[-1]+")"
				
			if modk.count("(")>modk.count(")"):
				modk = modk + ")"
			mods = extractSpecies(modk)
			modk = "lambda "+mods+" : "+modk
			for tt in mods.split(","):
				if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
					time_var  = tt
					break
		else:		
			mods = extractSpecies(modk)
			if reactions[x].getFast() == True:
				sphere = mods.split(",")
				for ss in sphere:
					fastRxnExp[ss] = modk	
					
			modk = "lambda "+mods+" :"+str(Factor)+"*"+modk
			for tt in mods.split(","):
				if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
					time_var  = tt
					break			
		react = react + modk	
		fftopofile.write(react+"\n")

	for y in range(len(AlgebrRules)):
		ss = varSubstitution(AlgebrRules[y],{},constant_par,constant_comp,RateRules)
		spss = extractSpecies(ss).split(",") 
		try:
			ss = eval("lambda "+",".join(spss)+" : "+ ss)
			sf = [ newSymbol(x) for x in spss ]
			ss = ss(*sf)
		except:
			pass
		AlgebrRules[y] = "Eq("+str(ss)+",0)"
	
	modifiersSp = {}
	for x in reactions:		
		for rx in reactions[x].getListOfModifiers():
			UseSpecies.add(rx.getSpecies())
			modifiersSp[rx.getSpecies()] = newSymbol(rx.getSpecies())
			fftopofile.write("0 NONE => " + rx.getSpecies()+", 0"+"\n")	
			NONE_added = True
		
	AlgebrRules = sympify(AlgebrRules)

	for x in EventAssign:
		if x not in species and x not in parameters and x not in compartments and x.find("status")==-1 and x.find("finish")==-1 and x.find("dtime")==-1 and x.find("delay")==-1:
			UseSpecies.add(x)
			fftopofile.write("0 NONE => " + x +", 1\n")	
			NONE_added = True			
					
	for x in RateRules:
		if x in parameters or x in compartments or x in species:
			UseSpecies.add(x)
			modk = extractFunction(RateRules[x],{},compartments,functions,functions_str)	

			mods = extractSpecies(RateRules[x])
			for tt in mods.split(","):
				if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
					time_var  = tt
					break			
			
			if modk == "":
				modk = RateRules[x]
			try:
				float(modk)
				if x in species:
					if HasOnlySUnits[x]:
						pass
					else:
						comp = compartments[species_comp[x]][0]
						modk = str(eval(modk+"*"+str(comp)))
				fftopofile.write("0 NONE"+" => "+ x +", "+modk+"\n")
			except:
				modk = modk.replace(x,"$"+x)
				modk = varSubstitution(modk,{},constant_par,constant_comp,RateRules)
				modk = "lambda "+extractSpecies(modk)+" : "+modk
				modk = modk.replace("$"+x,x)
				if x in species:
					if HasOnlySUnits[x]:
						pass
					else:
						comp = compartments[species_comp[x]][0]
						modk = modk+"*"+str(comp)
					fftopofile.write("0 NONE"+" => "+ x + ",1   :::::   " + modk+"\n") 
				else:
					fftopofile.write("0 NONE"+" => "+ x + ",1   :::::   " + modk+"\n") 
		else:
			pass
			#modk = RateRules[x]
			#fftopofile.write("0 NONE"+" => "+ x +", "+modk+"\n")
			
		NONE_added = True
				
	for x in AssignRules:
		if x in species or x in parameters or x in compartments:
			UseSpecies.add(x)
			fftopofile.write("0 NONE"+" => "+ x + ",0"+"\n")
		NONE_added = True
		ss = varSubstitution(AssignRules[x],Rparams,parameters,compartments,RateRules)
		sss = extractSpecies(ss)
		for tt in sss.split(","):
			if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
				time_var  = tt
				break	
	
	timevar = None
	for x in InitialAssign:
		ss = varSubstitution(InitialAssign[x],Rparams,parameters,compartments,RateRules)
		sss = extractSpecies(ss)
		for tt in sss.split(","):
			if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
				time_var  = tt
				timevar = 0
				break
		
	for x in nonConstant_comp:
		if x not in RateRules and x not in AssignRules and x not in AlgebrRules:
			fftopofile.write("0 NONE"+" => "+ x + ",0"+"\n")
		NONE_added = True
		
	for x in nonConstPar:
		if x not in RateRules and x not in AssignRules:
			fftopofile.write("0 NONE"+" => "+ x + ",0"+"\n")
		NONE_added = True			
		
	if len(reactions) == 0:
		for x in species:
			sp = species[x]
			if x in AssignRules:
				pass
			elif x in RateRules:
				pass	
			else:
				UseSpecies.add(x)
				modifiersSp[x] = newSymbol(x)
				fftopofile.write(x+" => 0 NONE, 0.0\n")
				
		for x in parameters:
			if x in InitialAssign and x not in RateRules:
				fftopofile.write(x+" => 0 NONE, 0.0\n")
				
		if len(nonConstPar) == 0:
			for x in parameters:
				fftopofile.write(x+" => 0 NONE, 0.0\n")
			
		NONE_added = True
		
	for x in species:
		if x not in UseSpecies:
			fftopofile.write(x+" => 0 NONE, 0.0\n")
			NONE_added = True
			
	if time_var != None:
		if timevar == None:
			fftopofile.write("0 NONE => "+time_var+", 1\n")
		elif timevar == 0:
			fftopofile.write("0 NONE => "+time_var+", 0\n")
		
	fftopofile.write("\n@Concentrations\n")		
	
	for x in species:
		sp = species[x]
		if x in fastRxnExp and sp.getBoundaryCondition() == False:
			ss = extractSpecies(fastRxnExp[x])
			spsum = ""
			sbsum = ""
			val = 0
			bal = 0
			syms = []
			noBoundary = True
			for s in ss.split(","):
				syms.append(newSymbol(s))
				if species[s].isSetInitialConcentration():
					vv = species[s].getInitialConcentration()
				elif species[s].isSetInitialAmount():
					vv = species[s].getInitialAmount()
					
				if species[s].getBoundaryCondition() == False:
					val = val + vv
					spsum = spsum + s +"+"
				else:
					noBoundary = False
					sbsum = s
					bal = bal + vv
					
			if noBoundary:
				spsum = spsum[0:-1]
				eqst = sympify(["Eq("+fastRxnExp[x]+",0)", "Eq("+spsum+","+str(val)+")"])
			else:
				eqst = sympify(["Eq("+fastRxnExp[x]+",0)", "Eq("+sbsum+","+str(bal)+")"])
			sol = solve(eqst,syms)
			if sol:
				val1 = sol[newSymbol(x)]
			else:
				val1 = 0

		elif sp.isSetInitialConcentration():
			val1 = sp.getInitialConcentration();
		elif sp.isSetInitialAmount():
			val1 = sp.getInitialAmount()
		else:
			val1 = ""
		if x in InitialAssign:
			ss = varSubstitution(InitialAssign[x],Rparams,parameters,compartments,RateRules)
			for ini in SpInitialConc:
				ss = ss.replace(ini,str(SpInitialConc[ini]))
			try:
				val1 = eval(ss)
			except:
				val1 = 0
				val2 = ss
		if x in AssignRules:
			ss = varSubstitution(AssignRules[x],Rparams,parameters,compartments,RateRules)
			sss = extractSpecies(ss)
			for tt in sss.split(","):
				if tt not in species and tt not in parameters and tt not in compartments and tt.strip()!="" and tt.find("delay")==-1:
					time_var  = tt
					break			
			try:
				val2 = eval(ss)
				if float(val2):
					if val1 == "":
						val1 = val2
			except:
				val2 = ss
				if val1 == "":
					val1 = 0

		if x not in EventAssign:
			cmod = ", "			
			if x in modifiersSp or x not in UseSpecies and x not in constant_species:
				ssphere = newSymbol(x)
				dd = solve(AlgebrRules,ssphere)
				if dd:
					cmod = cmod + "lambda "+extractSpecies(str(dd[ssphere]))+" : "+str(dd[ssphere])
			if cmod == ", ":
				cmod = ""

			if sp.getBoundaryCondition() == False and x not in AssignRules:
				fftopofile.write(sp.getId()+" , "+str(val1)+cmod+"\n");
			elif x in AssignRules:
				newFunct = extractFunction(str(val2),Rparams,compartments,functions,functions_str)
				if newFunct == "":
					newFunct = val2				
				fftopofile.write(sp.getId()+" ,"+str(val1)+", lambda "+extractSpecies(newFunct)+": "+str(newFunct)+str(cmod)+"\n");		
			elif sp.getBoundaryCondition() == True and x in RateRules:
				fftopofile.write(sp.getId()+" , "+str(val1)+cmod+"\n");
			elif sp.getBoundaryCondition() == True and x not in RateRules and x not in AssignRules and cmod.strip()!="" and x not in EventAssign:
				fftopofile.write(sp.getId()+" , "+str(val1)+cmod+"\n");
			else:
				fftopofile.write(sp.getId()+" , "+str(val1)+", lambda "+sp.getId()+": "+str(val1)+cmod+"\n");
		else:
			for ig in range(len(EventAssign[x])):
				cmod = ", "		
				cmod = cmod + EventAssign[x][ig]

				if sp.getBoundaryCondition() == False:
					fftopofile.write(sp.getId()+" , "+str(val1)+cmod+"\n");
				else:
					ss = cmod[1:].split(":")
					lamPa = ss[0].replace("lambda","").strip()
					cond = ss[1].split("if")[1].split("and")[0].strip()
					conPa = ss[1].split("if")[0].strip()
					fftopofile.write(sp.getId()+" , "+str(val1)+", lambda "+lamPa +": None if not "+cond+" else "+conPa+"\n");		
			
	for x in RateRules:
		if x in parameters and x not in EventAssign:
			fftopofile.write(x+","+str(parameters[x][0])+"\n")	
			NONE_added = True
		elif x in compartments:
			fftopofile.write(x+" , "+str(compartments[x][0])+"\n")	
			NONE_added = True
		else:
			pass
			#fftopofile.write(x+" , 1\n")	
			#NONE_added = True
			
	writtenEvent = {}
	for x in EventAssign:
		if x not in species:
			for ih in range(len(EventAssign[x])):
				if x in parameters:
					ss = str(parameters[x][0])
				elif x in nonConstant_comp:
					ss = str(nonConstant_comp[x][0])
				else:
					if x.find("finish")>-1 or x.find("delay")>-1:
						ss = "1"
					elif x.find("status")>-1:
						if IniValTrig[x.split("_")[1]] == True:
							ss = "1"
						elif IniValTrig[x.split("_")[1]] == False:
							ss = "0"
					else:
						ss = "0"
				Eev = x +","+ss+","+EventAssign[x][ih]+"\n"
				if Eev not in writtenEvent:
					writtenEvent[Eev] = True
					fftopofile.write(Eev)	
			
	for x in AssignRules:
		if x not in species:
			try:
				float(AssignRules[x])
				fftopofile.write(x +", "+AssignRules[x]+"\n")	
			except:
				ss = varSubstitution(AssignRules[x],Rparams,constant_par,compartments,RateRules)
				fftopofile.write(x +", 0, lambda "+extractSpecies(ss)+" : "+ss+"\n")	
				
	for x in InitialAssign:
		if x not in species and x not in RateRules:
			try:
				float(AssignRules[x])
				fftopofile.write(x +", "+InitialAssign[x]+"\n")	
			except:
				ss = varSubstitution(InitialAssign[x],Rparams,parameters,compartments,RateRules)
				if not molar:
					fftopofile.write(x +", 0, lambda "+extractSpecies(ss)+" : "+ss+"\n")
			
	for x in nonConstPar:
		if x not in RateRules and x not in AssignRules:
			ssphere = newSymbol(x)
			dd = solve(AlgebrRules,ssphere)	
			ss = str(parameters[x][0] if parameters[x][0] is not None else 0)
			try:
				fftopofile.write(x +","+ss+","+ "lambda "+extractSpecies(str(dd[ssphere]))+" : "+str(dd[ssphere])+"\n")	
			except:
				fftopofile.write(x +","+ss+"\n")
		NONE_added = True
		
	for x in nonConstant_comp:
		if x in InitialAssign and x not in RateRules and x not in AssignRules:
			fftopofile.write(x +"," +str(eval(str(nonConstant_comp[x][0])))+"\n" "\n")
			NONE_added = True
		elif x not in RateRules and x not in AssignRules and x not in EventAssign:
			ssphere = newSymbol(x)
			dd = solve(AlgebrRules,ssphere)
			fftopofile.write(x +"," +str(eval(str(nonConstant_comp[x][0])))+","+ "lambda "+extractSpecies(str(dd[ssphere]))+" : "+str(dd[ssphere])+"\n" "\n")
			NONE_added = True
			
	if len(reactions) == 0:
		if len(nonConstPar) == 0:
			for x in parameters:
				fftopofile.write(x+","+str(parameters[x][0])+"\n")
		
	if time_var != None:
		fftopofile.write(time_var+", 0\n")
		NONE_added = True
		
	if NONE_added:
		fftopofile.write("NONE , 1.0\n")
		
		
	fftopofile.close()