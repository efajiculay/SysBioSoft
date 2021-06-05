from sympy import *
import numpy as np
from BioSANS2020.gui_functs.scrollable_text import *


done_parsing = set()
def process(x):
	global done_parsing
	xx =  x.strip("*").strip().replace("*/","/")
	if xx[0] == "/":
		xx = "1"+xx
	if xx[-1] == "/":
		xx = xx.strip("/")
	val = sympify(xx.strip("*"))
	if val in done_parsing:
		return None
	else:
		done_parsing.add(val)
	return val

def propExt(expr,prop):
	ex = str(expr)
	open = 0
	close = 0
	diff = 0
	
	collect = ""
	for v in ex:
		if v in ["+","-"] and diff == 0:
			if len(collect)>0:
				d = process(collect)
				if d != None:
					prop.append(d)
				collect = ""
				open = 0
				close = 0
		elif v == "(":
			open = open + 1
			collect = collect + v
		elif v == ")":
			close = close + 1
			collect = collect + v
		elif (v.isnumeric() or v == ".") and (diff == 0):
			if len(collect)>=3:
				if collect[-2:] == "**" or collect[-1] not in ["*","/"]:
					collect = collect + v
			elif len(collect)>=1:
				if collect[-1] not in ["*","/"," "]:
					collect = collect + v
			else:
				pass
		else:
			collect = collect + v	
		diff = open - close

	d = process(collect)
	if d != None:
		prop.append(d)
		
def termExt(expr):
	term = []
	ex = str(expr)
	open = 0
	close = 0
	diff = 0
	
	collect = ""
	last_sign = ""
	for v in ex:
		if v in ["+","-"] and diff == 0:
			if len(collect)>0:
				term.append(last_sign+collect)
				collect = ""
				open = 0
				close = 0
			last_sign = v
		elif v == "(":
			open = open + 1
			collect = collect + v
		elif v == ")":
			close = close + 1
			collect = collect + v
		else:
			collect = collect + v	
		diff = open - close

	term.append(last_sign+collect)
	return term
		
def get_prop_stoich(dxdt):
	prop = []
	dAdt = []
	for expro in dxdt:
		expr = expand(sympify(expro))
		propExt(expr,prop)
		dAdt.append(expr)
		
	w = prop
	V = [[0 for x in range(len(w))] for y in range(len(dAdt))]
	
	for i in range(len(dAdt)):
		for j in range(len(w)):
			s = 0
			for xx in termExt(dAdt[i]):
				x = sympify(xx)/w[j]
				try:
					s = s + float(x)
				except:
					pass
			V[i][j] = s
	
	return Matrix(V), Matrix(w)

def print_stoich_prop(dxdt):
	global done_parsing
	done_parsing = set()
	print()
	V, w = get_prop_stoich(dxdt)
	for t in V*w:
		print(t)
	print()

	for c in np.array(V):
		print([round(y,4) for y in c])
	print()
		
	for c in np.array(w):
		print(c)

def transform_to_rxn(x, dxdt, xo, ks, items):
	global done_parsing
	
	if items:
		text = prepare_scroll_text(items)
		ffrint = lambda x: text.insert(INSERT,x+"\n")
	else:
		ffrint = lambda x: print(" ".join([str(y) for y in x]),end="")		
	
	done_parsing = set()
	V, w = get_prop_stoich(dxdt)
	S = np.around(np.array(V).astype(float),3)
	Rxn = []
	Ksn = set()
	ind = 0 

	for col in S.T:
		R = ""
		P = ""
		for i in range(len(col)):
			if col[i] != 0:
				if col[i]<0:
					R = R + str(abs(col[i]))+" "+x[i] +" "+ "+ "
				else:
					P = P + str(abs(col[i]))+" "+x[i] +" "+ "+ "
		if R.strip() == "":
			R = "0 NONE"
		if P.strip() == "":
			P = "0 NONE"
		
		inSp = []
		for s in w[ind].free_symbols:
			ss = str(s)
			if ss in x:
				inSp.append(ss)
			else:
				Ksn.add(ss)
		inSp = ",".join(inSp)
		Rxn.append(R.strip("+ ")+" => "+P.strip("+ ")+", 1 ::::: lambda "+inSp+" : "+str(w[ind]))
		ind = ind + 1

	ffrint("Function_Definitions:")
	for k in Ksn:
		ffrint(k.strip()+" = type actual value") if k.strip() not in ks else ffrint(k.strip()+" = "+ks[k.strip()])
	for c in x:
		ffrint(c.strip()+"_ini = type actual value") if c.strip() not in xo else ffrint(c.strip()+"_ini = "+xo[c.strip()])
		
	ffrint("")
	ffrint("#REACTIONS")
	for rx in Rxn:
		ffrint(rx)
		
	ffrint("")
	ffrint("@CONCENTRATION")
	for c in x:
		ffrint(c+" , "+c.strip()+"_ini")
	ffrint("NONE, 1")
	return text

def odedxdt_to_topo(mfile,items):		
	dd = open(mfile,"r")
	print()
	x = []
	xo = {}
	ks = {}
	dxdt = []
	last = ""
	for xx in dd:
		if last == "ODE_DECLARATIONS" and xx.strip() not in ["INI_CONCENTRATIONS:","RATE_CONSTANTS:"]:
			if xx.strip() != "":
				row = xx.split("=")
				x.append(row[0].strip())
				dxdt.append(row[1])
		elif last == "INI_CONCENTRATIONS" and xx.strip() not in ["ODE_DECLARATIONS:","RATE_CONSTANTS:"]:
			if xx.strip() != "":
				row = xx.split("=")
				xo[row[0].strip()] = row[1].strip()
		elif last == "RATE_CONSTANTS"  and xx.strip() not in ["ODE_DECLARATIONS:","INI_CONCENTRATIONS:"]:
			if xx.strip() != "":
				row = xx.split("=")
				ks[row[0].strip()] = row[1].strip()	
		elif xx.strip() == "ODE_DECLARATIONS:":
			last = "ODE_DECLARATIONS"
		elif xx.strip() == "INI_CONCENTRATIONS:":
			last = "INI_CONCENTRATIONS"
		elif xx.strip() == "RATE_CONSTANTS:":
			last = "RATE_CONSTANTS"
			
	#print_stoich_prop(dxdt)
	return transform_to_rxn(x,dxdt,xo,ks,items)

#odedxdt_to_topo("NewText.txt")