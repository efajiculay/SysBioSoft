import mglobals as globals2
from sbmlMath import *
import re
#from sympy import *
import numpy as np
import random as rm

reserve_events_words = {"t","time","status","status2","timer","finish","delay","dtime"}

INF = np.inf
NaN = np.nan
inf = INF
nan = NaN
pi = np.pi

def rateOf(x):
	global concP, RateSp, orasP, SiP
	indx = SiP.index(RateSp)
	if len(orasP)<2:
		return 0
	delt = orasP[-1] - orasP[-2]
	rate = (concP[-1][indx] - concP[-2][indx])/delt
	return rate

def delay(x,y):
	global orasP, concP, SiP, actualSp, concP2 
	if len(orasP)>0:
		last = orasP[-1]
		first = orasP[0]
		if last - y>=0:
			delt = orasP[-1]-orasP[-2]
			inc = int(y/delt)
			i = inc
			while abs(abs(last-orasP[-i])/y-1.0)>1.0e-2:
				if last-orasP[-i]<y:
					i = i + 1
				else:
					i = i - 1
			return concP[-inc-1][SiP.index(actualSp)]
		else:
			try:
			#if True:
				if len(orasP)>1:
					delt = orasP[-1] - orasP[-2]
				else:
					delt = 0
				if actualSp in globals2.modified:
					for yh in range(len(globals2.modified[actualSp])):
						spvar, pfunc = globals2.modified[actualSp][yh]
						spvar = [c.strip() for c in spvar]
						if "t" in spvar or "time" in spvar:
							try:
								return pfunc(last+delt-y)
							except:
								sprv = []
								for c in range(len(spvar)):
									sprv.append(concP2[spvar[c].strip()])
								return pfunc(*sprv)
				else:
					return 0
			except Exception as err:
				return concP[0][SiP.index(actualSp)]
	else:
		return concP[SiP.index(actualSp)]

def apply_rules(conc, yconc, oras=[], spconc =[], Si = []):	
	global orasP, concP, SiP, actualSp, concP2, RateSp 
	orasP = oras
	concP = spconc
	concP2 = conc
	SiP = Si
	tuples_pfunct = []
	for x in globals2.modified:
		for y in range(len(globals2.modified[x])):
			spvar, pfunc = globals2.modified[x][y]
			actualSp = spvar[0].strip()
			RateSp = spvar[-1].strip()
			try:
			#if True:
				sprv = []
				for c in range(len(spvar)):
					sprv.append(conc[spvar[c].strip()])
				suby = pfunc(*sprv)
				try:
					suby = sympy.re(suby.evalf())
				except:
					pass
				if not type(suby) == tuple:
					if suby != None:
						yconc[x] = suby
				else:
					tuples_pfunct.append([suby[1],x,suby[0],spvar,pfunc,len(suby)])
			#"""
			except Exception as err:
				suby = pfunc()
				try:
					suby = sympy.re(suby.evalf())
				except:
					pass
				if not type(suby) == tuple:
					if suby != None:
						yconc[x] = suby
				else:
					tuples_pfunct.append([suby[1],x,suby[0],spvar,pfunc,len(suby)])
			#"""
				
	current_key = None
	update_sp = []
	Perst = True
	for val in tuples_pfunct:
		if val[5] ==3:
			Perst = False
			
	while len(tuples_pfunct)>0:
		tuples_pfunct.sort(key = lambda a: a[0]) 
		last = 0
		for ih in range(len(tuples_pfunct)):
			if tuples_pfunct[last] == tuples_pfunct[ih]:
				pass
			else:
				bslice = tuples_pfunct[last:ih]
				rm.shuffle(bslice)
				tuples_pfunct[last:ih] = bslice
				last = ih	
		
		row = tuples_pfunct.pop()
		if current_key == None:
			current_key = row[0]
		if row[2] != None:
			x = row[1]
			yconc[x] = row[2]
			if current_key == row[0]:
				update_sp.append(x)
			if len(tuples_pfunct) == 0 or current_key != row[0]:
				if Perst:
					for xo in update_sp:
						if xo.split("_")[0] not in reserve_events_words:
							conc[xo] = yconc[xo]
				else:
					for xo in update_sp:
						if row[5]==3:
							conc[xo] = yconc[xo]
				current_key = row[0]
				update_sp = [x]
				
				for tup in tuples_pfunct:
					spvar, pfunc = tup[3:5]
					sprv = []
					for c in range(len(spvar)):
						sprv.append(conc[spvar[c].strip()])
					suby = pfunc(*sprv)
					tup[0] = suby[1]
					tup[2] = suby[0]
	
	for x in globals2.modified:
		conc[x] = yconc[x]

def get_globals(rfile):

	globals2.modified = {}
	globals2.PropModified = {}
	setattr(globals2, 'tCheck', [])

	try:
		with open(rfile,"r") as file:
			rows = []
			last = ""
			for row in file:
				if last == "Function_Definitions":
					if row.strip()!="" and row[0] != "#":
						exec(row.strip(),globals())
					elif row[0] == "#":
						last = "#"
				elif last == "#":
					if row.strip()!="" and row[0]!="@":
						rows.append(row)
					elif row.strip()!="" and row[0]!="#":
						last = "@"
				elif last == "@":
					if row.strip()!="" and row[0]!="@":
						cvar = row.split(",") 
						if len(cvar)>=3: 
							cc = ",".join(cvar[2:])
							try:
								res = re.findall("\w{1,3}\(t,.+\)\s",cc)
								if res:
									try:
										rr = float(str(res[0]).split(",")[1].replace(")",""))
									except:
										rr = eval(str(res[0]).split(",")[1].replace(")",""))						
									globals2.tCheck.append(rr)
							except:
								pass
							cc2 = cc.split(":")[0].replace("lambda","").split(",")
							
							delayPos = cc.replace("delay_","_yaled").find("delay")
							if delayPos>=0:
								dfirst = ""
								for ib in range(delayPos+6,len(cc)):
									dvar = cc[ib]
									if dvar == ",":
										break
									else:
										dfirst = dfirst + dvar
								dfirst = dfirst.strip()
								cc2new = [dfirst] + [ x for x in cc2 if x.strip() != dfirst.strip() ]
								cc2 = cc2new
								gcc = cc.split(":")[1]
								cc = "lambda "+",".join(cc2)+":"+gcc
								
							RatePos = cc.find("rateOf")
							if RatePos>=0:
								dfirst = ""
								for ib in range(RatePos+7,len(cc)):
									dvar = cc[ib]
									if dvar == ")":
										break
									else:
										dfirst = dfirst + dvar
								dfirst = dfirst.strip()
								cc2new =  [ x for x in cc2 if x.strip() != dfirst.strip() ] + [dfirst]
								cc2 = cc2new
								gcc = cc.split(":")[1]
								cc = "lambda "+",".join(cc2)+":"+gcc								
								
							if cvar[0].strip() not in globals2.modified:
								globals2.modified[cvar[0].strip()] = [[cc2,eval(cc)]]
							else:
								globals2.modified[cvar[0].strip()].append([cc2,eval(cc)])
				elif row[0] == "#":
					last = "#"
					gg = row.split(",")[1:]
					try:
						for x in gg:
							xx = x.split("=") 
							globals2.settings[xx[0].strip()] = xx[1].strip()
					except:
						pass
				elif row[0]== "@":
					last = "@"
				elif row.strip() == "Function_Definitions:":
					last = "Function_Definitions"
			file.close()
			
		Rxn = len(rows)
		for ih in range(Rxn):
			col_row = rows[ih].split(":::::")
			row = col_row[0].strip().split(",")
			if len(col_row)>1:
				krow = col_row[1].strip().split(":::")
				if len(krow)==2:
					cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
					cc3 = krow[1].split(":")[0].replace("lambda","").split(",")
					globals2.PropModified["Prop_"+str(ih)] = [(cc2,eval(krow[0])),(cc3,eval(krow[1]))]			
				else:
					cc2 = krow[0].split(":")[0].replace("lambda","").split(",")
					globals2.PropModified["Prop_"+str(ih)] = [(cc2,eval(krow[0]))]
					
		globals2.tCheck.sort()
	except:
		pass
