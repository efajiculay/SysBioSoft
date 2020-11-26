import mglobals as globals2
from sbmlMath import *
import re
from sympy import *

reserve_events_words = {"t","time","status","status2","timer","finish","delay","dtime"}

def apply_rules(conc, yconc):	
	tuples_pfunct = []
	for x in globals2.modified:
		for y in range(len(globals2.modified[x])):
			spvar, pfunc = globals2.modified[x][y]
			try:
				sprv = []
				for c in range(len(spvar)):
					sprv.append(conc[spvar[c].strip()])
				suby = pfunc(*sprv)
				if not type(suby) == tuple:
					if suby != None:
						yconc[x] = suby
				else:
					tuples_pfunct.append([suby[1],x,suby[0],spvar,pfunc,len(suby)])
			except:
				suby = pfunc()
				if not type(suby) == tuple:
					if suby != None:
						yconc[x] = suby
				else:
					tuples_pfunct.append([suby[1],x,suby[0],spvar,pfunc,len(suby)])
				
	current_key = None
	update_sp = []
	Perst = True
	for val in tuples_pfunct:
		if val[5] ==3:
			Perst = False
			
	while len(tuples_pfunct)>0:
		tuples_pfunct.sort(key = lambda a: a[0]) 
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