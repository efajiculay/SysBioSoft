import numpy as np
from propensity import *
from recalculate_globals import *
import mglobals as globals2

def Gillespie(t,Sp,Ks,conc,Rr,Rp,V,rr,implicit=False,rfile=""):
	get_globals(rfile)
	tmax = t[-1]
	np.random.seed(int(rr*100))
	tnew = []
	S = V.T
	AllSp = [z for z in Sp ]
	Spc = [z for z in Sp if z not in reserve_events_words]
	Spc2 = [z for z in Sp if z in reserve_events_words]
	concz = { x:conc[x] for x in conc }
	yconc = { x:conc[x] for x in conc }
	apply_rules(concz, yconc)
	UpdateSp = [ AllSp.index(z) for z in Spc  ]
	
	Z = [[concz[z] for z in Spc ]]
	tc = 0
	tnew.append(tc)	
	if not implicit:
		while tc<tmax:
			D = propensity_vec(Ks,concz,Rr,Rp)
			alp = np.sum(D)
			r1 = np.random.uniform()
			while r1 == 0:
				r1 = np.random.uniform()
			P = np.cumsum([ d/alp for d in D ])
			dt = (1/alp)*(np.log(1/r1))
			if np.isnan(dt) or np.isinf(dt):
				break
			r2 = np.random.uniform()
			for i in range(len(P)):
				if r2<=P[i]:
					Allpos = True
					for x in range(len(Spc)):
						holder = Z[-1][x]+S[i][UpdateSp[x]]
						if holder>=0 or Spc[x] in globals2.modified:
							concz[Spc[x]] = holder
						else:
							Allpos = True
							break
					if Allpos:
						for x in range(len(Spc2)):
							concz[Spc2[x]] = concz[Spc2[x]] + dt
						apply_rules(concz, yconc)
						Z.append([concz[x] for x in Spc])
						tc = tc + dt
						tnew.append(tc)
					else:
						for x in range(len(Spc)):
							concz[Spc[x]] = Z[-1][x]
					break
	else:
		Zc = []
		tindex = 0
		index = 0
		tchlen = len(globals2.tCheck)
		while tc<tmax:
			D = propensity_vec(Ks,concz,Rr,Rp)
			alp = np.sum(D)
			r1 = np.random.uniform()
			while r1 == 0:
				r1 = np.random.uniform()
			P = np.cumsum([ d/alp for d in D ])
			dt = (1/alp)*(np.log(1/r1))
			
			if index != tchlen:
				if tc + dt>=globals2.tCheck[index]:
					dt = globals2.tCheck[index] - tc
					index = index + 1
					
			if np.isnan(dt) or np.isinf(dt):
				Zc.append(Z[-1])
				while t[tindex]!=tmax:
					Zc.append(Z[-1])
					tindex = tindex + 1
				break
			r2 = np.random.uniform()
			for i in range(len(P)):
				if r2<=P[i]:
					Allpos = True
					for x in range(len(Spc)):
						holder = Z[-1][x]+S[i][UpdateSp[x]]
						if holder>=0 or Spc[x] in globals2.modified:
							concz[Spc[x]] = holder
						else:
							Allpos = False
							break
					if Allpos:
						for x in range(len(Spc2)):
							concz[Spc2[x]] = concz[Spc2[x]] + dt
							
						tc = tc + dt
						if "t" in Spc2:
							concz["t"] = tc 
						elif "time" in Spc2:
							concz["time"] = tc 
						else:
							pass						
													
						apply_rules(concz, yconc)
						Z.append([concz[x] for x in Spc])
						#tc = tc + dt
						try:
							if tc == t[tindex]:
								Zc.append(Z[-1])
								tindex = tindex + 1
							else:
								while tc>t[tindex]:
									Zc.append(Z[-2])
									tindex = tindex + 1
						except:
							pass
						tnew.append(tc)
					else:
						for x in range(len(Spc)):
							concz[Spc[x]] = Z[-1][x]
					break
		tnew = t
		Z = Zc
	return (tnew, np.array(Z))