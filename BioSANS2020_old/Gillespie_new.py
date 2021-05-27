#import sys
#import os
#sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
from BioSANS2020.propensity import *
from BioSANS2020.recalculate_globals import *
import BioSANS2020.mglobals as globals2

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
			#dtOri = 0
			if index != tchlen:
				if tc + dt>=globals2.tCheck[index]:
					#dtOri = dt
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

						try:
							while tc>=t[tindex]:
								Zc.append(Z[-2]) if tc>t[tindex] else Zc.append(Z[-1])
								tindex = tindex + 1
						except:
							pass
						#if dtOri>0:
							#tc = tc - dt + dtOri
						tnew.append(tc)
					else:
						for x in range(len(Spc)):
							concz[Spc[x]] = Z[-1][x]
					break
		tnew = t
		Z = Zc
	return (tnew, np.array(Z))