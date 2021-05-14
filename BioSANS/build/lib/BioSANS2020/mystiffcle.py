import numpy as np
from propensity import *
from recalculate_globals import *
import mglobals as globals2
	
def cle_model(Sp,Ks,conc,Rr,Rp,V,dt,delX,reg=False):
	D = propensity_vec(Ks,conc,Rr,Rp)
	G = np.sqrt(D)
	N = np.random.randn(len(D)).reshape(len(D),1)
	fx = np.matmul(V,D)
	gx = np.matmul(V,G*N)  
	if not reg:
		h1 = np.abs(fx)+1.0e-30
		h2 = np.abs(gx)+1.0e-30
		dt = min(delX*min( np.min(1/h1),np.min(1/h2) ),dt)
	sqdt = np.sqrt(dt)
	
	ind = 0
	for sp in Sp:
		key = sp.strip().split("_")[0]
		if key in reserve_events_words:	
			gx[ind] = 0
		ind = ind + 1
	
	fd = fx*dt+gx*sqdt   
	return [fd.reshape(len(Sp)),dt]
	
def cle_calculate(t,Sp,Ks,sconc,Rr,Rp,V,delX=10,rr=1,implicit=False,rfile=""):
	get_globals(rfile)
	np.random.seed(int(rr*100))
	tnew = []
	dt = t[-1]-t[-2]
	yconc = { x:sconc[x] for x in sconc }
	conc = { x:sconc[x] for x in sconc }
	apply_rules(conc, yconc)
	S = [[conc[z] for z in Sp]]
	
	if not implicit:
		tnow = t[0]
		tnew.append(tnow)
		while tnow<t[-1]:
			m = np.nan_to_num(cle_model(Sp,Ks,conc,Rr,Rp,V,dt,delX))
			mm = m[0].reshape(1,len(m[0]))[0] 
			tnow = tnow+m[1]
			ind = 0
			for sp in Sp:
				conc[sp] = S[-1][ind] + mm[ind]
				ind = ind+1
			apply_rules(conc, yconc)
			S.append([conc[z] for z in Sp])
			for sp in Sp:
				conc[sp] = max(0, conc[sp])
			tnew.append(tnow)
	else:
		tnow = t[0]
		tnew.append(tnow)
		tindex = 1
		C = [conc[z] for z in Sp]
		while tnew[-1]<t[-1]:
			m = np.nan_to_num(cle_model(Sp,Ks,conc,Rr,Rp,V,dt,delX))
			tnow = tnow+m[1]
			if tnow>t[tindex]:
				tnow = tnow-m[1]
				dt2 = t[tindex]-tnow
				m = np.nan_to_num(cle_model(Sp,Ks,conc,Rr,Rp,V,dt2,delX,True))
				mm = m[0].reshape(1,len(m[0]))[0] 
				tnow = tnow+m[1]
				
				ind = 0
				for sp in Sp:
					conc[sp] = C[ind] + mm[ind]
					ind = ind+1
				apply_rules(conc, yconc)
				C = [conc[z] for z in Sp]
				for sp in Sp:
					conc[sp] = max(0, conc[sp])
				S.append(C)
					
					
				tnew.append(tnow)
				tindex = tindex + 1
			else:
				mm = m[0].reshape(1,len(m[0]))[0] 
				ind = 0
				for sp in Sp:
					conc[sp] = C[ind] + mm[ind]
					ind = ind+1
				apply_rules(conc, yconc)
				C = [conc[z] for z in Sp]
				for sp in Sp:
					conc[sp] = max(0, conc[sp])
	return (tnew, np.array(S))
	
def cle2_calculate(t,Sp,Ks,sconc,Rr,Rp,V,delX=1,rr=1,rfile=""):
	get_globals(rfile)
	np.random.seed(int(rr*100))
	div = max(1,int(1/delX))
	dt = (t[-1]-t[-2])/div
	yconc = { x:sconc[x] for x in sconc }
	conc = { x:sconc[x] for x in sconc }
	apply_rules(conc, yconc)
	S = [[conc[z] for z in Sp]]	
	tnow = t[0]
	tnew = [tnow]
	dv = 0
	
	while abs(tnow-t[-1])>1.0e-10:
		m = np.nan_to_num(cle_model(Sp,Ks,conc,Rr,Rp,V,dt,delX,True))
		mm = m[0].reshape(1,len(m[0]))[0] 
		tnow = tnow+m[1]
		ind = 0
		for sp in Sp:
			conc[sp] = S[-1][ind] + mm[ind]
			ind = ind+1
		apply_rules(conc, yconc)
		if dv == div-1:
			S.append([conc[z] for z in Sp])
			tnew.append(tnow)
			dv = 0
		else:
			dv = dv + 1
		for sp in Sp:
			conc[sp] = max(0, conc[sp])
		
	return (tnew, np.array(S))