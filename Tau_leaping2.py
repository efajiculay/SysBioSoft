import numpy as np
from propensity import *
from recalculate_globals import *

def Tau_leaping2(t,Sp,Ks,conc,Rr,Rp,V,rr,delX=1,implicit=False,rfile=""):
	get_globals(rfile)
	tmax = t[-1]
	np.random.seed(int(rr*100))
	tnew = []
	V2 = V**2
	S = V.T
	
	AllSp = [z for z in Sp ]
	Spc = [z for z in Sp if z not in reserve_events_words]
	Spc2 = [z for z in Sp if z in reserve_events_words]
	concz = { x:conc[x] for x in conc }
	yconc = { x:conc[x] for x in conc }
	apply_rules(concz, yconc)
	UpdateSp = [ AllSp.index(z) for z in Spc  ]
	
	Vcri = np.where(V<0)
	Vncr = np.array(list(range(len(V[0]))))
	Vncr = Vncr[~np.isin(Vncr,Vcri[1])]
	
	Z = [[concz[z] for z in Spc]]
	tc = 0
	tnew.append(tc)
			
	gi = []
	for i in range(len(Spc)):
		keys = Sp[Spc[i]]
		maxrlen = 0
		rlen = 0
		for key in keys:
			if Spc[i] in Rr[key]:
				maxrlen = max(maxrlen,len(Rr[key]))
			elif Spc[i] in Rp[key]:
				maxrlen = max(maxrlen,len(Rp[key]))
		if maxrlen == 1:
			for key in keys:
				if Spc[i] in Rr[key]:
					rlen = max(rlen,Rr[key][Spc[i]])
				elif Spc[i] in Rp[key]:
					rlen = max(rlen,Rp[key][Spc[i]])
			if rlen == 1:				
				gi.append(lambda x : 1)
			else:
				gi.append(lambda x : 2)
		elif maxrlen == 2:
			try:
				gi.append(lambda x : 2+1/(x-1))
			except:
				gi.append(lambda x : 2)
		else:
			gi.append(lambda x : 1)	

	if not implicit:			
		while tc<tmax:
			
			if len(Vcri[1]) > 0:
				Lcri = set()
				Lncr = set()
				for x in range(len(Vncr)):
					Lncr.add(x)

				for x in range(len(Vcri[0])):
					i, j = [ Vcri[0][x], Vcri[1][x] ] 
					if abs(Z[-1][i]/V[i,j])<10:
						Lcri.add(j)
					else:
						Lncr.add(j)
			else:
				Lcri = set()
				
			D = propensity_vec(Ks,concz,Rr,Rp)
			D[D<0] = 0
			if len(Lcri) == len(V[0]):
				alp = np.sum(D)		   
				dt = (1/alp)*(np.log(1/np.random.uniform())) 
			elif len(Lcri)== 0:
				alp = np.sum(D)
				uuj = np.matmul(V,D)
				sig = np.matmul(V2,D)  
				exigi = np.array( [Z[-1][j]*(1/gi[j](Z[-1][j])) for j in range(len(Spc))] )*0.03*delX
				exigi1 = np.maximum(exigi,np.full(len(Spc),1))
				dt = min(np.min(exigi1/np.abs(uuj)),np.min(exigi1*exigi1/np.abs(sig)))
			else:
				Lcri = np.array(list(Lcri))
				Lncr = np.array(list(Lncr))
				alpc = np.sum(D[Lcri])		   
				dtc = (1/alpc)*(np.log(1/np.random.uniform())) 
				
				alp = np.sum(D[Lncr])
				uuj = np.matmul(V[:,Lncr],D[Lncr])
				sig = np.matmul(V2[:,Lncr],D[Lncr])  
				exigi = np.array( [Z[-1][j]*(1/gi[j](Z[-1][j])) for j in range(len(Spc))] )*0.03*delX
				exigi1 = np.maximum(exigi,np.full(len(Spc),1))
				dt = min(np.min(exigi1/np.abs(uuj)),np.min(exigi1*exigi1/np.abs(sig)),dtc)			
				
			K = np.random.poisson(D*dt)
			Allpos = True
			cc = {}
			bb = np.sum(K*S[:,UpdateSp],0)
			for x in range(len(Spc)):
				holder = Z[-1][x]+bb[x]
				if holder>=0 or Spc[x] in globals2.modified:
					cc[Spc[x]] = holder
				else:
					Allpos = False
					break
			if Allpos:
				for x in range(len(Spc)):
					concz[Spc[x]] = cc[Spc[x]]
				for x in range(len(Spc2)):
					concz[Spc2[x]] = concz[Spc2[x]] + dt
				apply_rules(concz, yconc)
				Z.append([concz[x] for x in Spc])
				tc = tc + dt
				tnew.append(tc)	
	else:
		tindex = 1
		C = [concz[z] for z in Spc]
		while t[tindex]<tmax:
			
			if len(Vcri[1]) > 0:
				Lcri = set()
				Lncr = set()
				for x in range(len(Vncr)):
					Lncr.add(x)

				for x in range(len(Vcri[0])):
					i, j = [ Vcri[0][x], Vcri[1][x] ] 
					if abs(Z[-1][i]/V[i,j])<10:
						Lcri.add(j)
					else:
						Lncr.add(j)
			else:
				Lcri = set()
				
			D = propensity_vec(Ks,concz,Rr,Rp)
			D[D<0] = 0
			if len(Lcri) == len(V[0]):
				alp = np.sum(D)		   
				dt = (1/alp)*(np.log(1/np.random.uniform())) 
			elif len(Lcri)== 0:
				alp = np.sum(D)
				uuj = np.matmul(V,D)
				sig = np.matmul(V2,D)  
				exigi = np.array( [C[j]*(1/gi[j](C[j])) for j in range(len(Spc))] )*0.03*delX
				exigi1 = np.maximum(exigi,np.full(len(Spc),1))
				dt = min(np.min(exigi1/np.abs(uuj)),np.min(exigi1*exigi1/np.abs(sig)))
			else:
				Lcri = np.array(list(Lcri))
				Lncr = np.array(list(Lncr))
				alpc = np.sum(D[Lcri])		   
				dtc = (1/alpc)*(np.log(1/np.random.uniform())) 
				
				alp = np.sum(D[Lncr])
				uuj = np.matmul(V[:,Lncr],D[Lncr])
				sig = np.matmul(V2[:,Lncr],D[Lncr])  
				exigi = np.array( [C[j]*(1/gi[j](C[j])) for j in range(len(Spc))] )*0.03*delX
				exigi1 = np.maximum(exigi,np.full(len(Spc),1))
				dt = min(np.min(exigi1/np.abs(uuj)),np.min(exigi1*exigi1/np.abs(sig)),dtc)			
			if tc + dt>t[tindex]:		
				dt = max(0,t[tindex]-tc)
				K = np.round(np.random.poisson(D*dt))
				Allpos = True
				cc = {}
				bb = np.sum(K*S[:,UpdateSp],0)
				for x in range(len(Spc)):
					holder = C[x]+bb[x]
					if holder>=0 or Spc[x] in globals2.modified:
						cc[Spc[x]] = holder
					else:
						Allpos = False
						break
				if Allpos:
					for x in range(len(Spc)):
						concz[Spc[x]] = cc[Spc[x]]
					for x in range(len(Spc2)):
						concz[Spc2[x]] = concz[Spc2[x]] + dt
					apply_rules(concz, yconc)
					Z.append([concz[x] for x in Spc])
					C = Z[-1]
					tc = tc + dt
					tindex = tindex + 1
			else:
				dt = max(0,dt)
				K = np.round(np.random.poisson(D*dt))
				Allpos = True
				cc = {}
				bb = np.sum(K*S[:,UpdateSp],0)
				for x in range(len(Spc)):
					holder = C[x]+bb[x]
					if holder>=0 or Spc[x] in globals2.modified:
						cc[Spc[x]] = holder
					else:
						Allpos = False
						break
				if Allpos:
					for x in range(len(Spc)):
						concz[Spc[x]] = cc[Spc[x]]
					for x in range(len(Spc2)):
						concz[Spc2[x]] = concz[Spc2[x]] + dt
					apply_rules(concz, yconc)
					C = [ concz[x] for x in Spc ]
					tc = tc + dt	
					
		if len(t)!= len(Z):
			while len(t)!= len(Z):
				Z.append(C)
		tnew = t
	return (tnew,np.array(Z))