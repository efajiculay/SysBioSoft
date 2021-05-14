from scipy.integrate import odeint
from propensity import *
from recalculate_globals import *
from LNAapprox import lna_ss_jacobian, LNA_model_ss
import mglobals as globals2
import numpy as np

def LNA_ode_model(z,t,Sp,Ks,Rr,Rp,V,molar=False):  
	Spc = [s for s in Sp]
	conc = { Spc[x]: z[x] for x in range(len(Spc)) }
	if not molar:
		D = propensity_vec(Ks,conc,Rr,Rp,True)
	else:
		D = propensity_vec_molar(Ks,conc,Rr,Rp,True)

	dxdt = np.matmul(V,D).reshape(len(z))
	for x in globals2.conBoundary:
		ind = Spc.index(x)
		dxdt[ind] = 0
	return dxdt 
	
def LNA_cov_model(A,B,C):
	return np.matmul(A,C)+np.matmul(C,np.transpose(A))+B

def LNA_non_steady_state_old(conc,t,Sp,Ks,Rr,Rp,V,molar=True,rfile="",delX=10): 
	get_globals(rfile)
	z = [conc[a] for a in Sp]
	zoftime = odeint(LNA_ode_model,z,t, args=(Sp,Ks,Rr,Rp,V,molar)) 
	
	half = []
	SiNew = []
	Sps = [ s for s in Sp ]
	lenSps = len(Sps)
	for i in range(lenSps):
		s = Sps[i]
		for j in range(i,lenSps):
			p = Sps[j]
			SiNew.append("cov("+s+","+p+")")		
			half.append(lenSps*i+j)

	C = np.zeros((len(z),len(z)))
	Cov = [[x for x in C.flatten()[half]]]
	dt = t[-1]-t[-2]
	tnew = [0]
	
	for S in zoftime[1:]:
		ind = 0
		for sp in Sp:
			conc[sp] = S[ind]
			ind = ind+1		
		AA = lna_ss_jacobian(LNA_model_ss, S, Sp, V, Ks, Rr, Rp)
		f = propensity_vec_molar(Ks,conc,Rr,Rp,True)
		BB = np.matmul(np.matmul(V,np.diag(f.flatten())),V.T)
		
		A = np.nan_to_num(AA)
		B = np.nan_to_num(BB)
		
		fx = LNA_cov_model(A,B,C)
		C = C + fx*dt
		Cov.append([x for x in C.flatten()[half]])
		
	return [np.array(Cov),SiNew,t]
	
def LNA_non_steady_state(conc,t,Sp,Ks,Rr,Rp,V,molar=True,rfile="",delX=10): 
	get_globals(rfile)
	z = [conc[a] for a in Sp]
	#zoftime = odeint(LNA_ode_model,z,t, args=(Sp,Ks,Rr,Rp,V,molar)) 
	
	half = []
	SiNew = []
	Sps = [ s for s in Sp ]
	lenSps = len(Sps)
	for i in range(lenSps):
		s = Sps[i]
		for j in range(i,lenSps):
			p = Sps[j]
			SiNew.append("cov("+s+","+p+")")		
			half.append(lenSps*i+j)

	C = np.zeros((len(z),len(z)))
	Cov = [[x for x in C.flatten()[half]]]
	dt = t[-1]-t[-2]
	tnew = [0]
	
	#for S in zoftime[1:]:
	S = np.array(z)
	while tnew[-1]<t[-1]:
		ind = 0
		for sp in Sp:
			conc[sp] = S[ind]
			ind = ind+1		
		AA = lna_ss_jacobian(LNA_model_ss, S, Sp, V, Ks, Rr, Rp)
		f = propensity_vec_molar(Ks,conc,Rr,Rp,True)
		BB = np.matmul(np.matmul(V,np.diag(f.flatten())),V.T)
		
		A = np.nan_to_num(AA)
		B = np.nan_to_num(BB)
		
		fxCovr = LNA_cov_model(A,B,C)
		fxConc = LNA_ode_model(S,t,Sp,Ks,Rr,Rp,V,molar)
		
		h1 = np.abs(fxCovr)+1.0e-30
		dt = max(min(delX*np.min(1/h1),dt),1.0e-4)	
		tnew.append(tnew[-1]+dt)	
		S = S + fxConc*dt
		C = C + fxCovr*dt
		Cov.append([x for x in C.flatten()[half]])
		
	return [np.array(Cov),SiNew,tnew]
	
def LNA_non_steady_state2(conc,t,Sp,Ks,Rr,Rp,V,molar=True,rfile="",delX=10): 
	z = [conc[a] for a in Sp]
	Cov, Si, tnew = LNA_non_steady_state(conc,t,Sp,Ks,Rr,Rp,V,molar,rfile,delX)
	zoftime = odeint(LNA_ode_model,z,tnew, args=(Sp,Ks,Rr,Rp,V,molar))		

	SiNew = []
	Sps = [ s for s in Sp ]
	lenSps = len(Sps)
	FFdiv = []
	
	for x in Si:
		SiNew.append(x.replace("cov","FF"))
	
	for S in zoftime:
		row = []
		for i in range(lenSps):
			s = Sps[i]
			for j in range(i,lenSps):
				p = Sps[j]
				Sij = S[i]*S[j]
				row.append(np.sqrt( Sij if Sij != 0 else 1 ))	
		FFdiv.append(row)
	
	return [Cov/np.array(FFdiv),SiNew,tnew]
		