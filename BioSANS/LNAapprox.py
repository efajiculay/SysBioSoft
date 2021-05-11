import numpy as np
from propensity import *
from scipy import linalg as LA
from scipy.optimize import fsolve
from scrollable_text import *

def rem_rowcol_zero(a):
	return a[:,~np.all(a == 0, axis=0)][~np.all(a == 0, axis=1)]

def lna_ss_jacobian(model, z, S, V, Ks, Rr, Rp):
	J = np.zeros((len(z), len(z)))
	div = 1
	new = 1000
	old = 2000
	check = 1000
	while check >1.0e-7:
		for i in range(len(z)):
			h = abs(z[i]*1.0e-8)/div
			h2 = 2*h
			z[i] = z[i] - h
			Ini = model(z, S, Ks, Rr, Rp, V)
			z[i] = z[i] + 2*h
			Out = model(z, S, Ks, Rr, Rp, V)
			z[i] = z[i] - h
			for j in range(len(Out)):
				J[j,i] = (Out[j]-Ini[j])/h2
		div = div*2
		old = new
		new = sum([J[k,k] for k in range(len(z))])
		check = abs(new-old)/abs(new)
	return np.array(J)

def LNA_model_ss(S,Sp,Ks,Rr,Rp,V):   
	conc = {}
	ind = 0
	for sp in Sp:
		conc[sp] = S[ind]
		ind = ind+1	
	D = propensity_vec_molar(Ks,conc,Rr,Rp,True)
	fx = np.matmul(V,D)
	return fx.reshape(len(Sp))

def LNA_steady_state(t,Sp,Ks,conc,Rr,Rp,V,items=None):
	ind = 0
	S = []
	for sp in Sp:
		S.append(conc[sp])	 
	S = fsolve(LNA_model_ss, tuple(S), xtol=1.0e-10, args=(Sp,Ks,Rr,Rp,V))
	ind = 0
	for sp in Sp:
		conc[sp] = S[ind]
		ind = ind+1		
	AA = lna_ss_jacobian(LNA_model_ss, S, Sp, V, Ks, Rr, Rp)
	f = propensity_vec_molar(Ks,conc,Rr,Rp,True)
	BB = np.matmul(np.matmul(V,np.diag(f.flatten())),V.T)
	
	AA = np.nan_to_num(AA)
	BB = np.nan_to_num(BB)
	CC = LA.solve_continuous_lyapunov(AA, -BB)   
	
	if items:
		text = prepare_scroll_text(items)
		fprint = lambda x: text.insert(INSERT," ".join([str(y) for y in x]))
	else:
		fprint = lambda x: print(" ".join([str(y) for y in x]),end="")
	
	fprint(["\nConcentrations\n\n"])
	ind = 0
	for sp in Sp:
		if sp[0]!="-":
			fprint([sp," = ",S[ind],"\n"])
			ind = ind+1	

	fprint(["\nCovariance\n\n"])
	i = 0
	for spi in Sp:
		j = 0
		for spj in Sp:
			if j>=i:
				val = CC[i,j]
				if str(val) not in {"None", "nan", "0.0"} and spi[0]!="-" and spj[0]!="-":			
					fprint([  " ".join(["Covr",spi,spj]).ljust(50),"=",val,"\n"  ])
			j = j+1
		i = i+1	
	fprint(["\nCorrelation\n\n"])
	i = 0
	for spi in Sp:
		j = 0
		for spj in Sp:
			if j>=i:
				val = CC[i,j]/np.sqrt(np.abs(CC[i,i]*CC[j,j]))
				if str(val) not in {"None", "nan", "0.0"}  and spi[0]!="-" and spj[0]!="-":				
					fprint([  " ".join(["Corr",spi,spj]).ljust(50),"=",val,"\n"  ])
			j = j+1
		i = i+1	

	fprint(["\nFano Factor\n\n"])
	ind = 0
	for sp in Sp:
		val = CC[ind,ind]/S[ind]
		if str(val) not in {"None", "nan", "0.0"}  and sp[0]!="-":
			fprint([  " ".join(["Fano Factor for",sp]).ljust(50),"=",val,"\n"  ])
		ind = ind+1	
			
	return [CC,S]
