import sys
import os
sys.path.append(os.path.abspath("BioSANS2020"))

from scipy.integrate import odeint
from BioSANS2020.propensity import *
from BioSANS2020.recalculate_globals import *
import BioSANS2020.mglobals as globals2

def ODE_model(z,t,Sp,Ks,Rr,Rp,V,molar=False):  
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
	
def ODE_int(conc,t,Sp,Ks,Rr,Rp,V,molar=False,rfile=""): 
	get_globals(rfile)
	z = [conc[a] for a in Sp]
	return odeint(ODE_model,z,t, args=(Sp,Ks,Rr,Rp,V,molar))  