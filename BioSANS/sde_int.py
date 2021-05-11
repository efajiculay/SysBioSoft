from scipy.integrate import odeint
from propensity import *
import sdeint
import mglobals as globals2

def sde_fx(z,t,Sp,Ks,Rr,Rp,V,molar=False):  
	Spc = [s for s in Sp]
	conc = { Spc[x]: z[x] for x in range(len(Spc)) }
	if not molar:
		D = propensity_vec(Ks,conc,Rr,Rp)
	else:
		D = propensity_vec_molar(Ks,conc,Rr,Rp)

	dxdt = np.matmul(V,D)
	for x in globals2.conBoundary:
		ind = Spc.index(x)
		dxdt[ind] = 0
	return dxdt.T[0]
	
def sde_gx(z,t,Sp,Ks,Rr,Rp,V,molar=False):  
	Spc = [s for s in Sp]
	conc = { Spc[x]: z[x] for x in range(len(Spc)) }
	if not molar:
		D = propensity_vec(Ks,conc,Rr,Rp)
	else:
		D = propensity_vec_molar(Ks,conc,Rr,Rp)
	G = np.sqrt(D)
	#N = np.random.randn(len(D)).reshape(len(D),1)	
	dxdt = np.matmul(V,G)
	for x in globals2.conBoundary:
		ind = Spc.index(x)
		dxdt[ind] = 0
	return dxdt
	
def SDE_int(conc,t,Sp,Ks,Rr,Rp,V,molar=False,ito=True): 
	z = [conc[a] for a in Sp]
	if ito:
		return sdeint.itoint(
			lambda z,t : sde_fx(z,t,Sp,Ks,Rr,Rp,V,molar),
			lambda z,t : sde_gx(z,t,Sp,Ks,Rr,Rp,V,molar),z,t
		)
	else:
		return sdeint.stratint(
			lambda z,t : sde_fx(z,t,Sp,Ks,Rr,Rp,V,molar),
			lambda z,t : sde_gx(z,t,Sp,Ks,Rr,Rp,V,molar),z,t
		)
