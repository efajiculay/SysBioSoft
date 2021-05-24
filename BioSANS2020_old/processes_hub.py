#import sys
#import os
#sys.path.append(os.path.abspath("BioSANS2020"))

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from matplotlib.cm import get_cmap
import matplotlib._color_data as mcd
from BioSANS2020.propensity import *
from BioSANS2020.mystiffcle import *
from BioSANS2020.Gillespie import *
from BioSANS2020.Tau_leaping import *
from BioSANS2020.Tau_leaping2 import *
from BioSANS2020.mytauleap import *
from BioSANS2020.Euler import *
from BioSANS2020.LNAapprox import *
from BioSANS2020.LNAapprox2 import *
from BioSANS2020.NetworkLoc import *
from BioSANS2020.ode_int import *
from BioSANS2020.LNAfunctionOfTime import *
from BioSANS2020.sde_int import *
from BioSANS2020.runge_kutta4 import *
from BioSANS2020.param_est import *
from BioSANS2020.Analytical_sol import *
from BioSANS2020.create_wxmaxima_command import *
from BioSANS2020.Conversion import *

from BioSANS2020.draw_figure import *
import BioSANS2020.mglobals as globals2
import BioSANS2020.proc_global as proc_global
from BioSANS2020.param_slider import *

def process_hub(
	t,Sp,Ksn,concn,Rr,Rp,V,Vm=1,miter=5,logx=False,logy=False,
	delX=10,normalize=False,method="CLE",mix_plot=True,save=True,
	out_fname="",plot_show=True,vary="",mult_proc=False,items=None,
	vary2="",implicit=False,rfile="",expDataFile=None
):

	Vv = []
	nz = []
	SiNew = []
	Si = [x for x in Sp]
	for row in range(V.shape[0]):
		if np.sum(np.abs(V[row,:])) != 0 and Si[row][0]!="-":
			Vv.append(list(V[row,:]))
			nz.append(row)	   
	V = np.array(Vv)
	Sp = {Si[z]:Sp[Si[z]] for z in nz}

	Si=[Si for Si in Sp]
	data = []
	if len(vary)>0:
		hold = vary.pop()
		miter = len(vary)
	if len(vary2)>0:
		r1,r2 = vary2[0]
		kval  = vary2[1]    
		miter = len(kval)
	if mult_proc:
		pool = multiprocessing.Pool(min(miter,round(0.9*multiprocessing.cpu_count())))
		CONCN = []
		KSNS = []
		if len(vary)>0:
			for j in range(miter):
				concn[hold] = vary[j]   
				CONCN.append({x : concn[x] for x in concn})
		else:
			CONCN = [concn for x in range(miter)]				
		if len(vary2)>0:
			for j in range(miter):  
				KSNS.append({x : Ksn[x][:] for x in Ksn})
				if type(r1) == int:
					KSNS[-1][r1][r2] = kval[j]				 
				else:
					for k in range(len(r1)):   
						KSNS[-1][r1[k]][r2[k]] = kval[k][j]
		else:
			KSNS = [Ksn for x in range(miter)]			 
			
#		if __name__ == '__main__':
		if 1 == 1:	 #always true, just use above  command on some OS
			rands = [ x*np.random.rand() for x in range(miter) ]
			if method == "Tau-leaping":
				results = [
					pool.apply_async( 
						Tau_leaping, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,rands[ih],delX,implicit,rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + (Si,) for result in results ]
			elif method == "Tau-leaping2":
				results = [
					pool.apply_async( 
						Tau_leaping2, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,rands[ih],delX,implicit,rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + (Si,) for result in results ]
			elif method == "Sim-TauLeap":
				results = [
					pool.apply_async( 
						Sim_TauLeap, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,rands[ih],delX,implicit,rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + (Si,) for result in results ]				
			elif method=="CLE":
				results = [
					pool.apply_async( 
						cle_calculate, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,delX,rands[ih],implicit,rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + (Si,) for result in results ]
			elif method=="CLE2":
				results = [
					pool.apply_async( 
						cle2_calculate, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,delX,rands[ih],rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + (Si,) for result in results ]
			elif method=="Euler-1":
				results = [
					pool.apply_async( 
						Euler_int, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,delX,False,None,implicit,False,rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + (Si,) for result in results ]	
			elif method in ["Euler-2","Euler-3"]:
				results = [
					pool.apply_async( 
						Euler_int, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,delX,False,None,implicit,True,rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + (Si,) for result in results ]				
			elif method=="LNA":
				save = False
				plot_show = False
				results = [
					pool.apply_async( 
						Euler_int, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,delX,True,items,False,rfile)
					) for ih in range(miter)
				]
				data = [ result.get() + [Si]  for result in results ]
			elif method=="Gillespie_":
				results = [
					pool.apply_async( 
						Gillespie, 
						args = (t,Sp,KSNS[ih],CONCN[ih],Rr,Rp,V,rands[ih],implicit,rfile)
					) for ih in range(miter)
				]				
				data = [ result.get() + (Si,) for result in results ]
			else:
				print("Multiprocessing not supported for your method of choice")
		pool.close()
	else:
		rands = [ x*np.random.rand() for x in range(miter) ]
		for j in range(miter):	  
			tnew = []
			rr = rands[j]
			if len(vary)>0:
				concn[hold] = vary[j]		
			if len(vary2)>0:
				if type(r1)==int:
					Ksn[r1][r2] = kval[j]				
				else:
					for k in range(len(r1)):
						Ksn[r1[k]][r2[k]] = kval[k][j]
			if method=="CLE":
				tnew, z = cle_calculate(t,Sp,Ksn,concn,Rr,Rp,V,delX,rr,implicit,rfile)
			elif method=="CLE2":
				tnew, z = cle2_calculate(t,Sp,Ksn,concn,Rr,Rp,V,delX,rr,rfile)
			elif method=="Gillespie_":
				tnew, z = Gillespie(t,Sp,Ksn,concn,Rr,Rp,V,rr,implicit,rfile)
			elif method == "Tau-leaping":
				tnew, z = Tau_leaping(t,Sp,Ksn,concn,Rr,Rp,V,rr,delX,implicit,rfile)
			elif method == "Tau-leaping2":
				tnew, z = Tau_leaping2(t,Sp,Ksn,concn,Rr,Rp,V,rr,delX,implicit,rfile)
			elif method == "Sim-TauLeap":
				tnew, z = Sim_TauLeap(t,Sp,Ksn,concn,Rr,Rp,V,rr,delX,implicit,rfile)
			elif method=="Euler-1":
				tnew, z = Euler_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,False,None,implicit,False,rfile)
			elif method in ["Euler-2","Euler-3"]:
				tnew, z = Euler_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,False,None,implicit,True,rfile)
			elif method=="Euler2-1":
				tnew, z = Euler2_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,False,None,implicit,False,rfile)
			elif method in ["Euler2-2","Euler2-3"]:
				tnew, z = Euler2_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,False,None,implicit,True,rfile)
			elif method=="LNA":
				plot_show = False
				save = False			
				tnew, z = Euler_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,LNAsolve=True,items=items,rfile=rfile)  
				tnew, z = Euler2_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,False,None,implicit,True,rfile)
			elif method=="LNA(t)":
				z, SiNew, tnew = LNA_non_steady_state(concn,t,Sp,Ksn,Rr,Rp,V,molar=True,rfile=rfile,delX=delX)
			elif method=="LNA2(t)":
				z, SiNew, tnew = LNA_non_steady_state2(concn,t,Sp,Ksn,Rr,Rp,V,molar=True,rfile=rfile,delX=delX)
			elif method=="LNA-vs":
				plot_show = False
				save = False			
				tnew, z = LNA_symbolic(Sp,Ksn,concn,Rr,Rp,V,items=items,molar=True,mode="Numeric")	
			elif method=="LNA-ks":
				plot_show = False
				save = False			
				tnew, z = LNA_symbolic(Sp,Ksn,concn,Rr,Rp,V,items=items,molar=True,mode="fofks")	
			elif method=="LNA-xo":
				plot_show = False
				save = False			
				tnew, z = LNA_symbolic(Sp,Ksn,concn,Rr,Rp,V,items=items,molar=True,mode="fofCo")					
			elif method=="LNA2":
				plot_show = False
				save = False			
				tnew, z = LNA_symbolic(Sp,Ksn,concn,Rr,Rp,V,items=items)
			elif method=="LNA3":
				plot_show = False
				save = False			
				tnew, z = LNA_symbolic(Sp,Ksn,concn,Rr,Rp,V,items=items,molar=True)				
			elif method=="NetLoc1":
				plot_show = False
				save = False
				tnew, z = NetLoc_symbolic(Sp,Ksn,concn,Rr,Rp,V,items=items,molar=True)	
			elif method=="NetLoc2":
				plot_show = False
				save = False
				tnew, z = NetLoc_symbolic(Sp,Ksn,concn,Rr,Rp,V,items=items,molar=True,numer=True)				
			elif method=="ODE-1":
				z = ODE_int(concn,t,Sp,Ksn,Rr,Rp,V,False,rfile)  
				tnew = t
			elif method in ["ODE-2", "ODE-3"]:
				z = ODE_int(concn,t,Sp,Ksn,Rr,Rp,V,True,rfile)  
				tnew = t
			elif method=="Itoint-1":
				z = SDE_int(concn,t,Sp,Ksn,Rr,Rp,V)  
				tnew = t
			elif method in ["Itoint-2", "Itoint-3"]:
				z = SDE_int(concn,t,Sp,Ksn,Rr,Rp,V,False)  
				tnew = t				
			elif method=="Stratint-1":
				z = SDE_int(concn,t,Sp,Ksn,Rr,Rp,V,True,False)  
				tnew = t
			elif method in ["Stratint-2", "Stratint-3"]:
				z = SDE_int(concn,t,Sp,Ksn,Rr,Rp,V,False,False)  
				tnew = t								
			elif method=="rk4-1":
				tnew, z = rungek4_int(concn,t,Sp,Ksn,Rr,Rp,V,False,delX,rfile)  		
			elif method in ["rk4-2","rk4-3"]:
				tnew, z = rungek4_int(concn,t,Sp,Ksn,Rr,Rp,V,True,delX,rfile)  		
			elif method=="rk4-1a":
				tnew, z = rungek4a_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,False,implicit,rfile)				
			elif method in ["rk4-2a","rk4-3a"]:
				tnew, z = rungek4a_int(t,Sp,Ksn,concn,Rr,Rp,V,delX,True,implicit,rfile)		
			elif method=="k_est1": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=True,TrueDataFil=expDataFile) 				
			elif method=="k_est2": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,TrueDataFil=expDataFile) 
			elif method=="k_est3": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=True,mode="DEvol",TrueDataFil=expDataFile) 
			elif method=="k_est4": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=False,mode="DEvol",TrueDataFil=expDataFile) 
			elif method=="k_est5": 
				setP = [logx, logy, normalize, items]
				z = ParamODE_int(concn,t,Sp,Ksn,Rr,Rp,V,True,rfile,setP)  
				plot_show = False
				save = False	
				tnew = t
			elif method=="k_est6": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=True,mode="NeldMead",TrueDataFil=expDataFile) 
			elif method=="k_est7": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=False,mode="NeldMead",TrueDataFil=expDataFile) 
			elif method=="k_est8": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=True,mode="Powell",TrueDataFil=expDataFile) 
			elif method=="k_est9": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=False,mode="Powell",TrueDataFil=expDataFile) 
			elif method=="k_est10": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=True,mode="L-BFGS-B",TrueDataFil=expDataFile) 
			elif method=="k_est11": 
				plot_show = False
				save = False				
				tnew, z = param_estimate(concn,t,Sp,Ksn,Rr,Rp,V,items=items,molar=False,mode="L-BFGS-B",TrueDataFil=expDataFile) 
			elif method=="Analyt": 
				plot_show = False
				save = False
				z = Analyt_soln(Sp,Ksn,concn,Rr,Rp,V,items=items,rfile=rfile)	
			elif method=="Analyt-ftx": 
				plot_show = False
				save = False
				z = Analyt_soln(Sp,Ksn,concn,Rr,Rp,V,items=items,rfile=rfile,mode="ftxo")
			elif method=="SAnalyt": 
				plot_show = False
				save = False	
				z = Analyt_soln(Sp,Ksn,concn,Rr,Rp,V,items=items,rfile=rfile,not_semi=False)
			elif method=="SAnalyt-ftk": 
				plot_show = False
				save = False	
				z = Analyt_soln(Sp,Ksn,concn,Rr,Rp,V,items=items,rfile=rfile,not_semi=False,mode="ftks")
			elif method=="Analyt2": 
				plot_show = False
				save = False				
				tnew, z = for_wxmaxima(Sp,Ksn,concn,Rr,Rp,V,items=items,rfile=rfile)	
			elif method=="topoTosbml":
				plot_show = False
				save = False				
				tnew, z = topo_to_sbml(Sp,Ksn,concn,Rr,Rp,V,Vm,items=items,molar=False)		
			elif method=="topoTosbml2":
				plot_show = False
				save = False				
				tnew, z = topo_to_sbml(Sp,Ksn,concn,Rr,Rp,V,Vm,items=items,molar=True)
			elif method=="topoTosbml3":
				plot_show = False
				save = False				
				tnew, z = topo_to_sbml(Sp,Ksn,concn,Rr,Rp,V,Vm,items=items,molar=None)
			data.append([tnew,z,Si])	

	nz = []
	reserve_events_words = {"t","time","status","status2","timer","finish","delay","dtime"}
	for row in range(V.shape[0]):
		key = Si[row].strip().split("_")[0]
		if key not in reserve_events_words:
			nz.append(row)	   
	Sp = {Si[z]:Sp[Si[z]] for z in nz}
	Si = [Si for Si in Sp]		

	if len(SiNew)>0:
		Si = SiNew
		Sp = SiNew
		nz = range(len(Si))
				            
	if save:
		fname = out_fname+"_"+method+".dat"
		file = open(fname,"w")
		file.write("time\t"+"\t".join(Si)+"\n")
		for traj in data:
			x = traj[0]
			y = traj[1]
			xlen = len(x)
			for ih in range(xlen):
				file.write(str(x[ih])+"\t"+"\t".join([str(yi) for yi in y[ih][nz]])+"\n")
		file.close()
				
	if plot_show:
		
		Splen = len(Sp)
		if Splen<=10:
			col = ['C'+str(i) for i in range(10)]
		elif Splen>10 and Splen <=40:
			col = [ x for x in get_cmap("tab20").colors  ] + [ x for x in get_cmap("tab20b").colors  ]
		else:
			col = [ name for name in mcd.CSS4_COLORS ]
		if mix_plot:
			plt.figure(figsize=(9.68,3.95))
			plt.xlabel("time (sec)")
			plt.ylabel("conc") 	
			if len(SiNew)>0:
				plt.ylabel("cov" if SiNew[0][0:2]!="FF" else "FF") 
			lines = []			
			if logx:
				plt.xscale('log')	
			if logy:
				plt.yscale('log') 	
			
			for j in range(miter):
				for i in nz:
					if normalize:
						line = plt.plot(data[j][0],data[j][1][:,i]/(np.max(data[j][1][:,i])+1.0e-30),color=col[i])
					else:
						line = plt.plot(data[j][0],data[j][1][:,i],color=col[i])
			if len(Si)<=10:
				plt.legend(Si)
			else:
				plt.legend(Si,fontsize='xx-small')
			plt.tight_layout()
			plt.savefig(out_fname+"_"+method+".jpg")
			fig = plt.gcf() 
			lines.append(line)
			globals2.plotted.append([plt.gca(),fig,lines])
			fig_canvas_agg = draw_figure(items,fig)		
			#plt.close()
		else:
			lines = []	
			Ssi = []
			for i in nz:
				Ssi.append(Si[i])
				plt.figure(figsize=(9.68,3.95))
				plt.xlabel("time (sec)")
				plt.ylabel("conc") 
				if len(SiNew)>0:
					plt.ylabel("cov" if SiNew[0][0:2]!="FF" else "FF") 
				if logx:
					plt.xscale('log')
				if logy:
					plt.yscale('log')  
				for j in range(miter):  
					if normalize:
						line = plt.plot(data[j][0],data[j][1][:,i]/(np.max(data[j][1][:,i])+1.0e-30),color=col[i])
					else:
						line = plt.plot(data[j][0],data[j][1][:,i],color=col[i])
					lines.append(line)
				plt.legend([Si[i]])
				plt.tight_layout()
				plt.savefig(out_fname+"_"+method+"_"+str(i)+".jpg")
				fig = plt.gcf() 
				globals2.plotted.append([plt.gca(),fig,lines])
				fig_canvas_agg = draw_figure(items,fig)	
				#plt.close()
	return data