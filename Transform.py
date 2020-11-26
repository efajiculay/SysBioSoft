from sample_points import *
from scrollable_text import *
from draw_figure import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np
import mglobals as globals

def normalize(v):
	norm = np.linalg.norm(v)
	if norm == 0: 
		return v
	return v / norm

def calc_cross_corr(edata,items):
	if len(edata[0])!=3:
		ndata, Si = sample_points(edata)
		t = ndata[0]
		N = len(ndata)
		dlen = len(Si)
		
		ddm = [0]*len(Si)
		for i in range(1,N):
			dd = ndata[i]
			for j in range(dlen):
				ddm[j] = dd[j] + ddm[j]
				
		ddm = [ ddm[j]/(N-1) for j in range(dlen)]
		vvv = [ [ [0] for x in range(dlen) ] for x in range(dlen) ]
		for i in range(1,N):
			dd = ndata[i]
			for j in range(dlen):
				for k in range(dlen):
					vvv[j][k] = vvv[j][k]+np.correlate(normalize(dd[j][-1000:]-ddm[j][-1000:]),normalize(dd[k][-1000:]-ddm[k][-1000:]),"full")
		vvv = np.array(vvv)/(N-1)	
		for j in range(dlen):
			for k in range(dlen):
				plt.figure(figsize=(9.68,3.95))
				plt.xlabel("lag")	
				plt.ylabel(Si[j]+"-"+Si[k]+" correlation")
				x = int(len(vvv[j,k])/2)
				x = [ -y for y in list(range(x)) ][::-1]+[0]+list(range(x))
				line=plt.plot(x,vvv[j,k])
				fig = plt.gcf()
				plt.tight_layout()
				globals.plotted.append([plt.gca(),fig,[line]])
				fig_canvas_agg = draw_figure(items,fig)		
				plt.close()

def calc_covariance2(edata):
	ndata, Si = sample_points(edata)
	t = ndata[0]
	N = len(ndata)
	dlen = len(Si)
	
	ddm = [0]*len(Si)
	for i in range(1,N):
		dd = np.array(ndata[i])[:,-100:]
		for j in range(dlen):
			ddm[j] = dd[j] + ddm[j]
			
	ddm = np.mean(np.array([ ddm[j]/(N-1) for j in range(dlen)]),1).reshape(dlen,1)
	
	vv = 0
	for i in range(1,N):
		dd = np.array(ndata[i])[:,-100:]-ddm
		ddd = np.zeros((dlen,dlen))
		for j in range(dlen):
			for k in range(dlen):
				ddd[j,k] = np.mean(dd[j]*dd[k])
		vv = vv+ddd	
	vv = vv/(N-1)
	print("covariance")
	print(Si)
	print(vv)
	print()
	vvv = np.zeros((dlen,dlen))
	for i in range(dlen):
		for j in range(dlen):
			vvv[i,j] = vv[i,j]/np.sqrt(vv[i,i]*vv[j,j])
	print("\ncorrelation")
	print(Si)
	print(vvv)		
					
def calc_covariance(edata,items,points=100):
	data, Si = edata
	N = len(data)
	ddm = 0
	for i in range(1,N):
		dd = data[i][-points:,1:]
		ddm = ddm + dd
		
	ddm = np.mean(ddm/(N-1),0)
	print(ddm)
	
	vv = 0
	dlen = len(Si)
	for k in range(1,N):
		dd = data[k][-points:,1:] - ddm
		ddd = np.zeros((dlen,dlen))
		for i in range(dlen):
			for j in range(dlen):
				ddd[i,j] = np.mean(dd[:,i]*dd[:,j])
		vv = vv+ddd	
	vv = vv/(N-1)
	
	text = prepare_scroll_text(items)
	text.insert(INSERT,"covariance\n\n")
	for i in range(dlen):
		for j in range(i,dlen):	
			if str(vv[i,j]) not in {"None", "nan", "0.0"}:
				label = Si[i]+" and "+Si[j]
				text.insert(INSERT,label.ljust(50)+" = "+str(vv[i,j])+"\n")
				
	vvv = np.zeros((dlen,dlen))
	text.insert(INSERT,"\ncorrelation\n\n")
	for i in range(dlen):
		for j in range(i,dlen):
			if str(vv[i,j]) not in {"None", "nan", "0.0"}:
				vvv[i,j] = vv[i,j]/np.sqrt(vv[i,i]*vv[j,j])
				label = Si[i]+" and "+Si[j]
				text.insert(INSERT,label.ljust(50)+" = "+str(vvv[i,j])+"\n")				
				
	dd = 0
	for i in range(N):
		dd = dd+(data[i][-points:,1:]-ddm)**2
	FF = np.mean(dd/(N-1),0)/ddm				
	text.insert(INSERT,"\nFano Factor\n\n")
	for i in range(dlen):
		if str(FF[i]) not in {"None", "nan"}:
			text.insert(INSERT,Si[i].ljust(50)+" = "+str(FF[i])+"\n")				
			
	globals.plot_i = globals.plot_i+1	
			
def fano_factor(edata,items,points=100):
	data, Si = edata
	N = len(data)
	ddm = 0
	for i in range(N):
		dd = data[i][-points:,1:]
		ddm = ddm + dd
		
	ddm = np.mean(ddm/(N-1),0)
	
	dd = 0
	for i in range(N):
		dd = dd+(data[i][-points:,1:]-ddm)**2
	FF = np.mean(dd/(N-1),0)/ddm
	
	text = prepare_scroll_text(items)
	dlen = len(Si)
	text.insert(INSERT,"Fano Factor\n\n")
	for i in range(dlen):
		if str(FF[i]) not in {"None", "nan"}:
			text.insert(INSERT,Si[i].ljust(50)+" = "+str(FF[i])+"\n")
	
	globals.plot_i = globals.plot_i+1	
	
def Prob_density(edata,items):
	data, Si = edata
	N = len(data)
	dlen = len(Si)
	
	dd = data[0][:,1:]
	for i in range(1,N):
		dd =  np.concatenate((dd, data[i][:,1:]), axis=0)
		
	for j in range(dlen):
		plt.figure(figsize=(9.68,3.95))
		plt.xlabel("conc("+Si[j]+")")
		plt.ylabel("Prob") 		
		line = plt.hist(dd[:,j],bins=100, density=True, orientation='vertical')
		plt.legend([Si[j]])
		plt.tight_layout()
		fig = plt.gcf()		
		globals.plotted.append([plt.gca(),fig,[line]])
		fig_canvas_agg = draw_figure(items,fig)			
		plt.close()
		
def Prob_density2(edata,items):
	data, Si = edata
	N = len(data)
	dlen = len(Si)
	
	dd = data[0][:,1:]
	tt = data[0][:,0].flatten()
	tlen = len(tt)
	for i in range(1,N):
		dd =  np.concatenate((dd, data[i][:,1:]), axis=0)
		tt = np.concatenate((tt, data[i][:,0]))
			
	for j in range(dlen):
		r = dd
		Zs, Xs, Ys = np.histogram2d(tt, dd[:,j].flatten(), bins=(100, 100))
		plt.figure(figsize=(9.68,3.95))
		plt.xlabel("time")
		plt.ylabel("conc("+Si[j]+")")
		cntr = plt.contour(Xs[1:],Ys[1:],Zs.T,list(range(int(np.min(Zs)),int(np.max(Zs))+1,2)),cmap="Spectral_r")	
		fig = plt.gcf()
		fig.colorbar(cntr)
		plt.tight_layout()
		globals.plotted.append([plt.gca(),fig,[cntr]])
		fig_canvas_agg = draw_figure(items,fig)		
		plt.close()
		
def Prob_density3(edata,items,bins=50):
	if len(edata[0])!=3:
		ndata, Si = sample_points(edata)
		t = ndata[0]
		N = len(ndata)
		dlen = len(Si)
		
		ddd = [ [] for x in range(dlen) ]
		for k in range(1,N):
			dd = ndata[k]
			for j in range(dlen):
				ddd[j].append(dd[j])
		ddd = np.array(ddd)

		vvv = {}
		for i in range(dlen):
			for j in range(len(t)):
				vvv[(i,j)] = np.histogram(ddd[i][:,j],bins=bins)
		
		tlen = len(t)
		strt = int(tlen/10)
		for i in range(dlen):
			lines = []
			fig = plt.figure(figsize=(9.68,3.95))
			ax = fig.gca(projection='3d')
			for j in range(strt,tlen,int(tlen/5)):
				xs = vvv[(i,j)][1][1:]
				ys = vvv[(i,j)][0]
				width = (xs[1]-xs[0])
				line = ax.bar(xs, ys, zs=t[j], zdir='x', width=max(width,0.8), alpha=1)
				lines.append(line)
			ax.set_xlabel('time')
			ax.set_ylabel("conc("+Si[i]+")")
			ax.set_zlabel('freq')
			ax.view_init(elev=40, azim=-120)
			globals.plotted.append([ax,fig,lines])
			fig_canvas_agg = draw_figure(items,fig)	
			plt.close()			
			
def Ave_traj(edata,items):
	if len(edata[0])!=3:
		ndata, Si = sample_points(edata)
		t = ndata[0]
		N = len(ndata)
		dlen = len(Si)
		
		ddm = [0]*len(Si)
		for i in range(1,N):
			dd = ndata[i]
			for j in range(dlen):
				ddm[j] = dd[j] + ddm[j]
				
		ddm = [ ddm[j]/(N-1) for j in range(dlen)]

		plt.figure(figsize=(9.68,3.95))
		plt.xlabel("time (sec)")
		plt.ylabel("conc")		
		lines = []
		for j in range(dlen):
			line=plt.plot(t,ddm[j],label=Si[j])
			lines.append(line)
		plt.legend()
		fig = plt.gcf()
		plt.tight_layout()
		globals.plotted.append([plt.gca(),fig,lines])
		fig_canvas_agg = draw_figure(items,fig)		
		plt.close()