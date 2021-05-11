import matplotlib.pyplot as plt
import multiprocessing
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d

from draw_figure import *

def plot_traj(data,Si,items,plotted,mix_plot=True,logx=False,logy=False,normalize=False,SiTicked=None):
	miter = len(data)
	col = ['C'+str(i) for i in range(100)]
	if SiTicked == None:
		SiTicked = range(len(Si))
	if mix_plot:
		plt.figure(figsize=(9.68,3.95))
		plt.xlabel("time (sec)")
		plt.ylabel("conc") 	
		lines = []			
		if logx:
			plt.xscale('log')	
		if logy:
			plt.yscale('log') 						
		for j in range(miter):
			
			for i in SiTicked:
				if normalize:
					line = plt.plot(data[j][:,0],data[j][:,i+1]/(np.max(data[j][:,i+1])+1.0e-30),col[i])
				else:
					line = plt.plot(data[j][:,0],data[j][:,i+1],col[i])	
		plt.legend([Si[i] for i in SiTicked])
		plt.tight_layout()
		fig = plt.gcf() 
		lines.append(line)
		plotted.append([plt.gca(),fig,lines])
		fig_canvas_agg = draw_figure(items,fig)			
		plt.close()
	else:
		lines = []	
		for i in range(len(Si)):
			plt.figure(figsize=(9.68,3.95))
			plt.xlabel("time")
			plt.ylabel("conc") 
			if logx:
				plt.xscale('log')
			if logy:
				plt.yscale('log')  
			for j in range(miter):  
				if normalize:
					line = plt.plot([ t/3600 for t in data[j][0] ],data[j][1][:,i]/(np.max(data[j][1][:,i])+1.0e-30),col[i])
				else:
					line = plt.plot([ t/3600 for t in data[j][0] ],data[j][1][:,i],col[i])
				lines.append(line)
				
			plt.legend([Si[i] for i in SiTicked])
			plt.tight_layout()
			fig = plt.gcf() 
			plotted.append([plt.gca(),fig,lines])
			fig_canvas_agg = draw_figure(items,fig)	
			plt.close()

def plot_traj2(data,Si,items,plotted,logx=False,logy=False,normalize=False,xlabel="conc",ylabel="conc",zlabel="conc",trange=None):
	miter = len(data)
	col = ['C'+str(i) for i in range(100)]
	plt.figure(figsize=(9.68,3.95))
	plt.xlabel(xlabel)
	plt.ylabel(ylabel) 	
	if trange:
		tran = str(trange).split(":")
		tran = slice(int(tran[0]),int(tran[1]))
	else:
		tran = slice(0,-1)
	if zlabel == "None":
		lines = []			
		if logx:
			plt.xscale('log')	
		if logy:
			plt.yscale('log') 						
		for j in range(miter):
			i1 = -1 if xlabel == "time" else Si.index(xlabel)
			i2 = -1 if ylabel == "time" else Si.index(ylabel)
			if normalize:
				line = plt.plot(data[j][tran,i1+1],data[j][tran,i2+1]/(np.max(data[j][tran,i2+1])+1.0e-30))
			else:
				line = plt.plot(data[j][tran,i1+1],data[j][tran,i2+1])	
		plt.tight_layout()
		fig = plt.gcf() 
		lines.append(line)
		plotted.append([plt.gca(),fig,lines])
		fig_canvas_agg = draw_figure(items,fig)			
		plt.close()
	else:
		ax = plt.axes(projection ='3d') 
		for j in range(miter):
			i1 = -1 if xlabel == "time" else Si.index(xlabel)
			i2 = -1 if ylabel == "time" else Si.index(ylabel)
			i3 = -1 if zlabel == "time" else Si.index(zlabel)
			ax.plot3D(data[j][tran,i1+1], data[j][tran,i2+1], data[j][tran,i3+1]) 
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.set_zlabel(zlabel)
		plt.show()	
	
