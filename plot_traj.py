import matplotlib.pyplot as plt
import multiprocessing

from draw_figure import *

def plot_traj(data,Si,items,plotted,mix_plot=True,logx=False,logy=False,normalize=False):
	miter = len(data)
	col = ['C'+str(i) for i in range(100)]
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
			
			for i in range(len(Si)):
				if normalize:
					line = plt.plot(data[j][:,0],data[j][:,i+1]/(np.max(data[j][:,i+1])+1.0e-30),col[i])
				else:
					line = plt.plot(data[j][:,0],data[j][:,i+1],col[i])	
		plt.legend(Si)
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
			plt.legend([Si[i]])
			plt.tight_layout()
			fig = plt.gcf() 
			plotted.append([plt.gca(),fig,lines])
			fig_canvas_agg = draw_figure(items,fig)	
			plt.close()