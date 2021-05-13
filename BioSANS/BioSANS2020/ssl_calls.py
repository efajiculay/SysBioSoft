import numpy as np
import matplotlib.pyplot as plt

def load_data_traj(file_name):
	with open(file_name,"r") as f:
		data = []
		dd = []
		row1 = str(f.readline()).strip()
		Si = row1.split("\t")[1:]
		for row in f:
			cols = [ float(x) for x in row.split("\t") ]
			if cols[0] == 0.0 and len(dd)>0:
				data.append(np.array(dd))
				dd = []
			dd.append(cols)
		data.append(np.array(dd))
	current_data = (data, Si)
	return current_data

def calc_average_conc_at_tend(edata,points=100):
	data, Si = edata
	N = len(data)
	ddm = 0
	for i in range(1,N):
		dd = data[i][-points:,1:]
		ddm = ddm + dd
		
	ddm = np.mean(ddm/(N-1),0)
	dlen = len(Si)	
	print("\nSpecies Concentration\n")
	for i in range(dlen):
		if str(ddm[i]) not in {"None", "nan", "0.0"}:
			print(Si[i].ljust(50)+" = "+str(ddm[i]))					  

def calc_covariance(edata,points=100):
	data, Si = edata
	N = len(data)
	ddm = 0
	for i in range(1,N):
		dd = data[i][-points:,1:]
		ddm = ddm + dd
		
	ddm = np.mean(ddm/(N-1),0)
	
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
	
	print("covariance\n")
	for i in range(dlen):
		for j in range(i,dlen):	
			if str(vv[i,j]) not in {"None", "nan", "0.0"}:
				label = Si[i]+" and "+Si[j]
				print(label.ljust(50)+" = "+str(vv[i,j]))
				
	vvv = np.zeros((dlen,dlen))
	print("\ncorrelation\n")
	for i in range(dlen):
		for j in range(i,dlen):
			if str(vv[i,j]) not in {"None", "nan", "0.0"}:
				vvv[i,j] = vv[i,j]/np.sqrt(vv[i,i]*vv[j,j])
				label = Si[i]+" and "+Si[j]
				print(label.ljust(50)+" = "+str(vvv[i,j]))				
				
	FF = 0
	CV = 0
	for i in range(N):
		dd = 0
		dd = dd+(data[i][-points:,1:]-ddm)**2
		FF = FF + np.mean(dd,0)/ddm				
		CV = CV + np.mean(np.sqrt(dd),0)/ddm
	FF = FF/(N-1)
	CV = CV/(N-1)
	print("\nFano Factor\n")
	for i in range(dlen):
		if str(FF[i]) not in {"None", "nan", "0.0"}:
			print(Si[i].ljust(50)+" = "+str(FF[i]))				
			
	print("\nCoefficient of variation\n")
	for i in range(dlen):
		if str(CV[i]) not in {"None", "nan", "0.0"}:
			print(Si[i].ljust(50)+" = "+str(CV[i]))						 

def calc_covariance_per_traj(edata,points=100):
	data, Si = edata
	N = len(data)
	ddm = 0
	for i in range(1,N):
		dd = data[i][-points:,1:]
		ddm = ddm + dd
		
	ddm = np.mean(ddm/(N-1),0)
	
	vv = []
	ff = []
	cv = []
	dlen = len(Si)
	for k in range(1,N):
		dd = data[k][-points:,1:] - ddm
		ddd = np.zeros((dlen,dlen))
		fff = np.zeros(dlen)
		cvv = np.zeros(dlen)
		for i in range(dlen):
			for j in range(dlen):
				ddd[i,j] = np.mean(dd[:,i]*dd[:,j])
			fff[i] = ddd[i,i]/ddm[i]					
			cvv[i] = np.sqrt(ddd[i,i])/ddm[i]			
		vv.append(ddd)	
		ff.append(fff)	
		cv.append(cvv)			
	vv = np.array(vv) 
	ff = np.array(ff)   
	cv = np.array(cv)	
	mvv = np.mean(vv,axis=0)
	svv = np.std(vv,axis=0) 
	mff = np.mean(ff,axis=0) 
	sff = np.std(ff,axis=0)  
	mcv = np.mean(cv,axis=0) 
	scv = np.std(cv,axis=0)	  
	print("\ncovariance\n")
	for i in range(dlen):
		for j in range(i,dlen):	
			if str(mvv[i,j]) not in {"None", "nan", "0.0"}:
				label = Si[i]+" and "+Si[j]
				print(label.ljust(50)+" = "+str(mvv[i,j]),"std =",svv[i,j])
				
	plt.figure(figsize=(9.68,3.95))
	plt.xlabel("covariance(P)")
	plt.ylabel("prob") 	 
	plt.hist(vv[:,dlen-1,dlen-1],bins=50, density=True, orientation='vertical')   
	plt.axvline(mvv[dlen-1,dlen-1]) 
	plt.tight_layout()
	plt.savefig(Fname+"_"+Mname+"_covariance(P).jpg")
	plt.close()	
				
	print("\nFano Factor\n")
	for i in range(dlen):
		if str(mff[i]) not in {"None", "nan", "0.0"}:
			print(Si[i].ljust(50)+" = "+str(mff[i]),"std =",sff[i])				
			
	plt.figure(figsize=(9.68,3.95))
	plt.xlabel("Fano factor(P)")
	plt.ylabel("prob") 	 
	plt.hist(ff[:,dlen-1],bins=50, density=True, orientation='vertical')   
	plt.axvline(mff[dlen-1]) 
	plt.tight_layout()
	plt.savefig(Fname+"_"+Mname+"_Fano factor(P).jpg")
	plt.close()  
			
	print("\nCoefficient of variation\n")
	for i in range(dlen):
		if str(mcv[i]) not in {"None", "nan", "0.0"}:
			print(Si[i].ljust(50)+" = "+str(mcv[i]),"std =",scv[i])				
			
	plt.figure(figsize=(9.68,3.95))
	plt.xlabel("Coeff. of Var.(P)")
	plt.ylabel("prob") 	 
	plt.hist(cv[:,dlen-1],bins=50, density=True, orientation='vertical')   
	plt.axvline(mcv[dlen-1]) 
	plt.tight_layout()
	plt.savefig(Fname+"_"+Mname+"_CoeffOfVar(P).jpg")
	plt.close()		

def calc_covariance_bootsrap(edata,points=100,Msamp=1000):
	data, Si = edata
	N = len(data)
	bmvv = []
	bsvv = []
	bmff = []
	bsff = []
	bmcv = []
	bscv = [] 
	mlist = list(range(1,N))
	np.random.shuffle(mlist)
	for ih in range(1000):
		ddm = 0
		curr_list = np.random.choice(mlist,Msamp)
		for i in curr_list:
			dd = data[i][-points:,1:]
			ddm = ddm + dd

		ddm = np.mean(ddm/Msamp,0)
	
		vv = []
		ff = []
		cv = []
		dlen = len(Si)
		for k in curr_list:
			dd = data[k][-points:,1:] - ddm
			ddd = np.zeros((dlen,dlen))
			fff = np.zeros(dlen)
			cvv = np.zeros(dlen)
			for i in range(dlen):
				for j in range(dlen):
					ddd[i,j] = np.mean(dd[:,i]*dd[:,j])
				fff[i] = ddd[i,i]/ddm[i]					
				cvv[i] = np.sqrt(ddd[i,i])/ddm[i]					
			vv.append(ddd)	
			ff.append(fff)	
			cv.append(cvv)	
		vv = np.array(vv) 
		ff = np.array(ff) 
		cv = np.array(cv)		  
		mvv = np.mean(vv,axis=0)
		svv = np.std(vv,axis=0) 
		mff = np.mean(ff,axis=0) 
		sff = np.std(ff,axis=0)
		mcv = np.mean(cv,axis=0) 
		scv = np.std(cv,axis=0)
		bmvv.append(mvv)
		bsvv.append(svv)
		bmff.append(mff)
		bsff.append(sff)
		bmcv.append(mcv)
		bscv.append(scv)		
		
	bmvv = np.array(bmvv)
	bsvv = np.array(bsvv)
	bmff = np.array(bmff)
	bsff = np.array(bsff)
	bmcv = np.array(bmcv)
	bscv = np.array(bscv)	
	
	mvv = np.mean(bmvv,axis=0)	
	svv = np.std(bmvv,axis=0) 
	mff = np.mean(bmff,axis=0)	 
	sff = np.std(bmff,axis=0)
	mcv = np.mean(bmcv,axis=0)	 
	scv = np.std(bmcv,axis=0)	
	
	print("\ncovariance\n")
	for i in range(dlen):
		for j in range(i,dlen):	
			if str(mvv[i,j]) not in {"None", "nan", "0.0"}:
				label = Si[i]+" and "+Si[j]
				print(label.ljust(50)+" = "+str(mvv[i,j]),"std =",svv[i,j])
				
	plt.figure(figsize=(9.68,3.95))
	plt.xlabel("covariance(P)")
	plt.ylabel("prob") 	 
	plt.hist(bmvv[:,dlen-1,dlen-1],bins=50, density=True, orientation='vertical')   
	plt.axvline(mvv[dlen-1,dlen-1]) 
	plt.tight_layout()
	plt.savefig(Fname+"_"+Mname+"_covariance(P)_bootsrap.jpg")
	plt.close()	
				
	print("\nFano Factor\n")
	for i in range(dlen):
		if str(mff[i]) not in {"None", "nan", "0.0"}:
			print(Si[i].ljust(50)+" = "+str(mff[i]),"std =",sff[i])				
			
	plt.figure(figsize=(9.68,3.95))
	plt.xlabel("Fano factor(P)")
	plt.ylabel("prob") 	 
	plt.hist(bmff[:,dlen-1],bins=50, density=True, orientation='vertical')   
	plt.axvline(mff[dlen-1]) 
	plt.tight_layout()
	plt.savefig(Fname+"_"+Mname+"_Fano factor(P)_bootsrap.jpg")
	plt.close()  
			
	print("\nCoefficient of Variation\n")
	for i in range(dlen):
		if str(mcv[i]) not in {"None", "nan", "0.0"}:
			print(Si[i].ljust(50)+" = "+str(mcv[i]),"std =",scv[i])				
			
	plt.figure(figsize=(9.68,3.95))
	plt.xlabel("Coeff. of Var.")
	plt.ylabel("prob") 	 
	plt.hist(bmcv[:,dlen-1],bins=50, density=True, orientation='vertical')   
	plt.axvline(mcv[dlen-1]) 
	plt.tight_layout()
	plt.savefig(Fname+"_"+Mname+"_CoeffofVar(P)_bootsrap.jpg")
	plt.close() 
	
def Prob_density_wtime(edata,Fname,Mname):
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
		ddd = dd[:,j].flatten()
		csq = int(np.sqrt(len(ddd)))
		Zs, Xs, Ys = np.histogram2d(tt, ddd, bins=(csq, csq))
		plt.figure(figsize=(9.68,3.95))
		plt.xlabel("time")
		plt.ylabel("conc("+Si[j]+")")
		delCon = int(np.max(Zs) - np.min(Zs))/500
		if Si[j] == "P":
			plt.yscale('log')
		delCon = int(np.max(Zs) - np.min(Zs))/500
		try:
			cntr = plt.contour(Xs[1:],Ys[1:],Zs.T,list(range(int(np.min(Zs)),int(np.max(Zs))+1,int(delCon))),cmap="Spectral_r")
		except:
			cntr = plt.contour(Xs[1:],Ys[1:],Zs.T,100,cmap="Spectral_r")	   
		fig = plt.gcf()
		fig.colorbar(cntr)
		plt.tight_layout()
		plt.savefig(Fname+"_"+Mname+"_"+Si[j]+"_prob_density_wtime.jpg")
		plt.close()
		
def Prob_density_tslice(edata,bins=50):
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
			plt.savefig(Fname+"_"+Si[j]+"_prob_density_tslice.jpg")
			plt.close()			
			
def Prob_density(edata,Fname):
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
		ddd = dd[:,j]
		line = plt.hist(ddd,bins=int(np.sqrt(len(ddd))), density=True, orientation='vertical')
		plt.legend([Si[j]])
		plt.xscale('log')
		plt.tight_layout()
		fig = plt.gcf()		
		plt.savefig(Fname+"_"+Si[j]+"_prob_density.jpg")
		plt.close()
		