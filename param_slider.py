from scipy.integrate import odeint
from propensity import *
from recalculate_globals import *
import mglobals as globals2
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
from tkinter import filedialog

def load_data():
	try:
		file = filedialog.askopenfilename(title = "Select file")
		with open(file,"r") as f:
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
			return (data, Si)
	except:
		return None

def ParamODE_model(z,t,Sp,Ks,Rr,Rp,V,molar=False):  
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
	
def update_range(dk,valc):
	dk.valmin = valc[0]
	dk.valmax = valc[1]
	dk.ax.set_xlim(dk.valmin,dk.valmax)
	
def submit(text,dk):
	valc = [eval(x) for x in text.split(",")]
	update_range(dk,valc)
	
def update(val,Ks,KKglobals,fig,Si,l,Sp,Rr,Rp,V,molar,t,z):
	cc = 0
	for ih in range(len(Ks)):
		if len(Ks[ih]) == 1:
			Ks[ih][0] = KKglobals[cc].val
			cc = cc + 1
		else:
			Ks[ih][0] = KKglobals[cc].val
			cc = cc + 1
			Ks[ih][1] = KKglobals[cc].val
			cc = cc + 1
	res = odeint(ParamODE_model,z,t, args=(Sp,Ks,Rr,Rp,V,molar))
	for ih in range(len(z)):
		l[ih].set_ydata(res[:,ih])
		l[ih].set_label(Si[ih])
	fig.canvas.draw_idle()
	
def ParamODE_int(conc,t,Sp,Ks,Rr,Rp,V,molar=False,rfile="",setP=[]): 
	global CKglobals, KKglobals

	Spc = [s for s in Sp]
	z = [conc[a] for a in Sp]
	Cmiss = []
	for i in range(len(z)):
		if z[i]<0:
			Cmiss.append(Spc[i])	
	Kmiss = []
	for i in range(len(Ks)):
		if len(Ks[i]) == 1:
			if Ks[i][0]<0:
				Kmiss.append((i,0))
				Ks[i][0] = 2
		elif len(Ks[i]) == 2:
			if Ks[i][0]<0:
				Kmiss.append((i,0))
				Ks[i][0] = 2
			if Ks[i][1]<0:
				Kmiss.append((i,1))	
				Ks[i][0] = 2
	npar = len(Cmiss) + len(Kmiss)

	data2 = load_data()	
	if data2:	
		Si2 = data2[1]
		t = data2[0][0][:,0]
		data = data2[0][0][:,1:]

	get_globals(rfile)
	z = [conc[a] for a in Sp]
	Si=[Si for Si in Sp]
	
	nz = []
	reserve_events_words = {"t","time","status","status2","timer","finish","delay","dtime"}
	for row in range(V.shape[0]):
		key = Si[row].strip().split("_")[0]
		if key not in reserve_events_words:
			nz.append(row)	   
	Sp = {Si[z]:Sp[Si[z]] for z in nz}
	Si = [Si for Si in Sp]		
	
	scount = 0
	for ih in range(len(Ks)):
		if len(Ks[ih]) == 1:
			scount = scount + 1
		else:
			scount = scount + 2
	
	res = odeint(ParamODE_model,z,t, args=(Sp,Ks,Rr,Rp,V,molar)) 
	fig, ax = plt.subplots(figsize=(12,4))
	plt.subplots_adjust(left=0.1, right=0.70)
	if data2:
		plt.ylim(0,np.max(data)*1.05)
	if setP[0] == True:
		plt.xscale("log")
	if setP[1] == True:
		plt.yscale("log")
	
	if data2:	
		for ih in range(len(Si2)):
			plt.plot(t, data[:,ih], lw=1,ls='--',label=Si2[ih]+"_True")[0]
	
	l = []
	for ih in nz:
		l.append(plt.plot(t, res[:,ih], lw=2,label=Si[ih])[0])
	plt.legend()
		
	axcolor = 'lightgoldenrodyellow'
	KKglobals = []
	CKglobals = []
	cc = 0
	for ih in range(len(Ks)):
		strt1 = round(abs(Ks[ih][0]*0.10),2)
		ends1 = round(abs(Ks[ih][0]*1.90),2)
		if len(Ks[ih]) == 1:
			sk = plt.axes([0.73, 0.93-0.05*cc, 0.1, 0.04], facecolor=axcolor)
			sd = Slider(sk, 'kf'+str(ih+1), strt1, ends1, valinit=Ks[ih][0], valstep=Ks[ih][0]/100)
			KKglobals.append(sd)
			CKglobals.append(TextBox(plt.axes([0.93, 0.93-0.05*cc, 0.08, 0.04]),'range',initial=str(strt1)+', '+str(ends1)))
			cc = cc + 1
		else:
			sk = plt.axes([0.73, 0.93-0.05*cc, 0.1, 0.04], facecolor=axcolor)
			sd = Slider(sk, 'kf'+str(ih+1), strt1, ends1, valinit=Ks[ih][0], valstep=Ks[ih][0]/100)
			KKglobals.append(sd)	
			CKglobals.append(TextBox(plt.axes([0.93, 0.93-0.05*cc, 0.08, 0.04]),'range',initial=str(strt1)+', '+str(ends1)))
			cc = cc + 1
			
			strt2 = round(abs(Ks[ih][1]*0.10),2)
			ends2 = round(abs(Ks[ih][1]*1.90),2)
			
			sk = plt.axes([0.73, 0.93-0.05*cc, 0.1, 0.04], facecolor=axcolor)
			sd = Slider(sk, 'kb'+str(ih+1), strt2, ends2, valinit=Ks[ih][1], valstep=Ks[ih][1]/100)
			KKglobals.append(sd)	
			CKglobals.append(TextBox(plt.axes([0.93, 0.93-0.05*cc, 0.08, 0.04]),'range',initial=str(strt2)+', '+str(ends2)))
			cc = cc + 1
	
	for ih in range(len(CKglobals)):
		KKglobals[ih].on_changed(lambda val : update(val,Ks,KKglobals,fig,Si,l,Sp,Rr,Rp,V,molar,t,z))
	
	FunnyEvents = []
	for ih in range(len(CKglobals)):	
		funny = "CKglobals["+str(ih)+"].on_submit(lambda v : submit(v,KKglobals["+str(ih)+"]))"  #Why do this - because ih lambda get is only the last value
		FunnyEvents.append(funny)
		
	for fun in FunnyEvents:
		exec(fun,globals())
	
	plt.show()
	return [0]