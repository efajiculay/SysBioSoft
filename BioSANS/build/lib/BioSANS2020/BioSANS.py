import os
from datetime import datetime
import tkinter as gui
from tkinter import ttk
from tkinter import filedialog
from pathlib import Path
from PIL import Image as Image2, ImageTk
import numpy as np
import time
import threading
from queue import Queue
from math import ceil as Myceil
from sys import platform

import mglobals as globals
import proc_global as proc_global
from prepare_canvas import *
from process import *
from plot_traj import *
from Transform import *
from process_sbml import process_sbml as sbml_to_topo

from topology_view import *
from new_file import *

globals.init()
if __name__ == '__main__':
	proc_global.init()

top = gui.Tk()
top.title("BioSANS 1.0")
top.geometry("1005x550")
top.resizable(False, False)

header = gui.Label(top, text="BioSANS")
header.configure(
	bg = "green",
	fg = "white",
	height = 1,
	width = 1005,
	font = "Helvetica 18 bold italic"
)
header.pack()

frame = gui.Frame(top)
frame.configure(
	bg = "light cyan", 
	borderwidth=2,
	height=500,
	width=1005
)
frame.pack()

footer = gui.Label(top, text="Biological Stochastic Simulation Algorithms")
footer.configure(
	bg = "green",
	fg = "white",
	width = 1005,
	font = "Helvetica 10 bold italic",
	anchor='w'
)
footer.pack()

file_name = {}
def load_data(items):
	global file_name
	file = filedialog.askopenfilename(title = "Select file")
	file_name["topology"] = file
	globals.toConvert = file
	if os.path.isfile(file):
		file_name['last_open'] = view_topo(file,items)
		
def create_file(items):
	global file_name
	try:
		os.mkdir("Temporary_folder",777)
	except:
		for item in Path("Temporary_folder").iterdir():
			if item.is_dir():
				pass
			else:
				item.unlink()
		
	file_name['last_open'] = new_file(items)
	file_name["topology"] = "Temporary_folder/temp.txt"
	
def save_file():
	global file_name
	file = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
	file_name["topology"] = file.name
	if file is None:
		return
	file.write(file_name['last_open'].get("0.0",END))
	file.close()
	
def load_data2(plot=False):
	t_o = time.time()
	global file_name, current_data
	file = filedialog.askopenfilename(title = "Select file")
	file_name["trajectory"] = file
	with open(file_name["trajectory"],"r") as f:
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
	if plot:
		plot_traj(data,Si,items,globals.plotted,mix_plot=True,logx=False,logy=False,normalize=False)
	current_data = (data, Si)
	print(time.time()-t_o)
	
def tload_data2(plot=False):
	if __name__ == '__main__':
		t = threading.Thread(target=load_data2, args=(plot,), daemon=False)
		t.start()
	
def load_image(wdata=False):
	t_o = time.time()
	global items, current_data
	canvas,scroll_x,scroll_y = items
	file = filedialog.askopenfilename(title = "Select file")
	load = Image2.open(file)
	render = ImageTk.PhotoImage(load)
	img = gui.Label(canvas, image=render)
	img.image = render
	canvas.create_window(0, 426*globals.plot_i, anchor='nw', window=img)
	canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
	canvas.configure(scrollregion=canvas.bbox("all"))		
	globals.plot_i = globals.plot_i+1
	
	if wdata:
		file = str(file).replace("jpg", "dat").replace("png", "dat")
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
		current_data = (data, Si)
	print(time.time()-t_o)
			
def eval2(x):
	try:
		return eval(x)
	except:
		try:
			par = str(x).lower().capitalize()
			return eval(par)
		except:
			return str(x)

current_data = None
def dict_trans(x1):
	x2 = x1.split(",")
	x3 = {}
	try:
		for x in x2:
			r = x.split("=")
			x3[r[0].strip()] = float(r[1])
	except:
		pass
	return x3
	
def convert(x, con):
	try:
		return con(x)
	except:
		return x
	
def range_trans(x1):
	x3 = []
	x2 = x1.split(",")
	x2 = [ convert(x2[x],float) for x in range(len(x2)) ]
	if len(x2)>=4:
		try:
			w = 'linspace'
			if len(x2) == 5:
				w = x2[4]
			if w == 'linspace':
				x4 = np.linspace(x2[1],x2[2],x2[3])
			if w == 'logspace':
				x4 = np.logspace(x2[1],x2[2],x2[2])
			if w == 'uniform':
				x4 = np.random.uniform(x2[1],x2[2],int(x2[3]))		
			if w == 'normal':
				x4 = np.random.normal(x2[1],x2[2],int(x2[3]))			
			if w == 'lognormal':
				x4 = np.random.lognormal(x2[1],x2[2],int(x2[3]))				
			for x in x4:
				x3.append(x)
			x3.append(x2[0])
		except:
			cc = ",".join([str(x) for x in x2[1:]])
			x3 = list(eval(cc))
			x3.append(x2[0])
	return x3
	
def range_prep(x1):
	x3 = []
	x2 = x1.split(",")
	x2 = [ convert(x2[x],float) for x in range(len(x2)) ]
	if len(x2)>=4:
		cc = x2[0].lower().split("f")
		r2 = 0
		if len(cc)<2:
			cc = x2[0].lower().split("b")
			r2 = 1
		r1 = int(cc[1])-1
		try:
			w = 'linspace'
			if len(x2) == 5:
				w = x2[4]
			if w == 'linspace':
				x4 = np.linspace(x2[1],x2[2],x2[3])
			if w == 'logspace':
				x4 = np.logspace(x2[1],x2[2],x2[2])
			if w == 'uniform':
				x4 = np.random.uniform(x2[1],x2[2],int(x2[3]))		
			if w == 'normal':
				x4 = np.random.normal(x2[1],x2[2],int(x2[3]))			
			if w == 'lognormal':
				x4 = np.random.lognormal(x2[1],x2[2],int(x2[3]))				
			for x in x4:
				x3.append(x)
		except:
			cc = ",".join([str(x) for x in x2[1:]])
			x3 = list(eval(cc))
		return [[r1,r2],x3]
	return x1

def mrun(par,E,defs2):
	global items,current_data
	defs = []
	for i in range(len(E)):
		try:
			val = eval2(E[i].get())
			defs.append(val)
		except:
			val = eval2(defs2[i].get())
			defs.append(val)
	defs[15] = dict_trans(E[15].get())
	defs[16] = range_trans(E[16].get())
	defs[17] = range_prep(E[17].get())
	if not defs[9] in ["Analyt","SAnalyt","Analyt-ftx","SAnalyt-ftk","k_est1","k_est2","k_est3","k_est4","k_est5","k_est6","k_est7","k_est8","k_est9","k_est10","k_est11","NetLoc1","NetLoc2"]:
		with open(defs[13]+"_"+defs[9]+"_params.dat","w") as f:
			f.write("\n".join([str(x) for x in defs]))
	if defs[9] in ["k_est1","k_est2","k_est3","k_est4","k_est6","k_est7","k_est8","k_est9","k_est10","k_est11","LNA2","LNA3","Analyt","Analyt-ftx","SAnalyt","SAnalyt-ftk","Analyt2","topoTosbml","topoTosbml2","topoTosbml3","LNA-vs","LNA-ks","LNA-xo","NetLoc1","NetLoc2"]:
		par.destroy()
	defs.append(items)
	current_data = tprocess(defs)
	
def tprocess(defs):
	if defs[9] == "k_est5":
		process(*defs)
	else:
		if __name__ == '__main__':
			t = threading.Thread(target=lambda: process(*defs), daemon=False)
			t.start()
			
def plot_trajD(current_data,items):
	try:
		data, Si = current_data
		par = gui.Toplevel()
		par.resizable(False, False)
		par.wm_title("Plot settings")
		
		L = [ 
			gui.Label(par, text="choose x-axis",fg="blue"),
			gui.Label(par, text="choose y-axis",fg="blue"),
			gui.Label(par, text="choose z-axis",fg="blue"),
			gui.Label(par, text="choose step-range",fg="blue")
		]
		[ L[i].grid(row = i, column = 0, sticky = gui.W, pady = 2) for i in range(4) ]
		optsel = [ [par,gui.StringVar(),"None","time"] + Si for i in range(3) ]
		optselVar = [ optsel[i][1] for i in range(3) ]
		[ optselVar[i].set("None") for i in range(3) ]
		E = [ gui.OptionMenu(*optsel[i]) for i in range(3) ]
		[ E[i].config(width = 14) for i in range(3) ]
		E.append(gui.Entry(par, bd =5))
		[ E[i].grid(row = i, column = 1, sticky = gui.W, pady = 2) for i in range(4) ]
		E[-1].insert(gui.END,"0:-1")
		
		B1 = ttk.Button(par, text="PLOT", command=lambda : plot_traj2( \
		data,Si,items,globals.plotted,logx=False,logy=False,normalize=False, \
		xlabel=optselVar[0].get() ,ylabel=optselVar[1].get(), zlabel=optselVar[2].get(), trange=E[-1].get()))
		B1.grid(row = 5, column = 0, sticky = gui.W, pady = 2) 
	except:
		gui.messagebox.showinfo("Trajectory not loaded yet", "Please load the trajectory. BioSANS save it into a file during your last run.")
	
def getChecked(L1,L,Si):
	checkSi = []
	for i in range(len(L1)):
		key = L1[i].get()
		if key != "0":
			checkSi.append(Si.index(key))
	return checkSi
	
def plot_trajD2(current_data,items):
	try:
		data, Si = current_data
		par = gui.Toplevel()
		par.resizable(False, False)
		par.wm_title("Plot species")
		
		Ls = len(Si)
		L1 = [gui.StringVar() for i in range(Ls)]
		[ L1[i].set("0")  for i in range(Ls)]
		L = [ gui.Checkbutton(par, text=Si[i],variable=L1[i], onvalue=Si[i], offvalue="0") for i in range(Ls) ]
		[ L[i].grid(row = i, column = 0, sticky = gui.W, pady = 2) for i in range(Myceil(Ls/2)) ]
		[ L[i].grid(row = i-Myceil(Ls/2), column = 1, sticky = gui.W, pady = 2) for i in range(Myceil(Ls/2),Ls) ]
		
		B1 = ttk.Button(par, text="PLOT", \
		command=lambda : \
				plot_traj(data,Si,items,globals.plotted,mix_plot=True,logx=False,logy=False,normalize=False,SiTicked=getChecked(L1,L,Si)) \
		)
		B1.grid(row = Myceil(Ls/2), column = 0, sticky = gui.W, pady = 2)
	except:
		gui.messagebox.showinfo("Trajectory not loaded yet", "Please load the trajectory. BioSANS save it into a file during your last run.")

def paramSet(method):
	global file_name, items

	if file_name["topology"] == "Temporary_folder/temp.txt":
		with open(file_name["topology"],"w") as ff:
			ff.write(file_name['last_open'].get("0.0",END))
	
	path = Path(file_name["topology"])
	ss = str(file_name["topology"]).split("/")
	ss = ss[-1] if len(ss)>1 else ""
	name = str(path.parent)+"/"+ss+"_"+datetime.now().strftime("%Y%m%d_%H%M%S")
	par = gui.Toplevel()
	par.resizable(False, False)
	par.wm_title("Parameter setting")
	opts = [
		"File name",
		"Number of iteration :",
		"File Units? :",
		"Volume (L) :",
		"end time (tn) :",
		"tau-scaler",
		"Normalized",
		"logx",
		"logy",
		"method",
		"tsteps",
		"mix_plot",
		"save",
		"out fname",
		"show plot",
		"modify Cini",
		"Cini range",
		"K-range",
		"mult proc",
		"Implicit"
	]
	defs = [file_name["topology"],1,gui.StringVar(),1.0,100,1.5,False,False,False,method,1000,True,True,name,True,"","","",False,False]
	defs[2].set('molecules')
	
	topfile = open(file_name["topology"],"r")
	for row in topfile:
		if row[0] == "#":
			gg = row.split(",")[1:]
			for x in gg:
				xx = [ g.strip() for g in x.split("=")  ] 
				if xx[0] == "Volume":
					defs[3] = xx[1]
				elif xx[0] == "tend":
					defs[4] = xx[1]
				elif xx[0] == "FileUnit":
					defs[2].set(xx[1])
				elif xx[0] == "logx":
					defs[7] = xx[1]
				elif xx[0] == "Normalized":
					defs[6] = xx[1]
				elif xx[0] == "steps":
					defs[10] = xx[1]					

	oplen = len(opts)
	L = [ gui.Label(par, text=opts[i],fg="blue") for i in range(oplen) ]
	[ L[i].grid(row = i, column = 0, sticky = gui.W, pady = 2) for i in range(10) ]
	[ L[10+i].grid(row = i, column = 2, sticky = gui.W, pady = 2) for i in range(10) ]
	E = [ gui.Entry(par, bd =5) if i!=2 else gui.OptionMenu(par, defs[2],'molecules', 'molar', 'moles') for i in range(oplen) ]
	[ E[i].grid(row = i, column = 1, sticky = gui.W, pady = 2) for i in range(10) ]
	[ E[10+i].grid(row = i, column = 3, sticky = gui.W, pady = 2) for i in range(10) ]
	[ E[i].insert(gui.END,str(defs[i])) if i!=2 else None for i in range(oplen) ]
	E[9].configure(state="disable")
	E[2].config(width = 14)
	if method == "ODE": 
		E[5].configure(state="disable")
	if method == "ODE2": 
		E[5].configure(state="disable")
	elif method == "Gillespie_": 
		E[5].configure(state="disable")	
	elif method == "CLE": 
		E[5].delete(0, gui.END)
		E[5].insert(gui.END,str(10))	
	
	B1 = ttk.Button(par, text="RUN", command=lambda : mrun(par,E,defs))
	B1.grid(row = oplen, column = 0, sticky = gui.W, pady = 2) 
	if method in ["k_est1","k_est2","k_est3","k_est4","k_est6","k_est7","k_est8","k_est9","k_est10","k_est11","LNA2","LNA3","LNA-vs","LNA-ks","LNA-xo","NetLoc1","NetLoc2","Analyt","SAnalyt-ftk","SAnalyt","Analyt-ftx","Analyt2","topoTosbml","topoTosbml2","topoTosbml3"]:
		B1.invoke()
	
if __name__ == "__main__":

	menubut1 = gui.Menubutton(frame,text=" File/Model ",activebackground="#f2f20d",activeforeground="red",bg="#00cc00",fg="white" if platform.lower() != "darwin" else "green")
	menubut1.menu =  gui.Menu ( menubut1, tearoff = 1 )
	menubut1["menu"] =  menubut1.menu
	LoadMenu = gui.Menu(frame,tearoff = 1 )
	LoadMenu.add_command ( label="Topology/File",command=lambda: load_data(items),background="white",foreground="Blue"  )
	LoadMenu.add_command ( label="Trajectory file",command=tload_data2,background="white",foreground="Blue" )
	LoadMenu.add_command ( label="Traj. w/ plot",command=lambda: tload_data2(True),background="white",foreground="Blue" )
	LoadMenu.add_command ( label="Image of plot",command=load_image,background="white",foreground="Blue" )
	LoadMenu.add_command ( label="Image w/ data",command=lambda: load_image(True),background="white",foreground="Blue" )
	menubut1.menu.add_cascade(label="Open", menu=LoadMenu)	
	menubut1.menu.add_command ( label="New File",command=lambda: create_file(items))
	menubut1.menu.add_command ( label="Save File",command=lambda: save_file() )
	menubut1.place(x=2,y=5)

	menubut2 = gui.Menubutton(frame,text="Propagation",activebackground="#f2f20d",activeforeground="red",bg="#00cc00",fg="white" if platform.lower() != "darwin" else "green")
	menubut2.menu =  gui.Menu ( menubut2, tearoff = 1 )
	menubut2["menu"] =  menubut2.menu
	AnalMenu = gui.Menu(frame,tearoff = 1 )
	AnalMenu.add_command ( label="Pure Symbolic :f(t,xo,k)",command=lambda: paramSet("Analyt"),background="white",foreground="Blue"  )
	AnalMenu.add_command ( label="Semi-Symbolic :f(t)",command=lambda: paramSet("SAnalyt"),background="white",foreground="Blue"  )
	AnalMenu.add_command ( label="Semi-Symbolic :f(t,xo)",command=lambda: paramSet("Analyt-ftx"),background="white",foreground="Blue"  )
	AnalMenu.add_command ( label="Semi-Symbolic :f(t,k)",command=lambda: paramSet("SAnalyt-ftk"),background="white",foreground="Blue"  )
	AnalMenu.add_command ( label="For wxmaxima",command=lambda: paramSet("Analyt2"),background="white",foreground="Blue"  )
	menubut2.menu.add_cascade(label="Analytical soln.", menu=AnalMenu)
	
	ODEMenu = gui.Menu(frame,tearoff = 1 )
	ODEMenu.add_command ( label="Molecules(micro)",command=lambda: paramSet("ODE-1"),background="white",foreground="Blue"  )
	ODEMenu.add_command ( label="Molar(macro)",command=lambda: paramSet("ODE-2"),background="white",foreground="Blue"  )
	ODEMenu.add_command ( label="Mole(macro)",command=lambda: paramSet("ODE-3"),background="white",foreground="Blue"  )
	menubut2.menu.add_cascade(label="ODE int", menu=ODEMenu)
	
	Rungek4 = gui.Menu(frame,tearoff = 1 )
	Rungek4.add_command ( label="Molecules(micro)",command=lambda: paramSet("rk4-1"),background="white",foreground="Blue"  )
	Rungek4.add_command ( label="Molar(macro)",command=lambda: paramSet("rk4-2"),background="white",foreground="Blue"  )
	Rungek4.add_command ( label="Mole(macro)",command=lambda: paramSet("rk4-3"),background="white",foreground="Blue"  )
	menubut2.menu.add_cascade(label="RK4-fix-interval", menu=Rungek4)
	
	Rungek4a = gui.Menu(frame,tearoff = 1 )
	Rungek4a.add_command ( label="Molecules(micro)",command=lambda: paramSet("rk4-1a"),background="white",foreground="Blue"  )
	Rungek4a.add_command ( label="Molar(macro)",command=lambda: paramSet("rk4-2a"),background="white",foreground="Blue"  )
	Rungek4a.add_command ( label="Mole(macro)",command=lambda: paramSet("rk4-3a"),background="white",foreground="Blue"  )
	menubut2.menu.add_cascade(label="RK4-tau-adaptive", menu=Rungek4a)	
	
	EulrTau = gui.Menu(frame,tearoff = 1 )
	EulrTau.add_command ( label="Molecules(micro)",command=lambda: paramSet("Euler-1"),background="white",foreground="Blue" )
	EulrTau.add_command ( label="Molar(macro)",command=lambda: paramSet("Euler-2"),background="white",foreground="Blue" )
	EulrTau.add_command ( label="Mole(macro)",command=lambda: paramSet("Euler-3"),background="white",foreground="Blue" )
	menubut2.menu.add_cascade(label="Euler (tau-adaptive-1)", menu=EulrTau)
	
	EulrTau2 = gui.Menu(frame,tearoff = 1 )
	EulrTau2.add_command ( label="Molecules(micro)",command=lambda: paramSet("Euler2-1"),background="white",foreground="Blue" )
	EulrTau2.add_command ( label="Molar(macro)",command=lambda: paramSet("Euler2-2"),background="white",foreground="Blue" )
	EulrTau2.add_command ( label="Mole(macro)",command=lambda: paramSet("Euler2-3"),background="white",foreground="Blue" )
	menubut2.menu.add_cascade(label="Euler (tau-adaptive-2)", menu=EulrTau2)	
	
	#Itoints = gui.Menu(frame,tearoff = 1 )
	#Itoints.add_command ( label="Molecules",command=lambda: paramSet("Itoint-1"),background="white",foreground="Blue" )
	#Itoints.add_command ( label="Molar",command=lambda: paramSet("Itoint-2"),background="white",foreground="Blue" )
	#Itoints.add_command ( label="Mole",command=lambda: paramSet("Itoint-3"),background="white",foreground="Blue" )
	#menubut2.menu.add_cascade(label="Itoint", menu=Itoints)
	
	#Stratint = gui.Menu(frame,tearoff = 1 )
	#Stratint.add_command ( label="Molecules",command=lambda: paramSet("Stratint-1"),background="white",foreground="Blue" )
	#Stratint.add_command ( label="Molar",command=lambda: paramSet("Stratint-2"),background="white",foreground="Blue" )
	#Stratint.add_command ( label="Mole",command=lambda: paramSet("Stratint-3"),background="white",foreground="Blue" )
	#menubut2.menu.add_cascade(label="Stratint", menu=Stratint)

	cletauA = gui.Menu(frame,tearoff = 1 )
	cletauA.add_command ( label="Molecules(micro)",command=lambda: paramSet("CLE"),background="white",foreground="Blue"  )	
	menubut2.menu.add_cascade(label="CLE (tau-adaptive)", menu=cletauA)
	
	cletauA2 = gui.Menu(frame,tearoff = 1 )
	cletauA2.add_command ( label="Molecules(micro)",command=lambda: paramSet("CLE2"),background="white",foreground="Blue"  )	
	menubut2.menu.add_cascade(label="CLE (cle-fixIntvl)", menu=cletauA2)
	
	TauLMenu = gui.Menu(frame,tearoff = 1 )
	TauLMenu.add_command ( label="Tau-leapingV1-micro",command=lambda: paramSet("Tau-leaping"),background="white",foreground="Blue" )
	TauLMenu.add_command ( label="Tau-leapingV2-micro",command=lambda: paramSet("Tau-leaping2"),background="white",foreground="Blue" )
	TauLMenu.add_command ( label="Sim-TauLeap-micro",command=lambda: paramSet("Sim-TauLeap"),background="white",foreground="Blue" )
	menubut2.menu.add_cascade(label="Tau-leaping-micro", menu=TauLMenu)
	
	LNAMenu = gui.Menu(frame,tearoff = 1 )
	LNAMenu.add_command ( label="COV-time-dependent, Macroscopic",command=lambda: paramSet("LNA(t)"),background="white",foreground="Blue" )
	LNAMenu.add_command ( label="FF-time-dependent, Macroscopic",command=lambda: paramSet("LNA2(t)"),background="white",foreground="Blue" )
	LNAMenu.add_command ( label="Numeric, values",command=lambda: paramSet("LNA"),background="white",foreground="Blue" )
	LNAMenu.add_command ( label="Symbolic, Microscopic",command=lambda: paramSet("LNA2"),background="white",foreground="Blue" )
	LNAMenu.add_command ( label="Symbolic, Macroscopic",command=lambda: paramSet("LNA3"),background="white",foreground="Blue" )
	LNAMenu.add_command ( label="Symbolic, f(xo), Macroscopic",command=lambda: paramSet("LNA-xo"),background="white",foreground="Blue" )
	LNAMenu.add_command ( label="Symbolic, f(ks), Macroscopic",command=lambda: paramSet("LNA-ks"),background="white",foreground="Blue" )
	LNAMenu.add_command ( label="Symbolic, values, Macroscopic",command=lambda: paramSet("LNA-vs"),background="white",foreground="Blue" )
	menubut2.menu.add_cascade(label="Linear Noise Appx.", menu=LNAMenu)
	
	GilMenu = gui.Menu(frame,tearoff = 1 )
	GilMenu.add_command ( label="Direct method",command=lambda: paramSet("Gillespie_"),background="white",foreground="Blue"  )
	#GilMenu.add_command ( label="First Rxn Method",command=lambda: print("Not implemented yet"),background="white",foreground="Blue"  )
	#GilMenu.add_command ( label="Next Rxn Method",command=lambda: print("Not implemented yet"),background="white",foreground="Blue"  )
	#GilMenu.add_command ( label="Optimized Direct",command=lambda: print("Not implemented yet"),background="white",foreground="Blue"  )
	menubut2.menu.add_cascade(label="Gillespie", menu=GilMenu)
	
	ParEsMenu = gui.Menu(frame,tearoff = 1 )
	ParEsMenu.add_command ( label="Nelder-Mead (NM), Macroscopic",command=lambda: paramSet("k_est6"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="Nelder-Mead (NM), Microscopic",command=lambda: paramSet("k_est7"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="Powell, Macroscopic",command=lambda: paramSet("k_est8"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="Powell, Microscopic",command=lambda: paramSet("k_est9"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="L-BFGS-B, Macroscopic",command=lambda: paramSet("k_est10"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="L-BFGS-B, Microscopic",command=lambda: paramSet("k_est11"),background="white",foreground="Blue" )		
	ParEsMenu.add_command ( label="NM-Diff. Evol., Macroscopic",command=lambda: paramSet("k_est3"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="NM-Diff. Evol., Microscopic",command=lambda: paramSet("k_est4"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="Parameter slider/scanner",command=lambda: paramSet("k_est5"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="MCEM, Macroscopic",command=lambda: paramSet("k_est1"),background="white",foreground="Blue" )
	ParEsMenu.add_command ( label="MCEM, Microscopic",command=lambda: paramSet("k_est2"),background="white",foreground="Blue" )
	menubut2.menu.add_cascade(label="Estimate Params", menu=ParEsMenu)
	menubut2.place(x=95,y=5)
	
	menubut3 = gui.Menubutton(frame,text="    Analysis    ",activebackground="#f2f20d",activeforeground="red",bg="#00cc00",fg="white" if platform.lower() != "darwin" else "green")
	menubut3.menu =  gui.Menu ( menubut3, tearoff = 1 )
	menubut3["menu"] =  menubut3.menu
	menubut3.menu.add_command ( label="Covariance",command=lambda: calc_covariance(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Fano Factor",command=lambda: fano_factor(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Cross correlation",command=lambda: calc_cross_corr(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Probability density",command=lambda: Prob_density(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Freq. Dist w/r to t",command=lambda: Prob_density2(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Hist. slice of time",command=lambda: Prob_density3(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Average of traj.",command=lambda: Ave_traj(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Phase portrait",command=lambda: plot_trajD(current_data,items),background="white",foreground="Blue"  )
	menubut3.menu.add_command ( label="Plot Data",command=lambda: plot_trajD2(current_data,items),background="white",foreground="Blue"  )
	ConvMenu = gui.Menu(frame,tearoff = 1 )
	ConvMenu.add_command ( label="Topo(molecules) to SBML",command=lambda: paramSet("topoTosbml"),background="white",foreground="Blue"  )	
	ConvMenu.add_command ( label="Topo(molar) to SBML",command=lambda: paramSet("topoTosbml2"),background="white",foreground="Blue"  )	
	ConvMenu.add_command ( label="Topo(no unit) to SBML",command=lambda: paramSet("topoTosbml3"),background="white",foreground="Blue"  )	
	ConvMenu.add_command ( label="SBML to Topo",command=lambda: sbml_to_topo(globals.toConvert),background="white",foreground="Blue"  )	
	NetLMenu = gui.Menu(frame,tearoff = 1 )
	NetLMenu.add_command ( label="Symbolic, Macroscopic",command=lambda: paramSet("NetLoc1"),background="white",foreground="Blue" )
	NetLMenu.add_command ( label="Numeric, Macroscopic",command=lambda: paramSet("NetLoc2"),background="white",foreground="Blue" )
	menubut3.menu.add_cascade(label="Network Localization", menu=NetLMenu)
	menubut3.menu.add_cascade(label="Convert", menu=ConvMenu)	
	menubut3.place(x=189,y=5)	

	frame1 = gui.Frame(frame, height = 435, width = 972, bg='#8c8c8c', borderwidth=2)
	frame1.place(x=2,y=35) 
	items = prepare_frame_for_plot(frame1,972,435)

	top.mainloop()