import sys, os
import glob
import matplotlib.pyplot as plt
from func_timeout import func_timeout, FunctionTimedOut 
import time

sys.path.append(os.path.abspath("../../"))

from BioSANS2020.prepcodes.process import *
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.propagation.recalculate_globals import *

set_file = {}
"""
with open("console-result_summary.dat","r") as f:
	for row in f:
		rr = row.replace(".\\","").replace("-----------------------","").split()
		try:
			set_file[rr[0]] = float(rr[-1])
		except:
			pass
"""
		
cdirs = glob.glob("./0*")
start =  1

def run_this(ih):
	globals2.init()

	cdir = cdirs[ih]
	method = "rk4-3" #"Euler-3"#
	FileIn = "moles"
	molar = False
	relative = 0.0001
	absolute = 0.0001
	settings = open(cdir+"/"+cdir+"-settings.txt","r")
	NoAmountNoConc = 0
	variables = ""
	for x in settings:
		row = x.split(":")
		if row[0].strip() == "duration":
			tend = float(row[1])
		elif row[0].strip() == "steps":
			steps = int(row[1])
		elif row[0].strip() == "absolute":
			absolute = float(row[1])	
		elif row[0].strip() == "relative":
			relative = float(row[1])	
		elif row[0].strip() == "amount":
			if row[1].strip()=="":
				method = "rk4-2" #"Euler-2"#
				FileIn = "molar"
				molar = True
				NoAmountNoConc = NoAmountNoConc + 1
		elif row[0].strip() == "concentration":
			if row[1].strip()=="":
				method = "rk4-3" #"Euler-2"#
				FileIn = "moles"
				molar = False
				NoAmountNoConc = NoAmountNoConc + 1 
		elif row[0].strip() == "variables":
				variables = row[1].strip().split(",")
				variables = [x.strip() for x in variables]
			
	if NoAmountNoConc == 2:
		method = "rk4-3" #"Euler-3"#
		FileIn = "moles"
		molar = False
		
	#print(absolute, relative)
	if os.path.isfile(cdir+"/"+cdir+"-sbml-l3v1.xml"):
		sbml = cdir+"/"+cdir+"-sbml-l3v1.xml"
		sbml_to_topo(sbml,molar,variables)
		topo = cdir+"/"+cdir+"-sbml-l3v1.xml.topo"
	elif os.path.isfile(cdir+"/"+cdir+"-sbml-l2v4.xml"):
		sbml = cdir+"/"+cdir+"-sbml-l2v4.xml"
		sbml_to_topo(sbml,molar,variables)
		topo = cdir+"/"+cdir+"-sbml-l2v4.xml.topo"	
	elif os.path.isfile(cdir+"/"+cdir+"-sbml-l3v2.xml"):
		sbml = cdir+"/"+cdir+"-sbml-l3v2.xml"
		sbml_to_topo(sbml,molar,variables)
		topo = cdir+"/"+cdir+"-sbml-l3v2.xml.topo"

	outf = cdir+"/"+cdir
	topfile = open(topo,"r")
	Volume = 1
	for row in topfile:
		if row[0] == "#":
			gg = row.split(",")[1:]
			for x in gg:
				xx = [ g.strip() for g in x.split("=")  ] 
				if xx[0] == "Volume":
					Volume = xx[1]	
					break
			break
	
	options = [0.1,0.01,0.001,0.0001]
	sbml2 = sbml.split('.\\')[-1]
	if sbml2 in set_file:
		if set_file[sbml2] == 0.0001:
			options = [0.1,0.01,0.001,0.0001]
		else:
			options = [set_file[sbml2]]
	for delx in options:
		if os.path.isfile(topo):
			process(
				rfile    	= topo,
				miter		= 1,
				inMolar		= FileIn,
				Vm 			= Volume,
				tn			= tend ,
				delX		= delx,
				normalize	= False,
				logx		= False,
				logy		= False,
				method		= method,
				tlen		= steps,
				mix_plot	= True,
				save		= True,
				out_fname	= outf,
				plot_show	= False,
				Cinput		= {},
				vary 		= "",
				mult_proc	= False,
				implicit    = True,
				items		= None
			)
			pass
		
		dataOri = np.array(np.genfromtxt(outf+"-results.csv",delimiter=',',names=True))
		dataSim = np.array(np.genfromtxt(outf+"_"+method+".dat",delimiter='\t',names=True))	
		simLabel = dataSim.dtype.names
		oriLabel = dataOri.dtype.names
		
		dataOri = np.array( [ list(x) for x in dataOri ] )
		dataSim = np.array( [ list(x) for x in dataSim ] )
		
		idOrder = []
		for i_ori in range(len(oriLabel)):
			try:
				i_sim = simLabel.index(oriLabel[i_ori])
				idOrder.append(i_sim)
			except:
				i_sim = simLabel.index(oriLabel[i_ori].lower())
				idOrder.append(i_sim)
		dataSim = dataSim[:,idOrder]
		
		F = 1
		simD = dataSim[:,1:]
		simD[np.isnan(simD)] = 1000
		MeanTest = np.sum( np.abs(dataOri[:,1:] - simD) >=  absolute * F + relative * F * np.abs(dataOri[:,1:]),1)
		Wrong = np.sum(MeanTest>0)
		
		lenData = len(dataOri)
		if Wrong == 0:
			print(cdir+"-sbml-l3v1.xml","-----------------------correct---------------------",delx)
			out_file_summary.write(cdir+"-sbml-l3v1.xml"+"\t"+"-----------------------correct---------------------"+"\t"+str(delx)+"\n")
			return True	
			
	if Wrong > 0:
		print(cdir+"-sbml-l3v1.xml", Wrong/lenData,"-----------------------Wrong---------------------",delx)
		out_file_summary.write(cdir+"-sbml-l3v1.xml"+"\t"+str(Wrong/lenData)+"\t"+"-----------------------Wrong---------------------"+"\t"+str(delx)+"\n")
		return False
	return True	


out_file_summary = open("Newconsole-result_summary.dat","w") 

Wrong = 0
for ih in range(start-1,len(cdirs)):
	try:
	#if True:
		try: 
			sol = func_timeout(120,run_this, args=(ih,) ) 
			if sol == False:
				pass
				#break
		except FunctionTimedOut: 
			print("time out error")
	
	except Exception as e:
		print(e)
		out_file_summary.write(str(ih+1)+"-sbml-l3v1.xml"+"-------ERROR-------"+"\n")
		print(str(ih+1)+"-sbml-l3v1.xml"+"-------ERROR-------")
		pass
		#break
	
	#break
out_file_summary.close()