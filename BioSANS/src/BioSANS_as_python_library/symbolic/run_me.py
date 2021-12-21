import sys, os
from sympy import*
from func_timeout import func_timeout, FunctionTimedOut 
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("../../"))

from BioSANS2020.prepcodes.process import *
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.propagation.recalculate_globals import get_globals

Wrong = 0
start = 1
for ih in range(start-1,2):
	globals2.init(globals2)

	#method = "Analyt"
	method = "Analyt-ftx"
	#method = "ODE-2"
	FileIn = "molar" #uses macroscopic equations but still uses file units and no conversion
	relative = 0.0001
	absolute = 0.0001


	topo = "Test"+str(ih+1)+".txt"					
	outf = "symbolic_res"+str(ih+1)
	
	with open(topo,"r") as file:
		rows = []
		last = ""
		for row in file:
			if last == "Function_Definitions":
				if row.strip()!="" and row[0] != "#":
					exec(row.strip(),globals())
				elif row[0] == "#":
					last = "#"
			elif row.strip() == "Function_Definitions:":
				last = "Function_Definitions"
			else:
				pass
	
	Volume = 1
	
	try:
		result = func_timeout(600, lambda :  
			process(
				rfile    	= topo,
				miter		= 1,
				conc_unit		= FileIn,
				v_volms 			= Volume,
				tend			= 100 ,
				del_coef		= 1,
				normalize	= False,
				logx		= False,
				logy		= False,
				method		= method,
				tlen		= 200,
				mix_plot	= True,
				save		= True,
				out_fname	= outf,
				plot_show	= False,
				c_input		= {},
				vary 		= "",
				mult_proc	= False,
				implicit    = True,
				items		= 0
			)		
		)	
		for x in result[0][1]:
			print(x," = ",result[0][1][x])
		print("\n")
	except FunctionTimedOut: 
		pass		
	
	ode_res = "ode_res"+str(ih+1)+"_ODE-2.dat"
	dataOri = np.array(np.genfromtxt(ode_res,delimiter='\t',names=True))
	oriLabel = dataOri.dtype.names
	dataOri = np.array( [ list(x) for x in dataOri ] )
	time = dataOri[:,0]
	
	#if True:
	try:
		rs = result[0][1]
		#print(rs)
		dataSim = []
		for tx in time:
			t = tx
			row = [t]
			for x in oriLabel[1:]:
				key = x+"(t)"
				row.append( eval( rs[key] ) )
			dataSim.append(row)		
		dataSim = np.array(dataSim)	
		
		Wrong = 0
		F = 1
		for i_ori in range(1,len(oriLabel)):
			#Wrong = Wrong + np.sum(np.abs(dataOri[1:,i_ori]-dataSim[1:,i_ori]) >=  absolute * F + relative * F * np.abs(dataOri[1:,i_ori]))
			Wrong = Wrong + np.sum( (np.abs(dataOri[1:,i_ori]-dataSim[1:,i_ori])/dataOri[1:,i_ori])/len(dataOri[1:,i_ori])>0.05 )              #sum of time points with relative error > 5%
		Wrong = Wrong/len(oriLabel)                                                                                                            #average of wrong time points in all species
		
		if Wrong < 0.05*len(dataOri):                                                                                                          #check if wrong time points is less than 5% of data points
			print("Test"+str(ih+1),"-----------------------correct---------------------")
			#break	
		if Wrong >= 0.05*len(dataOri):
			print("Test"+str(ih+1), Wrong,"-----------------------Wrong---------------------")
			#break
	except:
		print("Test",ih+1," have error ")
		#break