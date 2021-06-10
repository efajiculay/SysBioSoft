import sys, os
from sympy import*
from func_timeout import func_timeout, FunctionTimedOut 
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("../../"))

from BioSANS2020.prepcodes.process import *
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.propagation.recalculate_globals import *
globals2.init()

if __name__ == '__main__':

	start = 1
	ends = 2
	for ih in range(start-1,ends):
		proc_global.init()

	#crucials = [14, 33, 34, 38, 44]
	#for val in crucials:
		#ih = val - 1
		#proc_global.init()	
		
		#Uncomment and comment to try each
		#method = "k_est6"   #NeldMead
		#method = "k_est8"   #Powell
		#method = "k_est10"  #L-BFGS-B
		method = "k_est1"    #MCEM
 		
		FileIn = "molar" #uses macroscopic equations but still uses file units and no conversion

		topo = "Test"+str(ih+1)+".txt"		
		EdataFile = "ode_res"+str(ih+1)+"_ODE-2.dat"

		#get_globals(topo)
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
		
		topfile = open(topo,"r")
		Volume = 1
		try:
			result = func_timeout(600, lambda :  
				process(
					rfile    	= topo,
					miter		= 1,
					inMolar		= FileIn,
					Vm 			= Volume,
					tn			= 50 ,
					delX		= 1,
					normalize	= False,
					logx		= False,
					logy		= False,
					method		= method,
					tlen		= 200,
					mix_plot	= True,
					save		= False,
					out_fname	= None,
					plot_show	= False,
					Cinput		= {},
					vary 		= "",
					mult_proc	= False,
					implicit    = True,
					items		= 0,
					expDataFile = EdataFile
				)		
			)	
			print("\n")
			for x in result[0][1]:
				print(x," = ",result[0][1][x])
				
			try:
				rs = result[0][1]
				
				AAD = 0
				for key in rs:
					print(eval(key))
					AAD = AAD + abs((eval(key)-rs[key])/eval(key))
				AAD = AAD/len(rs)
				
				if AAD < 0.05:                                                                                                 
					print("Test"+str(ih+1)," AAD =",AAD,"-----------------------correct---------------------")
					#break	
				else: 
					print("Test"+str(ih+1)," AAD =",AAD,"-----------------------Wrong-----------------------")
					#break
			except:
				print("Test",ih+1," have error ")				
				
		except FunctionTimedOut: 
			print("Time out error")			