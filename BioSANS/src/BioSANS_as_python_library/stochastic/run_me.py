import sys, os
import glob
from func_timeout import func_timeout, FunctionTimedOut 
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("../../"))

from BioSANS2020.prepcodes.process import *
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.propagation.recalculate_globals import get_globals
globals2.init(globals2)


if __name__ == '__main__':
	dirs = glob.glob("./0*")

	#start: 0
	#duration: 50
	#steps: 50
	#variables: X
	#absolute: 
	#relative: 
	#amount: X
	#concentration: 
	#output: X-mean, X-sd
	#meanRange: (-3, 3)
	#sdRange: (-5, 5)

	algMet = "Gillespie"
	TestMean = True
	TestStdv = True
	max_iter = 100

	start =  1
	ends = 2

	label = ("Mean" if TestMean else "") + ("Stdv" if TestStdv else "")+"_"+str(max_iter)+"_"+algMet+"_exact2"
	label = "MeanStdv_SSA_exact_testedtry"
	assess = lambda x : "Passed" if x else "Failed"

	test_setting = {
		"CLE"         : ["CLE",[1]],#200,100,50,10,5,3,
		"CLE2"        : ["CLE2",[1,0.1,0.001]],
		"Gillespie"   : ["Gillespie_",[1,1]],
		"Tau-leaping" : ["Tau-leaping",[1,0.2/0.3]],
		"Sim-TauLeap" : ["Sim-TauLeap",[1,0.2/0.3]],
		"Runge-kutta" : ["rk4-2",[0.001]],
	}
	outfile = open(label+"simulation_summary_"+algMet+"_"+str(max_iter)+".out","w")

	Passed = True
	for ih in range(start-1,ends):
		proc_global.init(proc_global)
		globals2.init(globals2)

		dir = dirs[ih]
		iters = max_iter
		method = test_setting[algMet][0]
		FileIn = "molecules"
		molar = False
		relative = 0.0001
		absolute = 0.0001
		settings = open(dir+"/"+dir+"-settings.txt","r")
		NoAmountNoConc = 0
		for x in settings:
			row = x.split(":")
			if row[0].strip() == "duration":
				tend = float(row[1])
			elif row[0].strip() == "steps":
				steps = float(row[1])
			elif row[0].strip() == "meanRange":
				MRange = eval(row[1])	
			elif row[0].strip() == "sdRange":
				SRange = eval(row[1])	

		sbml = dir+"/"+dir+"-sbml-l3v1.xml"
		if os.path.isfile(sbml):
			sbml_to_topo(sbml,molar)
			topo = dir+"/"+dir+"-sbml-l3v1.xml.topo"
		else:
			sbml = dir+"/"+dir+"-sbml-l2v4.xml"
			sbml_to_topo(sbml,molar)
			topo = dir+"/"+dir+"-sbml-l2v4.xml.topo"				
			
		outf = dir+"/"+dir+"_"+label
		outTrue = dir+"/"+dir
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
		topfile.close()
		
		out_file = open(sbml+"_simulation_summary"+label+".out","w")
		for scaler in test_setting[algMet][1]:
			if os.path.isfile(topo):
				
				process(
					rfile		= topo,
					miter		= iters,
					conc_unit		= FileIn,
					v_volms 			= Volume,
					tend			= tend ,
					del_coef		= scaler,
					normalize	= False,
					logx		= False,
					logy		= False,
					method		= method,
					tlen		= steps,
					mix_plot	= True,
					save		= True,
					out_fname	= outf,
					plot_show	= False,
					c_input		= {},
					vary 		= "",
					mult_proc	= True,
					implicit	= True,
					items		= None
				)
				
				dataSim = np.array(np.genfromtxt(outf+"_"+method+".dat",delimiter='\t',names=True))	
				simLabel = dataSim.dtype.names
				dataSim = np.array( [ list(x) for x in dataSim ] )

				dataOri = np.array(np.genfromtxt(outTrue+"-results.csv",delimiter=',',names=True))
				oriLabel = dataOri.dtype.names
				dataOri = np.array( [ list(x) for x in dataOri ] )
				#print(simLabel)
				#print(oriLabel)
				sim_order = [0]
				last_index = 1+int((len(oriLabel)-1)/2)
				for i_ori in range(1,last_index):
					key = oriLabel[i_ori].replace("mean","").replace("sd","").replace("-","")
					i_sim = simLabel.index(key)
					sim_order.append(i_sim)
				#print(sim_order)
				dataSim = dataSim[:,sim_order]

				dataSimSum = []
				for ig in range(1,iters):
					start = ig*int(steps+1)
					mends = start+int(steps+1)
					dataSimSum.append(dataSim[start:mends,:])
				dataSimSum  = np.array(dataSimSum)
				dataSimMean = np.mean(dataSimSum,0)
				dataSimStds = np.sqrt(np.var(dataSimSum,0))
				#print(dataSimMean,"\n")
				#print(dataSimStds)

				last = len(dataSimMean[0])
				meanTest = np.sqrt(float(iters))*(dataSimMean[:,1:last] - dataOri[:,1:last])/dataOri[:,last:]
				meanTest[np.isnan(meanTest)] = 1
				#print(meanTest,"\n")
							
				dataSimVars = []
				for ig in range(1,iters):
					start = ig*int(steps+1)
					mends = start+int(steps+1)
					dataSimVars.append(dataSim[start:mends,1:last]-dataOri[:,1:last])                     
				dataSimVar = (1/float(iters))*np.sum(np.power(dataSimVars,2),0)
				#print(dataSimVar)
				#print(np.power(div,2))
				stdsTest = np.sqrt(float(iters)/2)*(dataSimVar/np.power(dataOri[:,last:],2)-1)
				stdsTest[np.isnan(stdsTest)] = 1
				#print(stdsTest)

				ll = last-1
				PassedMean = True
				PassedStdv = True
				Mpr = 0
				Spr = 0
				for ig in range(1,len(meanTest)):
					if TestMean:
						if sum(meanTest[ig]>=MRange[0])==ll and sum(meanTest[ig]<=MRange[1]) == ll:
							pass
						else:
							PassedMean  = False
							Mpr = Mpr + 1
					if TestStdv:
						if sum(stdsTest[ig]>=SRange[0])==ll and sum(stdsTest[ig]<=SRange[1])==ll:
							pass
						else:
							PassedStdv  = False
							Spr = Spr + 1

				if PassedMean and PassedStdv:
					print("Passed -",sbml," mean = "+assess(PassedMean)+" stdev = "+assess(PassedStdv))
					outfile.write("Passed -"+sbml+" mean = "+assess(PassedMean)+" stdev = "+assess(PassedStdv)+"\n")
					break
		if not PassedMean or not PassedStdv :
			Mpr = 100*Mpr/len(meanTest)
			Spr = 100*Spr/len(meanTest)
			print("Failed -",sbml," mean = "+assess(PassedMean)+" stdev = "+assess(PassedStdv)+" % W-mean = "+str(Mpr)+" % W-stdv = "+str(Spr))
			outfile.write("Failed -"+sbml+" mean = "+assess(PassedMean)+" stdev = "+assess(PassedStdv)+" % W-mean = "+str(Mpr)+" % W-stdv = "+str(Spr)+"\n")
		out_file.write("\n")
		out_file.write(str(dataSimMean)+"\n")
		out_file.write("\n")
		out_file.write(str(dataSimStds)+"\n")
		out_file.close()
		#if not PassedMean or not PassedStdv :
			#break
	outfile.close()