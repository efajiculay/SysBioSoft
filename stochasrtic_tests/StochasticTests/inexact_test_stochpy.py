import sys, os
import glob
from func_timeout import func_timeout, FunctionTimedOut 
import stochpy
import matplotlib.pyplot as plt
import numpy as np

MainDIR = "C:/Users/efaji/Desktop/working/BIOSANS_WORK/StochasticTests/"


if __name__ == '__main__':
	dirs = glob.glob(MainDIR+"/0*")
	smod = stochpy.SSA()

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

	algMet = "stochpy-tauleaping"
	TestMean = True
	TestStdv = True
	max_iter = 100

	start =  1

	label = ("Mean" if TestMean else "") + ("Stdv" if TestStdv else "")+"_"+str(max_iter)+"_"+algMet+"_inexact"
	#label = "MeanStdv_SSA_exact_tested"
	assess = lambda x : "Passed" if x else "Failed"

	test_setting = "stochpy"
	outfile = open(MainDIR+label+"simulation_summary_"+algMet+"_"+str(max_iter)+".out","w")

	Passed = True
	for ih in range(start-1,len(dirs)):
	#for ihx in [21]:#, 22, 23, 24, 25, 26]: 
		#ih = ihx - 1
		try:
			dirW = dirs[ih]
			dirW = dirW.split("\\")[-1]
			iters = max_iter
			method = test_setting
			FileIn = "molecules"
			molar = False
			relative = 0.0001
			absolute = 0.0001
			settings = open(MainDIR+dirW+"/"+dirW+"-settings.txt","r")
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

			
			sbml = dirW+"-sbml-l2v4.xml"						
			outf = dirW+"/"+dirW+"_"+label
			outTrue = dirW+"/"+dirW
			Volume = 1
			
			out_file = open(MainDIR+dirW+"/"+dirW+sbml+"_simulation_summary"+label+".out","w")
			for scaler in [1]:
				
				hdir = MainDIR+dirW
				stochpy.SBML2PSC(sbmlfile=sbml,sbmldir=hdir, pscfile=None, pscdir=None)
				
				smod.Model(sbml+".psc")
				smod.DoStochSim(method="Tauleap",end=steps,mode="time",trajectories=max_iter,quiet=True) 

				smod.GetRegularGrid()
				time = smod.data_stochsim_grid.time  
				data = smod.data_stochsim_grid.species
				simLabel = smod.data_stochsim.species_labels
				
				SfileN = open(MainDIR+outf+"_"+method+".dat","w")
				SfileN.write("time"+"\t"+"\t".join(simLabel)+"\n")
				for x in range(max_iter):
					for j in range(len(time)):
						SfileN.write(str(time[j])+"\t"+"\t".join([str(data[i][x][j]) for i in range(len(data))])+"\n")
				SfileN.close()
			
				dataSim = np.array(np.genfromtxt(MainDIR+outf+"_"+method+".dat",delimiter='\t',names=True))	
				simLabel = dataSim.dtype.names
				dataSim = np.array( [ list(x) for x in dataSim ] )
				#print(dataSim)

				dataOri = np.array(np.genfromtxt(MainDIR+outTrue+"-results.csv",delimiter=',',names=True))
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
				#print(dataSimSum)
				dataSimMean = np.mean(dataSimSum,0)
				dataSimStds = np.sqrt(np.var(dataSimSum,0))
				#print(dataSimMean,"\n")
				#print(dataSimStds)

				last = len(dataSimMean[0])
				meanTest = dataSimMean[:,1:last] / dataOri[:,1:last]
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
				stdsTest = np.sqrt(dataSimVar)/dataOri[:,last:]
				stdsTest[np.isnan(stdsTest)] = 1
				#print(stdsTest)

				ll = last-1
				PassedMean = True
				PassedStdv = True
				MRange = [0.95, 1.05]
				SRange = [0.95, 1.05]
				Mpr = 0
				Spr = 0
				for ig in range(1,len(meanTest)):
					if TestMean:
						if sum(meanTest[ig]>=MRange[0])==ll and sum(meanTest[ig]<=MRange[1]) == ll:
							#PassedMean = True
							pass
						else:
							PassedMean  = False
							Mpr = Mpr + 1
					if TestStdv:
						if sum(stdsTest[ig]>=SRange[0])==ll and sum(stdsTest[ig]<=SRange[1])==ll:
							#PassedStdv = True
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
		except:
			pass
		break
	outfile.close()