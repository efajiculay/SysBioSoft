import sys, os
from sympy import*
from func_timeout import func_timeout, FunctionTimedOut 
import matplotlib.pyplot as plt

from BioSANS2020.prepcodes.process import *
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.propagation.recalculate_globals import *
globals2.init()

if __name__ == '__main__':

	process(
		rfile    	= "VGCN.dat",
		miter		= 10,
		inMolar		= "molar",
		Vm 			= 1.0e-19,
		tn			= 1 ,
		delX		= 1,
		normalize	= False,
		logx		= True,
		logy		= False,
		method		= "ODE-2",
		tlen		= 10000,
		mix_plot	= True,
		save		= True,
		out_fname	= "VCGN-result.txt",
		plot_show	= True,
		Cinput		= {},
		vary 		= "",
		mult_proc	= False,
		implicit    = False,
		items		= 0,
		expDataFile = None
	)	

	print("Propagation is done: Look at the files in folder")