import sys, os
from sympy import*
from func_timeout import func_timeout, FunctionTimedOut 
import matplotlib.pyplot as plt

from BioSANS2020.prepcodes.process import *
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.propagation.recalculate_globals import get_globals
globals2.init(globals2)

if __name__ == '__main__':

	process(
		rfile    	= "VGCN.dat",
		miter		= 10,
		conc_unit		= "molar",
		v_volms 			= 1.0e-19,
		tend			= 1 ,
		del_coef		= 1,
		normalize	= False,
		logx		= True,
		logy		= False,
		method		= "ODE-2",
		tlen		= 10000,
		mix_plot	= True,
		save		= True,
		out_fname	= "VCGN-result.txt",
		plot_show	= True,
		c_input		= {},
		vary 		= "",
		mult_proc	= False,
		implicit    = False,
		items		= 0,
		exp_data_file = None
	)	

	print("Propagation is done: Look at the files in folder")