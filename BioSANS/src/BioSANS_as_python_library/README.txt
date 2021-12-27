The following set of commands will update BioSANS by first uninstall old version followed by reinstallation. 
The updated or latest version in test.pypi.org will be installed.

	pip uninstall BioSANS2020
	pip install BioSANS2020 --upgrade
	
To run BioSANS, type any of the followig command and hit enter.

	BioSANS
	python -m BioSANS2020

How to use BioSANS as a library (python import)?
The following are commonly importe dlibraries

	import sys, os
	from sympy import*
	from func_timeout import func_timeout, FunctionTimedOut 
	import matplotlib.pyplot as plt

	from BioSANS2020.prepcodes.process import *
	from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
	from BioSANS2020.myglobal import mglobals as globals2
	from BioSANS2020.myglobal import proc_global as proc_global
	from BioSANS2020.propagation.recalculate_globals import get_globals
	globals2.init()
	
The following folders contains the examples codes that makes use of BioSANS as import

	run_me.py - "Contains the code that import BioSANS2020"

	parameter_estimation               # Performing parameter_estimation test and priting result
		ode_res1_ODE-1                 # sample given data
		ode_res1_ODE-2                 # sample given data
		Test1                          # test-case1
		Test2                          # test-case2
		
	semantic                           # Performing semantic test and priting result - the output is inside the numbered folder (e.g. 00001)
		Newconsole-result_summary      # Summary of test
		00001                          # test-case1 - contains SBML files, settings, and expected result, our result is named 00001_rk4-3.dat
		00002                          # test-case2 - contains SBML files, settings, and expected result, our result is named 00002_rk4-3.dat
	
	stochastic					       # Performing stochastic test and priting result - the output is inside the numbered folder (e.g. 00001)
		00001                          # test-case1 - contains SBML files, settings, and expected result, our result is named 00001_MeanStdv_SSA_exact_testedtry_Gillespie_.dat
		00002                          # test-case1 - contains SBML files, settings, and expected result, our result is named 00002_MeanStdv_SSA_exact_testedtry_Gillespie_.dat
		
	symbolic                           # Performing symbolic test and priting result
		ode_res1_ODE-1                 # sample given data
		ode_res1_ODE-2                 # sample given data
		Test1                          # test-case1
		Test2                          # test-case2
		
	propagation
		VGCN                           # topology file