import numpy
import gillespy2
from gillespy2.sbml import SBMLimport
import matplotlib.pyplot as plt
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver

# python3 -m pip install https://github.com/StochSS/GillesPy2/archive/main.zip --user --upgrade
		
model = SBMLimport.convert("00002-sbml-l3v1.xml")[0]
results = model.run(solver=BasicTauHybridSolver, number_of_trajectories=10)

for x in results:
	for k in x:
		print(k)