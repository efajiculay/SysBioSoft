import numpy
import gillespy2
from gillespy2.sbml import SBMLimport
import matplotlib.pyplot as plt
from gillespy2.solvers.numpy.tau_hybrid_solver import TauHybridSolver
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.numpy.tau_leaping_solver import TauLeapingSolver
from gillespy2.solvers.numpy import NumPySSASolver

# python3 -m pip install https://github.com/StochSS/GillesPy2/archive/main.zip --user --upgrade
		
model = SBMLimport.convert("00001-sbml-l2v1.xml")[0]
print(model)
model.timespan(numpy.linspace(0, 50, 100))
results = model.run(solver=SSACSolver, number_of_trajectories=10, t=50)


L = [k for k in results[0]]
for x in results:
	N = len(x[L[0]])
	#for i in range(N):
		#print([ x[k][i] for k in x ])
	[plt.plot(x[L[0]],x[k]) for k in L[1:] ]

plt.show()	
