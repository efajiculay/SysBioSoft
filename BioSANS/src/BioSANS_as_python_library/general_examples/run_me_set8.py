import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

# Example 8


modelA = """
	Function_Definitions:
	Ao = 100
	Co = 0

	#REACTIONS
	A => B, -1
	B => C, -1

	@CONCENTRATION
	A, Ao
	B, -1
	C, Co
"""

dataA = """
	time,A,C
	0,100,0
	1,60.65306595,5.775043762
	2,36.78794412,17.97900716
	3,22.313016,31.82710907
	4,13.5335283,45.00173948
	5,8.208499842,56.53020973
	10,0.673794694,88.56392495
	15,0.055308438,97.30571352
	20,0.004539994,99.38712194
	25,0.000372665,99.86228791
	30,3.06E-05,99.96919344
	35,2.51E-06,99.99311965
	40,2.06E-07,99.99846425
	45.5,1.32E-08,99.99970501
	50,1.38E-09,99.99992352
"""

if __name__ == '__main__':
	my_model = biosans.model(modelA).data(dataA)
	data = my_model.run("k_est6")
	my_model.clean()