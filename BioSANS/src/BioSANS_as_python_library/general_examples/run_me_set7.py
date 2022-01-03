import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

# Example 7


modelA = """
	Function_Definitions:
	Ao = 100
	Bo = 0
	kf1 = 0.5

	#REACTIONS
	A => B, -1

	@CONCENTRATION
	A, Ao
	B, Bo
"""

dataA = """
	time,A,B
	0.0,100.0,0.0
	5.0,8.208499847123571,91.79150015287648
	10.0,0.6737946943694972,99.32620530563057
	50.0,1.4844291096420276e-09,99.99999999851562
"""

if __name__ == '__main__':
	my_model = biosans.model(modelA).data(dataA)
	data = my_model.run("k_est6")
	my_model.clean()