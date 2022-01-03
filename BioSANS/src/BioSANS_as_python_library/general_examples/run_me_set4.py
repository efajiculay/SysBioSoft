import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

# Example 6

modelA = """
	FUNCTION_DEFINITIONS:
	Ao = 100
	Bo = 10
	kf1 = 1.5
	kb1 = 0.3

	#REACTIONS, Volume = 1, tend = 100, steps = 100, FileUnit = molar
	A <=> B, 1, 1 ::::: lambda A : kf1*A ::: lambda B : kb1*B

	@CONCENTRATION
	A, Ao
	B, Bo
"""


my_model = biosans.model(modelA).save_traj("AtoBrev_lambda.txt").plot()
data = my_model.run()
my_model.clean()