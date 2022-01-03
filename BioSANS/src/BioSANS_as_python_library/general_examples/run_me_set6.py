import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

# Example 6


modelA = """
	FUNCTION_DEFINITIONS:
	Ao = 1
	Bo = 1
	kf1 = 1
	kb1 = 1

	#REACTIONS, Volume = 1, tend = 100, steps = 100, FileUnit = molar
	A <=> B, kf1, kb1


	@CONCENTRATION
	A, Ao
	B, Bo
"""


my_model = biosans.model(modelA).save_traj("AtoBrev_lambda.txt").plot()
data = my_model.run("LNA3")
my_model.clean()