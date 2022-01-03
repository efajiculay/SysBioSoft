import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

# Example 5

modelA = """
	#REACTIONS, Volume = 1, tend = 100, steps = 100, FileUnit = molar
	A <=> B, 1, 1


	@CONCENTRATION
	A, 0
	B, 0
"""


my_model = biosans.model(modelA).save_traj("AtoBrev_lambda.txt").plot()

data = my_model.run("Analyt")
data = my_model.run("Analyt-ftx")
my_model.clean()


modelA = """
	FUNCTION_DEFINITIONS:
	Ao = 0
	Bo = 0
	kf1 = 0
	kb1 = 0

	#REACTIONS, Volume = 1, tend = 100, steps = 100, FileUnit = molar
	A <=> B, kf1, kb1


	@CONCENTRATION
	A, Ao
	B, Bo
"""


my_model = biosans.model(modelA).save_traj("AtoBrev_lambda.txt").plot()

data = my_model.run("Analyt")
data = my_model.run("Analyt-ftx")
my_model.clean()