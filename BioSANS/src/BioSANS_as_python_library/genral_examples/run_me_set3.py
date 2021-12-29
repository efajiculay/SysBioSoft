import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020.myglobal import proc_global
from BioSANS2020 import biosans_lib as biosans

# Example 5

modelA = """
	#REACTIONS
	A <=> B, 0.5, 0.3

	@CONCENTRATION
	A, 100
	B, 10
"""


if __name__ == '__main__':
	my_model = biosans.model(modelA, FileIn="molecules").save_traj("AtoBrev.txt").plot()
	proc_global.init(proc_global)
	data = my_model.run(method="Gillespie_", ntraj=30, mult_proc=True)
	my_model.clean()