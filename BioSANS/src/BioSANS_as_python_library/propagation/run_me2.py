import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

my_model = biosans.model(topo="VGCN.dat", Volume=1.0e-19, tend=1, steps=10000, save=False, out_fname="VCGN-result.txt")

data = my_model.run("ODE-2")

print(data)

