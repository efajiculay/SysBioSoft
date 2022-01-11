import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans
import numpy as np

################### PROPAGATION EXAMPLE #######################################

# Example 1

vcgn = """
	#REACTIONS
		X <=> 2 Y, 512, 512
	X + Y <=> Z  ,   8,   8
	Y + Z <=> X  ,   1,   1

	@CONCENTRATION
	X, 0.5
	Y, 0.5
	Z, 0.5
"""

my_model = biosans.model(vcgn)

print()
S, w, x, k = my_model.get_model_data()
print(S)
print()
print(w)
print()
print(x)
print()
print(k)
print()
print(np.dot(S,w))

print()
S, w, x, k = my_model.get_model_data(algeb=True)
print(S)
print()
print(w)
print()
print(x)
print()
print(k)
print()
print(np.dot(S,w))

my_model.clean()