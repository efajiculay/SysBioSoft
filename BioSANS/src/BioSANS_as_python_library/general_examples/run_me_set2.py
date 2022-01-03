import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

# Example 4

modelA = """
	#REACTIONS, Volume = 1, tend = 100, steps = 100, FileUnit = molar

	-Sa    => mA        , 100
	mA     => 0 phi     , 1
	mA     => A + mA    , 1
	A      => 0 phi     , 1

	-Sb    => mB        , 100
	mB     => 0 phi     , 1
	mB     => B + mB    , 1
	B      => 0 phi     , 1

	-Sc    => mC        , 100
	mC     => 0 phi     , 1
	mC     => C + mC    , 1
	C      => 0 phi     , 1  

	@CONCENTRATION

	mA    , 1
	mB    , 0
	mC    , 0

	A     , 0
	B     , 0
	C     , 0

	phi   , 0
	-Sa   , 0   ,lambda C : 1/(1+C**2)
	-Sb   , 0   ,lambda A : 1/(1+A**2)
	-Sc   , 0   ,lambda B : 1/(1+B**2)
"""
my_model = biosans.model(modelA).save_traj("repressilator.txt").plot()
data = my_model.run()
my_model.clean()
