import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans

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

my_model = biosans.model(vcgn).save_traj("vcgn.txt").plot(logx=True)
data = my_model.run(method="Euler-2",tend=1,step_size_scaler=0.001,steps=1000,implicit=False)
data = my_model.run(method="ODE-2",tend=1,steps=1000)
my_model.clean()


# Example 2

modelA = """
	#REACTIONS
	A <=> B, 0.5, 0.3

	@CONCENTRATION
	A, 100
	B, 10
"""
my_model = biosans.model(modelA).save_traj("AtoBrev.txt").plot()
data = my_model.run()
my_model.clean()


# Example 3

modelA = """
	Function_Definitions:
	alp = 0.15
	rx  = 0.12
	ry  = 0.7
	kx  = 0.5
	ky  = 0.1
	Yo  = 0.7
	Zo  = 0
	theta = lambda x,k,n : (x**n)/(x**n+k**n+sin(Zo))

	#REACTIONS
	0 NONE =>  X, 0
	0 NONE =>  timer_X, 1

	0 NONE =>  Y, 0  ::::: lambda X, Y    : -alp*Y + rx*theta(X,kx,2)
	0 NONE =>  Z, 0  ::::: lambda X, Y, Z : -alp*Z + rx*theta(X,kx,2)*ry*theta(Y,ky,2)


	@CONCENTRATION
	timer_X,   0, lambda timer_X  :  0 if timer_X >= 32 else timer_X
	X,         1, lambda timer_X  :  0 if timer_X >= 16 else 1
	Y,         Yo
	Z,         Zo
	NONE, 0 
"""
my_model = biosans.model(modelA).save_traj("feed_forward_loop.txt").plot()
data = my_model.run()
my_model.clean()