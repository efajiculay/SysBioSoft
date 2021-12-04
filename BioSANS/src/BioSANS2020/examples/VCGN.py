#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.prepcodes.process import *

#methods = ["CLE", "Tau_leaping", "Euler", "LNA", "Gillespie"]
#process(rfile="Reactions.dat",miter=10,inMolar=True,Vm = 1.0e-20,tn=1,delX=10,normalize=False,logx=False,logy=False,method="CLE",tlen=1000,save=True,plot_fname="")

"""
process(
	rfile="VGCN/VGCN.dat",
	miter=1,
	logx=True,
	#Vm=1.0e-20,
	inMolar=True,
	tn=1,
	#delX=10,
	#normalize=False
	method="odeint",
	tlen=100000,
	mix_plot=True,
	out_fname="VGCN/Faji"
)
"""

process(
    rfile="VGCN/VGCN.dat",
    miter=1,
    logx=True,
    Vm=1.0e-20,
    inMolar=True,
    tn=1,
    delX=10,
    # normalize=False,
    method="CLE",
    # tlen=100000,
    mix_plot=True,
    out_fname="VGCN/Faji"
)
