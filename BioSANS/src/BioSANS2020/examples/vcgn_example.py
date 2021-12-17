"""

An example usage of BioSANS as a python import

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.prepcodes.process import process

# methods = ["CLE", "tau_leaping", "Euler", "LNA", "Gillespie"]
# process(rfile="Reactions.dat",miter=10,conc_unit=True,v_volms = 1.0e-20,tend=1,
# del_coef=10,normalize=False,logx=False,logy=False,method="CLE",tlen=1000,
# save=True,plot_fname="")


process(
    rfile="VGCN/VGCN.dat",
    miter=1,
    logx=True,
    v_volms=1.0e-20,
    conc_unit=True,
    tend=1,
    del_coef=10,
    # normalize=False,
    method="CLE",
    # tlen=100000,
    mix_plot=True,
    out_fname="VGCN/Faji"
)
