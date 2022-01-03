import sys, os
sys.path.append(os.path.abspath("../../"))
from BioSANS2020 import biosans_lib as biosans
from BioSANS2020.myglobal import proc_global    # for "k_est1" or MCEM

# Example 8

modelA = """
    Function_Definitions:
    Bo = 0

    #REACTIONS
    A => B, -1
    B => C, -1

    @CONCENTRATION
    A, -1
    B, Bo
    C, -1
"""

dataA = """
    time,B
    0,0
    1,33.57189029
    2,45.23304872
    3,45.85987493
    4,41.46473222
    5,35.26129043
    10,10.76228036
    15,2.638978041
    20,0.608338062
    25,0.137339421
    30,0.03077597
    35,0.006877835
    40,0.00153554
    45,0.000342703
    50,7.65E-05
"""

if __name__ == '__main__':
    proc_global.init(proc_global)               # for "k_est1" or MCEM
    my_model = biosans.model(modelA).data(dataA)
    data = my_model.run("k_est1")
    my_model.clean()