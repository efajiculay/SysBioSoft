"""This example  make use of  the process  module  to  perform parameter
estimation given some topology and deterministic trajectories. Check the
topology files, those that will be estimate are with negative value.

The following are the list of argument in the process function.

Args:
	rfile (str): file name of BioSANS topology file.
	miter (int, optional): Number of iteration or trajectory samples
		for stochastic integration
	conc_unit (bool, optional): "mole","molar", or "molecules" - the
		unit used in any amount in topology file. Defaults to
		"molecules".
	v_volms (float, optional): the volume of compartment used in the
		simulation. Defaults to 1.0e-20.
	tend (float): trajectory simulation end time. Defaults to 1.
	del_coef (float, optional): factor for modifying time steps used
		in the integration/propagation of ODE. Defaults to 10.
	normalize (bool, optional): True  will  be normalized the y axis
		based on max value . Defaults to False.
	logx (bool, optional): If True, the x-axis will be in log scale.
		Defaults to False.
	logy (bool, optional): if True, the y-axis will be in log scale.
		Defaults to False.
	method (str, optional): Defaults to "CLE". Any of the option in
		the list of available method keywords is listed below;

		Stochastic (refer to section 10.2.4)

		1.	"CLE"            - Molecules(micro), tau-adaptive
		2.	"CLE2"           - Molecules(micro), cle-fixIntvl
		3.	"Gillespie_"     - Molecules(micro), Direct method
		4.	"Tau-leaping"    - Molecules(micro),
								Not swapping with Gillespie
		5.	"Tau-leaping2"   - Molecules(micro),
								Swapping with Gillespie
		6.	"Sim-TauLeap"    - Molecules(micro), Simplified,
								Swapping with Gillespie

		Deterministic (refer to section 10.2.1)

		7.	"Euler-1"        - Molecules(micro), tau-adaptive-1
		8.	"Euler-2"        - Molar (macro), tau-adaptive-1
		9.	"Euler-3"        - Mole (macro), tau-adaptive-1
		10.	"Euler2-1"	     - Molecules(micro), tau-adaptive-2
		11.	"Euler2-2"       - Molar (macro), tau-adaptive-2
		12.	"Euler2-3"       - Mole (macro), tau-adaptive-2
		13.	"ODE-1"          - Molecules(micro),
								using ode_int from scipy
		14.	"ODE-2"          - Molar(macro),
								using ode_int from scipy
		15.	"ODE-3"          - Mole(macro), using ode_int from scipy
		16.	"rk4-1"          - Molecules(micro), fix-interval
		17.	"rk4-2"          - Molar(macro), fix-interval
		18.	"rk4-3"          - Mole(macro), fix-interval
		19.	"rk4-1a"         - Molecules(micro), tau-adaptive
		20.	"rk4-2a"         - Molar(macro), tau-adaptive
		21.	"rk4-3a"         - Mole(macro), tau-adaptive

		Linear Noise Approximation (refer to 10.1.2 & 10.2.2)

		22.	"LNA"             - Numeric, values
		23.	"LNA-vs"          - Symbolic, values, Macroscopic
		24.	"LNA-ks"          - Symbolic, f(ks), Macroscopic
		25.	"LNA-xo"          - Symbolic, f(xo), Macroscopic
		26.	"LNA2"            - Symbolic, f(xo,ks), Microscopic
		27.	"LNA3"            - Symbolic, f(xo,ks), Macroscopic
		28.	"LNA(t)"          - COV-time-dependent, Macroscopic
		29.	"LNA2(t)"         - FF-time-dependent, Macroscopic

		Network Localization (refer to 10.1.3)

		30.	"NetLoc1"         - Symbolic, Macroscopic
		31.	"NetLoc2"         - Numeric, Macroscopic

		Parameter estimation (refer to 10.2.3)

		32.	"k_est1"          - MCEM, Macroscopic
		33.	"k_est2"          - MCEM, Microscopic
		34.	"k_est3"          - NM-Diff. Evol., Macroscopic
		35.	"k_est4"          - NM-Diff. Evol., Microscopic
		36.	"k_est5"          - Parameter slider/scanner
		37.	"k_est6"          - Nelder-Mead (NM), Macroscopic
		38.	"k_est7"          - Nelder-Mead (NM), Microscopic
		39.	"k_est8"          - Powell, Macroscopic
		40.	"k_est9"          - Powell, Microscopic
		41.	"k_est10"         - L-BFGS-B, Macroscopic
		42.	"k_est11"         - L-BFGS-B, Microscopic

		Symbolic/Analytical expression of species (refer to 10.1.1)

		43.	"Analyt"          - Pure Symbolic :f(t,xo,k)
		44.	"Analyt-ftx"      - Semi-Symbolic :f(t,xo)
		45.	"SAnalyt"         - Semi-Symbolic :f(t)
		46.	"SAnalyt-ftk"     - Semi-Symbolic :f(t,k)
		47.	"Analyt2"         - Creates commands for wxmaxima

	tlen (int, optional): number of integration steps reported in
		the final result. Defaults to 1000.
	mix_plot (bool, optional): If True, all species are plotted in
		one plot/figure. Defaults to True.
	save (bool, optional): If True, the resulting trajectory in the
		simulation will be saved as a file. Defaults to True.
	out_fname (str, optional): output filename. Defaults to "".
	plot_show (bool, optional): If True, an image of the plots will
		be created in the directory of the topology file.
	time_unit (str, optional): Defaults to "time (sec)".
	vary (str, optional): Varying initial concentration.
		Defaults to "".
	vary2 (str, optional): [description]. Defaults to "".
	mult_proc (bool, optional): If True, trajectories will be propa-
		gated on parallel. Defaults to False.
	implicit (bool, optional): True means report in time intervals
		similar to the input time intervals even if actual step is
		more or less. Defaults to False.
	items (tuple, optional): (canvas, scroll_x, scroll_y).
		Defaults to None.
	exp_data_file ([type], optional): Experimental data file contai-
		ning True or accepted trajectories. Defaults to None.
	c_input (dict, optional): [description]. Defaults to {}.

Returns:
	list: list of simulated trajecotry.
		data[j][0] - time for trajectory j
		data[j][1][:, i] - trajectories of each component i
"""

import sys
import os
from func_timeout import func_timeout, FunctionTimedOut
sys.path.append(os.path.abspath("../../"))
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.prepcodes.process import process

globals2.init(globals2)

if __name__ == '__main__':

    start = 1
    ends = 2
    for ih in range(start - 1, ends):
        proc_global.init(proc_global)

    #crucials = [14, 33, 34, 38, 44]
    # for val in crucials:
        #ih = val - 1
        # proc_global.init()

        # Uncomment and comment to try each
        method = "k_est6"   #NeldMead
        # method = "k_est8"   #Powell
        # method = "k_est10"  #L-BFGS-B
        # method = "k_est1"  # MCEM

        FileIn = "molar"  # uses macroscopic equations but still uses file units and no conversion

        topo = "Test" + str(ih + 1) + ".txt"
        EdataFile = "ode_res" + str(ih + 1) + "_ODE-2.dat"

        # get_globals(topo)
        with open(topo, "r") as file:
            rows = []
            last = ""
            for row in file:
                if last == "Function_Definitions":
                    if row.strip() != "" and row[0] != "#":
                        exec(row.strip(), globals())
                    elif row[0] == "#":
                        last = "#"
                elif row.strip() == "Function_Definitions:":
                    last = "Function_Definitions"
                else:
                    pass

        topfile = open(topo, "r")
        Volume = 1
        try:
            result = func_timeout(600, lambda:
                process(
                    rfile      = topo,
                    miter      = 1,
                    conc_unit  = FileIn,
                    v_volms    = Volume,
                    tend       = 50,
                    del_coef   = 1,
                    normalize  = False,
                    logx       = False,
                    logy       = False,
                    method     = method,
                    tlen       = 200,
                    mix_plot   = True,
                    save       = False,
                    out_fname  = None,
                    plot_show  = False,
                    c_input    = {},
                    vary       = "",
                    mult_proc  = False,
                    implicit   = True,
                    items      = 0,
                    exp_data_file = EdataFile
                )
            )
            print("\n")
            for x in result[0][1]:
                print(x, " = ", result[0][1][x])

            try:
                rs = result[0][1]

                AAD = 0
                for key in rs:
                    print(eval(key))
                    AAD = AAD + abs((eval(key) - rs[key]) / eval(key))
                AAD = AAD / len(rs)

                if AAD < 0.05:
                    print("Test" + str(ih + 1), " AAD =", AAD,
                          "-----------------------correct---------------------")
                    # break
                else:
                    print("Test" + str(ih + 1), " AAD =", AAD,
                          "-----------------------Wrong-----------------------")
                    # break
            except BaseException:
                print("Test", ih + 1, " have error ")

        except FunctionTimedOut:
            print("Time out error")
