"""

                 This module is the processes_hub module

The purpose of this module is to direct the information from the process
module into other modules. This module expects some to get some response
or data from thos emodule where its directs the information. The respon-
ses will be saved or plotted depending on the options provided to this
module.

The following are the list of function in this module.

1. process_hub


"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib._color_data as mcd

# from BioSANS2020.propagation.propensity import propensity_vec, \
#     propensity_vec_molar
from BioSANS2020.propagation.law_of_localization import law_loc_symbolic
from BioSANS2020.propagation.create_wxmaxima_command import for_wxmaxima
from BioSANS2020.propagation.stochastic.mystiffcle import \
    cle_calculate, cle2_calculate
from BioSANS2020.propagation.stochastic.gillespie_ssa import gillespie_ssa
from BioSANS2020.propagation.stochastic.tau_leaping import tau_leaping
from BioSANS2020.propagation.stochastic.tau_leaping2 import tau_leaping2
from BioSANS2020.propagation.stochastic.mytauleap import sim_tauLeap
# from BioSANS2020.propagation.stochastic.sde_int import SDE_int
from BioSANS2020.propagation.deterministic.euler_mod import euler_int, \
    euler2_int
# from BioSANS2020.propagation.deterministic.lna_approx import *
from BioSANS2020.propagation.deterministic.ode_int import ode_int
from BioSANS2020.propagation.deterministic.lna_function_of_time import \
    lna_non_steady_state, lna_non_steady_state2
from BioSANS2020.propagation.deterministic.runge_kutta4 \
    import rungek4_int, rungek4a_int
from BioSANS2020.propagation.symbolic.lna_approx2 import lna_symbolic2
from BioSANS2020.propagation.symbolic.analytical_sol import analyt_soln
from BioSANS2020.model.param_est.param_estimate import param_estimate

from BioSANS2020.model.fileconvert.convtopotosbml import topo_to_sbml
from BioSANS2020.model.param_est.param_slider import param_ode_int

from BioSANS2020.myglobal import mglobals as globals2
# from BioSANS2020.myglobal import proc_global

from BioSANS2020.gui_functs.draw_figure import draw_figure


def process_hub(
        time_var, sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
        v_volms=1, miter=1, logx=False, logy=False, del_coef=10,
        normalize=False, method="CLE", mix_plot=True,
        save=True, out_fname="", plot_show=True,
        time_unit="time (sec)", vary="", mult_proc=False,
        items=None, vary2="", implicit=False, rfile="",
        exp_data_file=None
):
    """This function  redirect all information from the process to other
    modules and expect the return from those  module for plotting or for
    saving into a file.

    Args:
        time_var (list): list of time points in the simulation
        sp_comp (dict): dictionary of appearance or position of species
            or component in the reaction tag of BioSANS topology file.

            For example;

                #REACTIONS
                A => B, -kf1    # negative means to be estimated
                B => C, kf2

            The value of sp_comp is

                sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}

                A appears in first reaction with index 0
                B appears in first and second reaction with index 0, 1
                C appears in second reaction with index 1
        ksn_dict (dict): dictionary of rate constant that appear in each
            reactions.

            For example;

                #REACTIONS
                A => B , 0.3        # first reaction
                B <=> C, 0.1, 0.2   # second reaction

            The value of ks_dict is

                ks_dict = {
                    0 : [0.3],      # first reaction
                    1 : [0.1, 0.2]  # second reaction
                }
        conc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation
        r_dict (dict): dictionary of reactant stoichiometry. For example

            r_dict = {
                0: {'A': 1},  # first reaction, coefficient of A is 1
                1: {'B': 1}   # second reaction, coefficient of B is 1
            }
        p_dict (dict): dictionary of product stoichiometry. For example

            p_dict = {
                0: {'B': 1},  # first reaction, coefficient of B is 1
                1: {'C': 1}   # second reaction, coefficient of C is 1
            }
        stoch_var (numpy.ndarray): stoichiometric matrix. For example

            stoch_var = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])
        v_volms (float, optional): volume of compartment. Defaults to 1.
        miter (int, optional): Number of iteration or trajectory samples
            for stochastic integration
        logx (bool, optional): If True, the x-axis will be in log scale.
            Defaults to False.
        logy (bool, optional): if True, the y-axis will be in log scale.
            Defaults to False.
        del_coef (float, optional): factor for modifying time steps used
            in the integration/propagation of ODE. Defaults to 10.
        normalize (bool, optional): True  will  be normalized the y axis
            based on max value . Defaults to False.
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
        mult_proc (bool, optional): If True, trajectories will be propa-
            gated on parallel. Defaults to False.
        items (tuple, optional): (canvas, scroll_x, scroll_y).
            Defaults to None.
        vary2 (str, optional): [description]. Defaults to "".
        implicit (bool, optional): True means report in time intervals
            similar to the input time intervals even if actual step is
            more or less. Defaults to False.
        rfile (str): file name of BioSANS topology file.
        exp_data_file ([type], optional): Experimental data file contai-
            ning True or accepted trajectories. Defaults to None.
    """

    stch_temp = []
    nzvar = []
    si_new = []
    slabels = list(sp_comp.keys())  # [xvar for xvar in sp_comp]
    for row in range(stoch_var.shape[0]):
        if np.sum(np.abs(stoch_var[row, :])) != 0 and slabels[row][0] != "-":
            stch_temp.append(list(stoch_var[row, :]))
            nzvar.append(row)
    stoch_var = np.array(stch_temp)
    sp_comp = {slabels[zvar]: sp_comp[slabels[zvar]] for zvar in nzvar}

    slabels = list(sp_comp.keys())
    data = []
    # if len(vary) > 0:
    if vary:
        hold = vary.pop()
        miter = len(vary)
    # if len(vary2) > 0:
    if vary2:
        r_index1, r_index2 = vary2[0]
        kval = vary2[1]
        miter = len(kval)
    if mult_proc:
        pool = multiprocessing.Pool(
            min(miter, round(globals2.CPU_MULT * multiprocessing.cpu_count())))
        concn_list = []
        ksns_list = []
        # if len(vary) > 0:
        if vary:
            for j in range(miter):
                concn[hold] = vary[j]
                concn_list.append({xvar: concn[xvar] for xvar in concn})
        else:
            concn_list = [concn for xvar in range(miter)]
        # if len(vary2) > 0:
        if vary2:
            for j in range(miter):
                ksns_list.append({xvar: ksn_dict[xvar][:]
                                  for xvar in ksn_dict})
                if isinstance(r_index1, int):
                    ksns_list[-1][r_index1][r_index2] = kval[j]
                else:
                    for k, _ in enumerate(r_index1):
                        ksns_list[-1][r_index1[k]][r_index2[k]] = kval[k][j]
        else:
            ksns_list = [ksn_dict for xvar in range(miter)]

            # if __name__ == '__main__':
        # if True:  # always true, just use above  command on some OS
        rands = [xvar * np.random.rand() for xvar in range(miter)]
        if method == "Tau-leaping":
            results = [
                pool.apply_async(
                    tau_leaping,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, rands[ih], del_coef,
                        implicit, rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        elif method == "Tau-leaping2":
            results = [
                pool.apply_async(
                    tau_leaping2,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, rands[ih], del_coef,
                        implicit, rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        elif method == "Sim-TauLeap":
            results = [
                pool.apply_async(
                    sim_tauLeap,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, rands[ih], del_coef,
                        implicit, rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        elif method == "CLE":
            results = [
                pool.apply_async(
                    cle_calculate,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, del_coef, rands[ih],
                        implicit, rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        elif method == "CLE2":
            results = [
                pool.apply_async(
                    cle2_calculate,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, del_coef, rands[ih],
                        rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        elif method == "Euler-1":
            results = [
                pool.apply_async(
                    euler_int,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, del_coef, False,
                        None, implicit, False, rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        elif method in ["Euler-2", "Euler-3"]:
            results = [
                pool.apply_async(
                    euler_int,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, del_coef, False,
                        None, implicit, True, rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        elif method == "LNA":
            save = False
            plot_show = False
            results = [
                pool.apply_async(
                    euler_int,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, del_coef, True,
                        items, False, rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + [slabels] for result in results]
        elif method == "Gillespie_":
            results = [
                pool.apply_async(
                    gillespie_ssa,
                    args=(
                        time_var, sp_comp, ksns_list[ih], concn_list[ih],
                        r_dict, p_dict, stoch_var, rands[ih], implicit,
                        rfile)
                ) for ih in range(miter)
            ]
            data = [result.get() + (slabels,) for result in results]
        else:
            print("Multiprocessing not supported for \
                your method of choice")
        pool.close()
    else:
        rands = [xvar * np.random.rand() for xvar in range(miter)]
        for j in range(miter):
            tnew = []
            r_seed = rands[j]
            # if len(vary) > 0:
            if vary:
                concn[hold] = vary[j]
            # if len(vary2) > 0:
            if vary2:
                if isinstance(r_index1, int):
                    ksn_dict[r_index1][r_index2] = kval[j]
                else:
                    for k, _ in enumerate(r_index1):
                        ksn_dict[r_index1[k]][r_index2[k]] = kval[k][j]
            if method == "CLE":
                tnew, zvar = cle_calculate(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, r_seed, implicit, rfile)
            elif method == "CLE2":
                tnew, zvar = cle2_calculate(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, r_seed, rfile)
            elif method == "Gillespie_":
                tnew, zvar = gillespie_ssa(
                    time_var, sp_comp, ksn_dict, concn, r_dict,
                    p_dict, stoch_var, r_seed, implicit, rfile)
            elif method == "Tau-leaping":
                tnew, zvar = tau_leaping(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, r_seed, del_coef, implicit, rfile)
            elif method == "Tau-leaping2":
                tnew, zvar = tau_leaping2(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, r_seed, del_coef, implicit, rfile)
            elif method == "Sim-TauLeap":
                tnew, zvar = sim_tauLeap(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, r_seed, del_coef, implicit, rfile)
            elif method == "Euler-1":
                tnew, zvar = euler_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, False, None, implicit, False,
                    rfile)
            elif method in ["Euler-2", "Euler-3"]:
                tnew, zvar = euler_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, False, None, implicit, True,
                    rfile)
            elif method == "Euler2-1":
                tnew, zvar = euler2_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, False, None, implicit, False,
                    rfile)
            elif method in ["Euler2-2", "Euler2-3"]:
                tnew, zvar = euler2_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, False, None, implicit, True,
                    rfile)
            elif method == "LNA":
                plot_show = False
                save = False
                tnew, zvar = euler_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, lna_solve=True, items=items,
                    rfile=rfile)
                tnew, zvar = euler2_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, False, None, implicit, True,
                    rfile)
            elif method == "LNA(t)":
                zvar, si_new, tnew = lna_non_steady_state(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, molar=True, rfile=rfile, del_coef=del_coef)
            elif method == "LNA2(t)":
                zvar, si_new, tnew = lna_non_steady_state2(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, molar=True, rfile=rfile, del_coef=del_coef)
            elif method == "LNA-vs":
                plot_show = False
                save = False
                tnew, zvar = lna_symbolic2(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, molar=True, mode="Numeric")
            elif method == "LNA-ks":
                plot_show = False
                save = False
                tnew, zvar = lna_symbolic2(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, molar=True, mode="fofks")
            elif method == "LNA-xo":
                plot_show = False
                save = False
                tnew, zvar = lna_symbolic2(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, molar=True, mode="fofCo")
            elif method == "LNA2":
                plot_show = False
                save = False
                tnew, zvar = lna_symbolic2(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items)
            elif method == "LNA3":
                plot_show = False
                save = False
                tnew, zvar = lna_symbolic2(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, molar=True)
            elif method == "NetLoc1":
                plot_show = False
                save = False
                tnew, zvar = law_loc_symbolic(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, molar=True)
            elif method == "NetLoc2":
                plot_show = False
                save = False
                tnew, zvar = law_loc_symbolic(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, molar=True, numer=True)
            elif method == "ODE-1":
                zvar = ode_int(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, False, rfile)
                tnew = time_var
            elif method in ["ODE-2", "ODE-3"]:
                zvar = ode_int(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, True, rfile)
                tnew = time_var
            # elif method == "Itoint-1":
                # zvar = SDE_int(concn, time_var, sp_comp, ksn_dict,
                               # r_dict, p_dict, stoch_var)
                # tnew = time_var
            # elif method in ["Itoint-2", "Itoint-3"]:
                # zvar = SDE_int(concn, time_var, sp_comp, ksn_dict, r_dict,
                               # p_dict, stoch_var, False)
                # tnew = time_var
            # elif method == "Stratint-1":
                # zvar = SDE_int(
                    # concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    # stoch_var, True, False)
                # tnew = time_var
            # elif method in ["Stratint-2", "Stratint-3"]:
                # zvar = SDE_int(
                    # concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    # stoch_var, False, False)
                # tnew = time_var
            elif method == "rk4-1":
                tnew, zvar = rungek4_int(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, False, del_coef, rfile)
            elif method in ["rk4-2", "rk4-3"]:
                tnew, zvar = rungek4_int(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, True, del_coef, rfile)
            elif method == "rk4-1a":
                tnew, zvar = rungek4a_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, False, implicit, rfile)
            elif method in ["rk4-2a", "rk4-3a"]:
                tnew, zvar = rungek4a_int(
                    time_var, sp_comp, ksn_dict, concn, r_dict, p_dict,
                    stoch_var, del_coef, True, implicit, rfile)
            elif method == "k_est1":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=True,
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est2":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, true_data_fil=exp_data_file,
                    rfile=rfile)
            elif method == "k_est3":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=True, mode="DEvol",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est4":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=False, mode="DEvol",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est5":
                set_p = [logx, logy, normalize, items]
                zvar = param_ode_int(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, True, rfile, set_p)
                plot_show = False
                save = False
                tnew = time_var
            elif method == "k_est6":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=True, mode="NeldMead",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est7":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=False, mode="NeldMead",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est8":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=True, mode="Powell",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est9":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=False, mode="Powell",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est10":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=True, mode="L-BFGS-B",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "k_est11":
                plot_show = False
                save = False
                tnew, zvar = param_estimate(
                    concn, time_var, sp_comp, ksn_dict, r_dict, p_dict,
                    stoch_var, items=items, molar=False, mode="L-BFGS-B",
                    true_data_fil=exp_data_file, rfile=rfile)
            elif method == "Analyt":
                plot_show = False
                save = False
                zvar = analyt_soln(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, rfile=rfile)
            elif method == "Analyt-ftx":
                plot_show = False
                save = False
                zvar = analyt_soln(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, rfile=rfile, mode="ftxo")
            elif method == "SAnalyt":
                plot_show = False
                save = False
                zvar = analyt_soln(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, rfile=rfile, not_semi=False)
            elif method == "SAnalyt-ftk":
                plot_show = False
                save = False
                zvar = analyt_soln(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, rfile=rfile, not_semi=False, mode="ftks")
            elif method == "Analyt2":
                plot_show = False
                save = False
                tnew, zvar = for_wxmaxima(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    items=items, rfile=rfile)
            elif method == "topoTosbml":
                plot_show = False
                save = False
                tnew, zvar = topo_to_sbml(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    v_volms, items=items, molar=False, rfile=rfile)
            elif method == "topoTosbml2":
                plot_show = False
                save = False
                tnew, zvar = topo_to_sbml(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    v_volms, items=items, molar=True, rfile=rfile)
            elif method == "topoTosbml3":
                plot_show = False
                save = False
                tnew, zvar = topo_to_sbml(
                    sp_comp, ksn_dict, concn, r_dict, p_dict, stoch_var,
                    v_volms, items=items, molar=None, rfile=rfile)
            data.append([tnew, zvar, slabels])

    nzvar = []
    reserve_events_words = {
        "t", "time", "status", "status2", "timer", "finish", "delay", "dtime"}
    for row in range(stoch_var.shape[0]):
        key = slabels[row].strip().split("_")[0]
        if key not in reserve_events_words:
            nzvar.append(row)
    sp_comp = {slabels[zvar]: sp_comp[slabels[zvar]] for zvar in nzvar}
    slabels = list(sp_comp.keys())

    # if len(si_new) > 0:
    if si_new:
        slabels = si_new
        sp_comp = si_new
        nzvar = range(len(slabels))

    if save:
        fname = out_fname + "_" + method + ".dat"
        file = open(fname, "w")
        file.write("time\t" + "\t".join(slabels) + "\n")
        for traj in data:
            xvar = traj[0]
            yvar = traj[1]
            xlen = len(xvar)
            for ihv in range(xlen):
                file.write(
                    str(xvar[ihv]) + "\t" + "\t"
                    .join([str(yi) for yi in yvar[ihv][nzvar]]) + "\n")
        file.close()

    if plot_show:

        yaxis_name = "conc"
        if method in ["rk4-1", "rk4-1a", "Euler-1", "Euler2-1", "ODE-1"]:
            yaxis_name = "molecules"
        elif method in ["rk4-2", "rk4-2a", "Euler-2", "Euler2-2", "ODE-2"]:
            yaxis_name = "moles/L"
        elif method in ["rk4-3", "rk4-3a", "Euler-3", "Euler2-3", "ODE-3"]:
            yaxis_name = "moles"
        elif method in ["CLE", "CLE2", "Tau-leaping",
                        "Tau-leaping2", "Sim-TauLeap"]:
            yaxis_name = "molecules"

        sp_len = len(sp_comp)
        if sp_len <= 10:
            col = ['C' + str(i) for i in range(10)]
        elif sp_len > 10 and sp_len <= 40:
            col = list(get_cmap("tab20").colors) \
                + list(get_cmap("tab20b").colors)
        else:
            col = list(mcd.CSS4_COLORS)
        if mix_plot:
            plt.figure(figsize=(9.68, 3.95))
            plt.xlabel(time_unit)
            plt.ylabel(yaxis_name)
            # if len(si_new) > 0:
            if si_new:
                plt.ylabel("cov" if si_new[0][0:2] != "FF" else "FF")
            lines = []
            if logx:
                plt.xscale('log')
            if logy:
                plt.yscale('log')

            for j in range(miter):
                for i in nzvar:
                    if normalize:
                        line = plt.plot(
                            data[j][0],
                            data[j][1][:, i] / (np.max(data[j][1][:, i])
                                                + 1.0e-30),
                            color=col[i])
                    else:
                        line = plt.plot(data[j][0], data[j]
                                        [1][:, i], color=col[i])
            if len(slabels) <= 10:
                plt.legend(slabels)
            else:
                plt.legend(slabels, fontsize='xx-small')
            plt.tight_layout()
            plt.savefig(out_fname + "_" + method + ".jpg")
            fig = plt.gcf()
            lines.append(line)
            globals2.PLOTTED.append([plt.gca(), fig, lines])
            # fig_canvas_agg = draw_figure(items, fig)
            draw_figure(items, fig)
            # plt.close()
        else:
            lines = []
            s_si = []
            for i in nzvar:
                s_si.append(slabels[i])
                plt.figure(figsize=(9.68, 3.95))
                plt.xlabel(time_unit)
                plt.ylabel(yaxis_name)
                # if len(si_new) > 0:
                if si_new:
                    plt.ylabel("cov" if si_new[0][0:2] != "FF" else "FF")
                if logx:
                    plt.xscale('log')
                if logy:
                    plt.yscale('log')
                for j in range(miter):
                    if normalize:
                        line = plt.plot(
                            data[j][0],
                            data[j][1][:, i] / (np.max(data[j][1][:, i])
                                                + 1.0e-30),
                            color=col[i])
                    else:
                        line = plt.plot(data[j][0], data[j]
                                        [1][:, i], color=col[i])
                    lines.append(line)
                plt.legend([slabels[i]])
                plt.tight_layout()
                plt.savefig(out_fname + "_" + method + "_" + str(i) + ".jpg")
                fig = plt.gcf()
                globals2.PLOTTED.append([plt.gca(), fig, lines])
                # fig_canvas_agg = draw_figure(items, fig)
                draw_figure(items, fig)
                # plt.close()
    return data
