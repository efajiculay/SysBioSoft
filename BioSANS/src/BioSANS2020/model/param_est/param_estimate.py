"""

                This is the param_est module

The following are the functions inside this module

1. load_data
2. custom_likelihood
3. ave_abs_dev
4. label_param
5. param_estimate


"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from tkinter import filedialog
import numpy as np
# from scipy.integrate import odeint
from scipy import optimize

from BioSANS2020.propagation.deterministic.ode_int import ODE_int
from BioSANS2020.model.param_est.my_mcem import run_mcem
from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.gui_functs.scrollable_text import INSERT
# from BioSANS2020.myglobal import proc_global


def load_data(file=None):
    """This function reads the trajectory file containing the data which
    should follow the following example format;

    time	A	B
    0.0	100.0	0.0
    0.25	88.24969025197632	11.750309748023732
    0.5	77.88007831231087	22.119921687689185
    0.75	68.72892784164061	31.27107215835944
    1.0	60.65306592491437	39.346934075085684
    1.25	53.526142785532	46.473857214468055
    1.5	47.236655135816875	52.76334486418318
    1.75	41.68620193454698	58.31379806545308
    2.0	36.78794415253036	63.21205584746969
    2.25	32.46524678349081	67.53475321650924
    ...

    The format  above have  a header where the first  column is time and
    all other columns are species or components. The rows are the values
    of measurements  corresponding to  the header. Each row is delimited
    by  tab character "\t". If  the data file is in excel, just copy the
    data  from  excel to a text editor, save  it with a  filename and it
    will already  be tab delimited. If there  are several  replicates of
    the  data, just append them to  the end without  the header and this
    function can still handle all replicates.

    time	A	B
    0.0	100.0	0.0                                 # first replicate
    ...                                             # continuation
    0.0	100.0	0.0                                 # second replicate
    ...                                             # continuation

    Args:
        file (string, optional): trajectory file. Defaults to None.

    Returns:
        tuple: tuple  of  data and  labels (data, labels). The data is a
            list of all  the trajectories in the file and the labels are
            the corresponding name of the columns in the trajectory
    """
    end_time = 0
    if file is None:
        file = filedialog.askopenfilename(title="Select file")
    with open(file, "r") as fvar:
        data = []
        ddvar = []
        row1 = str(fvar.readline()).strip()
        slabels = row1.split("\t")[1:]
        for row in fvar:
            cols = [float(xvar) for xvar in row.split("\t")]

            # if time decreases, it is a signal of new trajectory.
            if end_time > cols[0] and ddvar:
                data.append(np.array(ddvar))
                ddvar = []
            ddvar.append(cols)
            end_time = cols[0]
        data.append(np.array(ddvar))
        return (data, slabels)


def custom_likelihood(ks_par, args=None):
    """This  function returns  the negative of the sum of squared errors
    between the true value of trajectory and estimated trajectory. This
    serves  as the fitness/objective/cost function  with a maximum value
    of zero.

    Args:
        ks_par (list): list of parameters
        args (tuple): Tuple  of  (data,  conc, tvar,  sp_comp,  ks_dict,
            r_dict,  p_dict,  v_stoich,  c_miss,  k_miss, molar, rfile).
            Defaults  to None. See param_estimate function for  the pro-
            per definition of each variables.

    Returns:
        float: negative of the sum of squared error
    """
    data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict, \
        v_stoich, c_miss, k_miss, molar, rfile = args
    ind = 0
    slabels, data = data
    for row in k_miss:
        i, j = row
        ks_dict[i][j] = ks_par[ind]
        ind = ind + 1
    for spi in c_miss:
        conc[spi] = ks_par[ind]
        ind = ind + 1
    # print(conc)
    z_val = ODE_int(conc, tvar, sp_comp, ks_dict, r_dict,
                    p_dict, v_stoich, molar, rfile)
    spc1 = np.array([spi for spi in sp_comp])
    # print(z_val[0,np.isin(spc1,c_miss)])
    spc1 = ~np.isin(spc1, c_miss)
    spc2 = np.array(slabels)
    spc2 = ~np.isin(spc2, c_miss)
    return -np.sum((data[:, spc2] - z_val[:, spc1])**2)


def ave_abs_dev(ks_par, args=None):
    """This  function  returns  the  mean of relative absolute deviation
    between the true value of trajectory and estimated trajectory. This
    serves  as the fitness/objective/cost function  with a maximum value
    of zero.

    Args:
        ks_par (list): list of parameters
        args (tuple): Tuple  of  (data,  conc, tvar,  sp_comp,  ks_dict,
            r_dict,  p_dict,  v_stoich,  c_miss,  k_miss, molar, rfile).
            Defaults  to None. See param_estimate function for  the pro-
            per definition of each variables.

    Returns:
        float: mean of relative absolute deviation
    """
    data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict, \
        v_stoich, c_miss, k_miss, molar, rfile = args
    ind = 0
    slabels, data = data
    for row in k_miss:
        i, j = row
        ks_dict[i][j] = ks_par[ind]
        ind = ind + 1
    for spi in c_miss:
        conc[spi] = ks_par[ind]
        ind = ind + 1
    z_val = ODE_int(conc, tvar, sp_comp, ks_dict, r_dict,
                    p_dict, v_stoich, molar, rfile)
    spc1 = np.array([spi for spi in sp_comp])
    spc1 = ~np.isin(spc1, c_miss)
    spc2 = np.array(slabels)
    spc2 = ~np.isin(spc2, c_miss)
    denom = data[:, spc2]
    denom[np.abs(denom) < 1.0e-10] = 1
    return np.mean(np.abs((data[:, spc2] - z_val[:, spc1]) / denom))


def label_param(k_miss, c_miss, ks_par):
    """[summary]

    Args:
        k_miss (list): list of position of unkwon parameters in the
            reaction tag of BioSANS topology file. For example;

            #REACTIONS
            A => B , -1      # forward rate constant is unknown
            B <=> C, -1, -1  # forward and backward is unknown

            The negative values means rate constant is unknown for those
            reaction.

            The value of k_miss will be;

                k_miss = [(0,0), (1, 0), (1, 1)]

            The format of the tuple is;

                (reaction index, rate constant index)
                (0, 0) - first reaction, first rate constant
                (1, 0) - second reaction, first rate constant
                (1, 1) - second reaction, second rate constant

        c_miss (list): component/species without trajectory data. They
            can be set in a topology file with negative values as well.

            @CONCENTRATION
            A, Ao
            B, -1
            C, Co

            The value of c_miss will be;

                c_miss = ['B']

        ks_par (list): list of estimated parameter values for k_miss as
            it appears sequencially and list of initial value for c_miss
            as it appears sequencially.

            Example;

                ks_par = [
                    k1f, # numeric value for position (0, 0)
                    k2f, # numeric value for position (1, 0)
                    k2b, # numeric value for position (1, 1)
                    Bo   # numeric initial value for B
                ]

    Returns:
        dictionary: dictionary  of estimated  values  with the following
        format {label : estimate,...}
    """
    pars = {}
    count = 0
    lvar = ["kf", "kb"]
    for xvar in k_miss:
        i, j = xvar
        pars[lvar[j] + str(i + 1)] = ks_par[count]
        count = count + 1

    for xvar in c_miss:
        pars[xvar + "o"] = ks_par[count]
        count = count + 1
    return pars


def param_estimate(conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
                   v_stoich, items, molar=False, mode="MCEM",
                   true_data_fil=None, rfile=None):
    """[summary]

    Args:
        conc (dict): dictionary of initial concentration.

            For example;

                {'A': 100.0, 'B': -1.0, 'C': 0.0}
                negative means unknown or for estimation

        tvar (list): time stamp of trajectories i.e. [0, 0.1, 0.2, ...]
        sp_comp (dict): dictionary of appearance or position of species
            or component in the reaction tag of BioSANS topology file.

            For example;

                #REACTIONS
                A => B, -kf1
                B => C, kf2

            The value of sp_comp is

                sp_comp = {'A': {0}, 'B': {0, 1}, 'C': {1}}

                A appears in first reaction with index 0
                B appears in first and second reaction with index 0, 1
                C appears in second reaction with index 1

        ks_dict (dict): dictionary of rate constant that appears in each
            reactions. For the above reaction, the value of ks_dict is

            For example;

                #REACTIONS
                A => B , 0.3        # first reaction
                B <=> C, 0.1, 0.2   # second reaction

            The value of ks_dict is

                ks_dict = {
                    0 : [0.3],      # first reaction
                    1 : [0.1, 0.2]  # second reaction
                }

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
        v_stoich (numpy.ndarray): stoichiometric matrix. For example

            v_stoich = np.array([
                [   -1,           0   ]            # species A
                [    1,          -1   ]            # species B
                [    0,           1   ]            # species C
                  #1st rxn    2nd rxn
            ])

        items (list): list of [canvas, scroll_x, scroll_y]
        molar (bool, optional): If True, the units for any amount is in
            molar. Propensity will be macroscopic. Defaults to False.
        mode (str, optional): parameter estimation method keyword which
            is  one of "MCEM", "DEvol", "NeldMead", "Powell", "L-BFGS-B"
        true_data_fil (string, optional): file  name  of  the trajectory
            file containing experimental or observed trajectory.
            Defaults to None.
        rfile (string, optional): name of topology file where some
            parameters or components are negative indicating  they  have
            to be estimated. Defaults to None.

    Returns:
        tuple: (0, param_res) where param_res is a dictionary of estima-
        ted  values  with the following format {label : estimate,...}
    """

    data2 = load_data(true_data_fil)
    slabels = data2[1]
    spc = [spi for spi in sp_comp]
    z_val = [conc[a] for a in sp_comp]
    c_miss = []
    for i, _ in enumerate(z_val):
        if z_val[i] < 0:
            c_miss.append(spc[i])
    k_miss = []
    for i, _ in enumerate(ks_dict):
        if len(ks_dict[i]) == 1:
            if ks_dict[i][0] < 0:
                k_miss.append((i, 0))
        elif len(ks_dict[i]) == 2:
            if ks_dict[i][0] < 0:
                k_miss.append((i, 0))
            if ks_dict[i][1] < 0:
                k_miss.append((i, 1))
    npar = len(c_miss) + len(k_miss)

    tvar = data2[0][0][:, 0]
    data = (slabels, data2[0][0][:, 1:])
    # print(data)

    if items:
        text = prepare_scroll_text(items)

        def ffprint(xvar):
            return text.insert(INSERT, " ".join([str(y) for y in xvar]))
    else:
        def ffprint(xvar):
            return print(" ".join([str(y) for y in xvar]))

    param_res = {}
    if mode == "MCEM":
        fvar = 1
        k = 10
        ks_par, er_min = run_mcem(
            fvar, npar, fvar, k * fvar, positive_only=True,
            likelihood=custom_likelihood,
            arg=(data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
                 v_stoich, c_miss, k_miss, molar, rfile)
        )

        error = ave_abs_dev(
            ks_par, (data, conc, tvar, sp_comp, ks_dict, r_dict,
                     p_dict, v_stoich, c_miss, k_miss, molar, rfile)
        )

        count = 0
        best_er = 1e+100
        best_ks = ks_par

        rands = [(xvar + 1) * np.random.uniform(0, 1)
                 for xvar in range(200000)]
        ind = 0
        while error > 1.0e-6 and k <= 10000:
            ks_par, er_min = run_mcem(
                min(fvar, 10), npar, fvar, k * fvar, positive_only=True,
                likelihood=custom_likelihood,
                arg=(data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
                     v_stoich, c_miss, k_miss, molar, rfile), r_rand=rands[ind]
            )

            error = ave_abs_dev(
                ks_par, (data, conc, tvar, sp_comp, ks_dict, r_dict,
                         p_dict, v_stoich, c_miss, k_miss, molar, rfile))
            if error < best_er:
                best_er = error
                best_ks = ks_par

                param_res = label_param(k_miss, c_miss, best_ks)
                ffprint(["\n\npartial result =\n", "errors = ", error, er_min])
                ffprint(["\n"])
                for xvar in param_res:
                    ffprint(["\n" + str(xvar) + " = " + str(param_res[xvar])])

            fvar = fvar + 1
            count = count + 1
            if count == 10:
                k = k * 10
                fvar = 1
                count = 0
            ind = ind + 1
        param_res = label_param(k_miss, c_miss, best_ks)
        ffprint(["\n\nFinal result =\n", best_ks, error, er_min])
        ffprint(["\n"])
        for xvar in param_res:
            ffprint(["\n" + str(xvar) + " = " + str(param_res[xvar])])
    elif mode == "DEvol":
        x_ini = []
        for _ in k_miss:
            x_ini.append(np.random.uniform())
        for _ in c_miss:
            x_ini.append(np.random.uniform())
        ddvar = optimize.minimize(
            lambda ks_par: -
            custom_likelihood(
                ks_par, (data, conc, tvar, sp_comp, ks_dict, r_dict,
                         p_dict, v_stoich, c_miss, k_miss, molar, rfile)
            ),
            x_ini, method='Nelder-Mead', tol=1e-10,
            options={
                'maxiter': 10000,
                'adaptive': True})

        bounds = []
        counts = 0
        for _ in k_miss:
            val = ddvar.x[counts]
            bounds.append((max(0, 0.1 * val), max(0, 10 * val)))
            counts = counts + 1
        for _ in c_miss:
            val = ddvar.x[counts]
            bounds.append((max(0, 0.1 * val), max(0, 10 * val)))
            counts = counts + 1

        ddvar = optimize.differential_evolution(
            lambda ks_par: -custom_likelihood
            (ks_par, (data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
                      v_stoich, c_miss, k_miss, molar, rfile)),
            bounds, maxiter=100, popsize=10, tol=1.0e-10)

        ddvar = optimize.minimize(
            lambda ks_par: -custom_likelihood
            (ks_par, (data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
                      v_stoich, c_miss, k_miss, molar, rfile)),
            ddvar.x, method='Nelder-Mead', tol=1e-10,
            options={'maxiter': 10000, 'adaptive': True})

        param_res = label_param(k_miss, c_miss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for xvar in param_res:
            ffprint(["\n" + str(xvar) + " = " + str(param_res[xvar])])
    elif mode == "NeldMead":
        x_ini = []
        for _ in k_miss:
            x_ini.append(np.random.uniform())
        for _ in c_miss:
            x_ini.append(np.random.uniform())

        ddvar = optimize.minimize(
            lambda ks_par: -custom_likelihood(
                ks_par, (data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
                         v_stoich, c_miss, k_miss, molar, rfile)), x_ini,
            method='Nelder-Mead', tol=1e-10,
            options={'maxiter': 100000, 'adaptive': True})

        param_res = label_param(k_miss, c_miss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for xvar in param_res:
            ffprint(["\n" + str(xvar) + " = " + str(param_res[xvar])])
    elif mode == "Powell":
        x_ini = []
        for _ in k_miss:
            x_ini.append(np.random.uniform())
        for _ in c_miss:
            x_ini.append(np.random.uniform())

        ddvar = optimize.minimize(
            lambda ks_par: -custom_likelihood(
                ks_par,
                (data, conc, tvar, sp_comp, ks_dict, r_dict,
                 p_dict, v_stoich, c_miss, k_miss, molar, rfile)
            ),
            x_ini, method='Powell', tol=1e-10,
            options={'maxiter': 100000, 'ftol': 1e-10, 'xtol': 1e-10})

        param_res = label_param(k_miss, c_miss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for xvar in param_res:
            ffprint(["\n" + str(xvar) + " = " + str(param_res[xvar])])
    elif mode == "L-BFGS-B":
        x_ini = []
        for _ in k_miss:
            x_ini.append(np.random.uniform())
        for _ in c_miss:
            x_ini.append(np.random.uniform())

        ddvar = optimize.minimize(
            lambda ks_par: -custom_likelihood
            (ks_par, (data, conc, tvar, sp_comp, ks_dict, r_dict, p_dict,
                      v_stoich, c_miss, k_miss, molar, rfile)
             ), x_ini, method='L-BFGS-B', tol=1e-10,

            options={'maxiter': 100000, 'ftol': 1e-10, 'xtol': 1e-10})

        param_res = label_param(k_miss, c_miss, ddvar.x)
        ffprint(["\nFinal result =\n", ddvar])
        ffprint(["\n"])
        for xvar in param_res:
            ffprint(["\n" + str(xvar) + " = " + str(param_res[xvar])])
    return (0, param_res)
