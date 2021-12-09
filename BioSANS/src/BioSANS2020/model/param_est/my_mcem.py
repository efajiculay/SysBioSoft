"""

                 This module is the my_mcem module

This module is a simple implementation of monte-carlo expectation maximi
zation. The likelihood  function here is the joint probability distribu-
tion of the  error or  the difference between the true value and estima-
ted value. The error on each parameter estimate  is assumed to  be inde-
pendent and follows a normal  distribution. A uniform prior was set. The
log-likelihood  of  this function becomes the negative of the sum of the
squared  error  between  the  estimated value and true value. Parameters
were drawn  from  a log-normal distribution. In our  implementation, the
posterior  probability  is  calculated at each sampling  stage.  In  the
maximization step, the  ratio of posterior  probability between consecu-
tive draws (or the exponential of the difference between consecutive log
-likelihood) are  used to decide  whether to accept or reject the latest
values of the parameters based on a uniform random variable. After seve-
ral samplings  (decided programmatically in the algorithm), the mean and
standard deviation of the parameters  from all  the accepted  values are
calculated  and  used as the mean  and standard  deviation  for the next
sampling  stage.  This serves as the expectation step in our implementa-
tion because we  do not have a  close form for the parameters. The cycle
is repeated until  the calculated  mean and  standard deviations of each
parameters are no longer changing or  until  the maximum number of steps
is reached.

The following are the functions inside this module

1. log_likelihood
2. cost_value
3. exptn_maxtn
4. run_mcem


"""

# Metropolis Hasting Algorithm + Expectation Maximization Algorithm
# Personally coded by Erickson Fajiculay

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

# import time
# import math
# from scipy import linalg
# from scipy.optimize import fsolve
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint
# import random as rn

import warnings
import numpy as np
from BioSANS2020.myglobal import proc_global


warnings.filterwarnings('ignore')


def log_likelihood(ks_var, custom_function, args=None):
    """This function evaluates the loglikelihood based on the definitons
    provided in the custom_function. Here, it is the negative of the sum
    of squared errors between true and estimated value.

    Args:
        ks_var (list): list of parameter values
        custom_function (function): user/program defined objective/error
            function  to  optimized  my expectation  maximization  or EM
        args (tuple): Tuple  of  (data,  conc, tvar,  sp_comp,  ks_dict,
            r_dict,  p_dict,  v_stoich,  c_miss,  k_miss, molar, rfile).
            Defaults  to None. See param_estimate module for  the proper
            definition of variables.

    Returns:
        float: numeric value of custom_function(ks_var, args)
    """
    return custom_function(ks_var, args)


# def uniform_prior():
#     return 1.0


def cost_value(ks_var, custom_function, args=None):
    """This function  evaluates  the cost value or SSE of the definitons
    provided in the custom_function.

    Args:
        ks_var (list): list of parameter values
        custom_function (function): user/program defined objective/error
            function  to  optimized  my  expectation  maximization or EM
        args (tuple): Tuple  of  (data,  conc, tvar,  sp_comp,  ks_dict,
            r_dict,  p_dict,  v_stoich,  c_miss,  k_miss, molar, rfile).
            Defaults  to None. See param_estimate module for  the proper
            definition of variables.

    Returns:
        float: numeric value of custom_function(ks_var, args)
    """
    return abs(log_likelihood(ks_var, custom_function, args))


def exptn_maxtn(lst, seed_var, maxiter=50, inner_loop=1000, lenks=3,
                positive_only=False, likelihood=None, args=None, thr=1.0e-10):
    """This  function  performs  the  expectation  step and maximization
    steps  bounded by  maxiter and inner_loop  variables. The process is
    similar to the metrololis  hasting algorithm. There is a random walk
    performed as the  posterior probability is maximized. The parameters
    are updated in the randomwalk based on sampled values.

    Args:
        lst (multiprocessing.managers.ListProxy):
            multiprocessing.Manager() from proc_global module. This help
            in halting simulation when one of the chains already achived
            the defined tolerance.
        seed_var (float): random  seed  value picked  at random for each
            trajectory. They have been sampled from the calling program.
        maxiter (int, optional):  maximum  number  of expectation steps.
            Defaults to 50. This is the outer loop.
        inner_loop (int, optional): maximum number of maximization step.
            Defaults to 1000.
        lenks (int, optional): numbers  unknown  parameters to estimate.
        positive_only (bool, optional): constraint the parameter estima-
            tion to choose positive values. Defaults to False.
        likelihood (function, optional): user/program defined objective/
            error function  to  optimized
        args (tuple): Tuple  of  (data,  conc, tvar,  sp_comp,  ks_dict,
            r_dict,  p_dict,  v_stoich,  c_miss,  k_miss, molar, rfile).
            Defaults  to None. See param_estimate module for  the proper
            definition of variables.
        thr (function, optional): error tolerance. Defaults to 1.0e-10.

    Returns:
        tuple: (list of estimated parameter, minimum error)
    """

    np.random.seed(int(seed_var * 100))
    sd_var = np.random.uniform(0, 1)

    stds = [0] * lenks
    for i in range(lenks):
        if positive_only:
            stds[i] = abs(np.random.normal(sd_var, 0.01))
        else:
            stds[i] = abs(np.random.normal(sd_var, 100))

    kso = [0] * lenks
    for i in range(lenks):
        if positive_only:
            kso[i] = np.random.lognormal(0, stds[i])
        else:
            kso[i] = np.random.normal(0, stds[i])

    valo = [0] * lenks
    valp = [0] * lenks
    for i in range(lenks):
        valo[i] = 1
        valp[i] = 2

    ks_var = []
    for i in range(lenks):
        ks_var.append([])

    rs_var = []

    vmax = 1.0e+100
    iteration = 0
    test = 1000
    error = 1000
    # sigma = 1

    delta = int(np.ceil(inner_loop / (maxiter - 0.33 * maxiter)))
    inner_loop_orig = inner_loop
    inner_loop = delta
    while iteration < maxiter:
        if lst:
            return (valp, 1e+100)
        kis = []
        for i in range(lenks):
            kis.append([])

        r_big = vmax * log_likelihood(kso, likelihood, args)

        # Maximization step by metropolis hasting algorithm
        for i in range(inner_loop):
            if lst:
                break
            ksp = [0] * lenks
            for j in range(lenks):
                if positive_only:
                    ksp[j] = np.random.lognormal(
                        np.log(kso[j]), max(0, stds[j]))
                else:
                    ksp[j] = np.random.normal(kso[j], stds[j])

            l_big = vmax * log_likelihood(ksp, likelihood, args)

            u_rand = np.random.uniform(0, 1)
            check = np.exp(l_big - r_big)
            if check > u_rand:

                for j in range(lenks):
                    kis[j].append(ksp[j])
                    kso[j] = ksp[j]

                rs_var.append(l_big)
                r_big = l_big
            else:
                for j in range(lenks):
                    kis[j].append(kso[j])

        # Expectatation step by averaging
        for i in range(lenks):
            if len(kis[i]) > 1:
                for k in kis[i]:
                    ks_var[i].append(k)
                if positive_only:
                    stds[i] = np.std(np.log(kis[i]))
                else:
                    stds[i] = np.std(kis[i])
                valo[i] = valp[i]
                valp[i] = np.array(kis[i]).mean()
            else:
                stds[i] = abs(np.random.normal(2 * stds[i], 0.1 * stds[i]))

        test = 0
        for i in range(lenks):
            test = test + abs((valo[i] - valp[i]) / valo[i])

        if test < 1.0e-15:
            break

        iteration = iteration + 1
        inner_loop = min(inner_loop_orig, inner_loop + delta)

    error = min(error, cost_value(valp, likelihood, args))
    if error < thr:
        lst.append(1)  # early stopping of all chains
    return (valp, error)


def run_mcem(chains, n_pars, maxiter=5, inner_loop=5 * 1000,
             positive_only=False, likelihood=None, arg=None,
             r_rand=np.random.uniform(0, 1)):
    """[summary]

    Args:
        chains (integer): number of parallel runs or chains
        n_pars (integer): number of parameters
        maxiter (int, optional):  maximum  number  of expectation steps.
            Defaults to 50. This is the outer loop.
        inner_loop (int, optional): maximum number of maximization step.
            Defaults to 1000.
        positive_only (bool, optional): constraint the parameter estima-
            tion to choose positive values. Defaults to False.
        likelihood (function, optional): user/program defined objective/
            error function  to  optimized
        args (tuple): Tuple  of  (data,  conc, tvar,  sp_comp,  ks_dict,
            r_dict,  p_dict,  v_stoich,  c_miss,  k_miss, molar, rfile).
            Defaults  to None. See param_estimate module for  the proper
            definition of variables.
        r_rand (float, optional): np.random.uniform(0, 1) for seeding

    Returns:
        [type]: [description]
    """
    m_proc = proc_global.mp
    pool = m_proc.Pool(chains)
    np.random.seed(int(r_rand * 100))
    rands = [(xvar + 1) * np.random.uniform(0, 1) for xvar in range(chains)]

    thr = 1.0e-10
    results = [
        pool.apply_async(
            exptn_maxtn, args=(proc_global.lst, rands[ih], maxiter * (ih + 1),
                               inner_loop * (ih + 1), n_pars, positive_only,
                               likelihood, arg, thr)
        ) for ih in range(chains)
    ]

    ffvar = [result.get() for result in results]
    pool.close()

    ffvar = [result.get() for result in results]
    er_list = []
    for xvar in ffvar:
        er_list.append(xvar[1])

    er_min = min(er_list)
    ks_var = []
    for xvar in ffvar:
        if xvar[1] <= er_min:
            ks_var = xvar[0]
    return (ks_var, er_min)
