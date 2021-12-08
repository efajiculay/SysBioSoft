"""

                 This module is the my_mcem module

This module is a simple implementation of monte-carlo expectation maximi
zation. The likelihood function we used is the joint probability
distribution of the error between true value and estimated value, which
we assume to be independent and follow a normal distribution. A uniform
prior was assumed. The log-likelihood of this function becomes the
negative of the sum of the squared error between the estimated value and
true value. Parameters were drawn from a log-normal distribution. In our
implementation, the posterior probability is calculated at each sampling
stage. In the maximization step, the ratio of posterior probability
between consecutive draws (or the exponential of the difference between
consecutive log-likelihood) are used to decide whether to accept or
reject the latest values of the parameters based on a uniform random
variable. After several samplings (decided programmatically in the
algorithm), the mean and standard deviation of the parameters from all
the accepted values are calculated and used as the mean and standard
deviation for the next sampling stage. This serves as the expectation
step in our implementation because we do not have a close form for the
parameters. The cycle is repeated until the calculated mean and standard
deviations of each parameters are no longer changing or until the
maximum number of steps is reached.

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
    return custom_function(ks_var, args)


# def uniform_prior():
#     return 1.0


def cost_value(ks_var, custom_function, args=None):
    return abs(log_likelihood(ks_var, custom_function, args))


def mtrop_hstng(lst, sd_var, maxiter=50, inner_loop=1000, lenks=3,
                positive_only=False, likelihood=None, args=None, thr=1.0e-10):
    np.random.seed(int(sd_var * 100))
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
        lst.append(1)
    return (valp, error)


def run_mcem(chains, params, maxiter=5, inner_loop=5 * 1000,
             positive_only=False, likelihood=None, arg=None,
             r_rand=np.random.uniform(0, 1)):
    m_proc = proc_global.mp
    pool = m_proc.Pool(chains)
    np.random.seed(int(r_rand * 100))
    rands = [(xvar + 1) * np.random.uniform(0, 1) for xvar in range(chains)]

    thr = 1.0e-10
    results = [
        pool.apply_async(
            mtrop_hstng, args=(proc_global.lst, rands[ih], maxiter * (ih + 1),
                               inner_loop * (ih + 1), params, positive_only,
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
