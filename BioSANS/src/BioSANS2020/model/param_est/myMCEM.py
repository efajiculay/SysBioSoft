#!/usr/bin/env python
# coding: utf-8

# Metropolis Hasting Algorithm + Expectation Maximization Algorithm

# Personally coded by Erickson Fajiculay

#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.myglobal import proc_global as proc_global
import time
import math
from scipy import linalg
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random as rn
import numpy as np
import warnings
warnings.filterwarnings('ignore')


def log_likelihood(ks, custom_function, args=None):
    return custom_function(ks, args)


def prior(x, std=0):
    return 1.0


def cost_value(ks, custom_function, args=None):
    return abs(log_likelihood(ks, custom_function, args))


def MH_call(lst, sd, maxiter=50, inner_loop=1000, lenks=3, positive_only=False, likelihood=None, args=None, thr=1.0e-10):
    np.random.seed(int(sd*100))
    sd = np.random.uniform(0, 1)

    stds = [0]*lenks
    for i in range(lenks):
        if positive_only:
            stds[i] = abs(np.random.normal(sd, 0.01))
        else:
            stds[i] = abs(np.random.normal(sd, 100))

    kso = [0]*lenks
    for i in range(lenks):
        if positive_only:
            kso[i] = np.random.lognormal(0, stds[i])
        else:
            kso[i] = np.random.normal(0, stds[i])

    valo = [0]*lenks
    valp = [0]*lenks
    for i in range(lenks):
        valo[i] = 1
        valp[i] = 2

    Ks = []
    for i in range(lenks):
        Ks.append([])

    Rs = []

    vmax = 1.0e+100
    iteration = 0
    test = 1000
    error = 1000
    sigma = 1

    delta = int(np.ceil(inner_loop/(maxiter-0.33*maxiter)))
    inner_loop_orig = inner_loop
    inner_loop = delta
    while iteration < maxiter:
        if lst:
            return (valp, 1e+100)
        kis = []
        for i in range(lenks):
            kis.append([])

        R = vmax*log_likelihood(kso, likelihood, args)

        for i in range(inner_loop):
            if lst:
                break
            ksp = [0]*lenks
            for j in range(lenks):
                if positive_only:
                    ksp[j] = np.random.lognormal(
                        np.log(kso[j]), max(0, stds[j]))
                else:
                    ksp[j] = np.random.normal(kso[j], stds[j])

            L = vmax*log_likelihood(ksp, likelihood, args)

            U = np.random.uniform(0, 1)
            check = np.exp(L-R)
            if check > U:

                for j in range(lenks):
                    kis[j].append(ksp[j])
                    kso[j] = ksp[j]

                Rs.append(L)
                R = L
            else:
                for j in range(lenks):
                    kis[j].append(kso[j])

        for i in range(lenks):
            if len(kis[i]) > 1:
                [Ks[i].append(k) for k in kis[i]]
                if positive_only:
                    stds[i] = np.std(np.log(kis[i]))
                else:
                    stds[i] = np.std(kis[i])
                valo[i] = valp[i]
                valp[i] = np.array(kis[i]).mean()
            else:
                stds[i] = abs(np.random.normal(2*stds[i], 0.1*stds[i]))

        test = 0
        for i in range(lenks):
            test = test + abs((valo[i]-valp[i])/valo[i])

        if test < 1.0e-15:
            break

        iteration = iteration + 1
        inner_loop = min(inner_loop_orig, inner_loop + delta)

    error = min(error, cost_value(valp, likelihood, args))
    if error < thr:
        lst.append(1)
    return (valp, error)


def run_MCEM(chains, params, maxiter=5, inner_loop=5*1000, positive_only=False, likelihood=None, arg=None, rr=np.random.uniform(0, 1)):
    mp = proc_global.mp
    pool = mp.Pool(chains)
    np.random.seed(int(rr*100))
    rands = [(x+1)*np.random.uniform(0, 1) for x in range(chains)]

    thr = 1.0e-10
    results = [pool.apply_async(MH_call, args=(proc_global.lst, rands[ih], maxiter*(
        ih+1), inner_loop*(ih+1), params, positive_only, likelihood, arg, thr)) for ih in range(chains)]

    ffvar = [result.get() for result in results]
    pool.close()

    ffvar = [result.get() for result in results]
    er = []
    for x in ffvar:
        er.append(x[1])

    er_min = min(er)
    ks = []
    for x in ffvar:
        if x[1] <= er_min:
            ks = x[0]
    return (ks, er_min)
