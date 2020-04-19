#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:04:15 2018
function for drawing parameters
@author: stsmall
"""
from __future__ import print_function
from __future__ import division
import numpy as np
from collections import defaultdict
from functools import partial


def unif(low, high):
    """
    """
    return(np.random.uniform(low, high))


def unifint(low, high):
    """
    """
    return(np.random.randint(low, high+1))


def normint(mu, sigma):
    """
    """
    return(np.round(sigma*np.random.randn()+mu))


def lognormint(mu, sigma):
    """
    """
    pass
    return(None)


def beta(a, b):
    """
    """
#                    import matplotlib.pyplot as plt
#                    mu = .5
#                    var = .05
#                    assert 0 < mu < 1  # .5
#                    assert 0 < var < .25  # .1
#                    alpha = ((((1-mu)/var) - (1/mu))*mu)**2
#                    beta = alpha*((1/mu) - 1)
#                    s = np.random.beta(alpha,beta,size=100000)
#                    count, bins, ignored = plt.hist(s, 100, normed=True, align='mid')
#                    mode = (alpha - 1) / (alpha + beta - 2)
#                    mu = alpha / (alpha + beta)
#                    var = alpha*beta / ((alpha + beta)^2*(alpha+beta+1))
    return(np.random.beta(a, b))


def constant(c):
    """
    """
    return(c)


def drawParams(paramsFile):
    """Draw parameters according to distribution specified in paramsfile
    """
    parfx = []
    parlist = []
    demodict = defaultdict(list)
    with open(paramsFile, 'r') as par:
        for line in par:
            if line.startswith("#"):
                pass
            elif line.startswith("tbi"):
                x = line.split()
                if "U" in x[3]:
                    if "int" in x[3]:
                        # divergence, Ne
                        parlist.append("{}{}".format(x[1], x[2]))
                        low = int(x[4])
                        high = int(x[5])
                        parfx.append(partial(unifint, low, high+1))
                    elif "flt" in x[3]:
                        # inheritance
                        parlist.append("{}{}".format(x[1], x[2]))
                        low = float(x[4])
                        high = float(x[5])
                        parfx.append(partial(unif, low, high))
                elif "N" in x[3]:
                    if "int" in x[3]:
                        # more confidence in prior for divergence, Ne
                        parlist.append("{}{}".format(x[1], x[2]))
                        mu = int(x[4])
                        sigma = int(x[5])
                        parfx.append(partial(normint, mu, sigma))
                    elif "log" in x[1]:
                        # draws values from normal, then log tansforms (no 0s)
                        parlist.append("{}{}".format(x[1], x[2]))
                        mu = int(x[4])
                        sigma = int(x[5])
                        parfx.append(partial(lognormint(mu, sigma)))
                elif "Beta" in x[3]:
                    # more confidence on inheritance
                    parlist.append("{}{}".format(x[1], x[2]))
                    a = float(x[4])
                    b = float(x[5])
                    parfx.append(partial(beta, a, b))
                elif "C" in x[3]:
                    parlist.append("{}{}".format(x[1], x[2]))
                    if "." in x[4]:
                        c = float(x[4])
                    else:
                        c = int(x[4])
                    parfx.append(partial(constant, c))
                else:
                    print("not a distribution")
                    raise ValueError
            else:
                x = line.split()
                demodict[int(x[0])].append(["{}{}".format(x[1], x[2]), x[3], x[4]])
    return(parfx, parlist, demodict)
