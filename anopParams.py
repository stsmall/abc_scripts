#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:04:15 2018
function for drawing parameters
@author: stsmall
"""

import numpy as np
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
    with open(paramsFile, 'r') as par:
        for line in par:
            if line.startswith("#"):
                pass
            else:
                x = line.split()
                if "uniform_int" in x[1]:
                    # divergence, Ne
                    parlist.append(x[0])
                    low = int(x[2])
                    high = int(x[3])
                    parfx.append(partial(unifint, low, high))
                    # parfx.append(lambda: np.random.randint(low, high+1))
                elif "uniform_flt" in x[1]:
                    # inheritance
                    parlist.append(x[0])
                    low = float(x[2])
                    high = float(x[3])
                    parfx.append(partial(unif, low, high))
                    # parfx.append(lambda: np.random.uniform(low, high))
                elif "norm_int" in x[1]:
                    # more confidence in prior for divergence, Ne
                    parlist.append(x[0])
                    mu = int(x[2])
                    sigma = int(x[3])
                    parfx.append(partial(normint, mu, sigma))
                    # parfx.append(lambda: np.round(sigma*np.random.randn()+mu))
                elif "log_norm_int" in x[1]:
                    # draws values from normal, then log tansforms (no 0s)
                    parlist.append(x[0])
                    mu = int(x[2])
                    sigma = int(x[3])
                    parfx.append(partial(lognormint(mu, sigma)))
                elif "beta" in x[1]:
                    # more confidence on inheritance
                    parlist.append(x[0])
                    a = float(x[2])
                    b = float(x[3])
                    parfx.append(partial(beta, a, b))
                    # parfx.append(lambda: np.random.beta(a, b))
                elif "constant" in x[1]:
                    parlist.append(x[0])
                    if "." in x[2]:
                        c = float(x[2])
                    else:
                        c = int(x[2])
                    parfx.append(partial(constant, c))
                else:
                    print("not a distribution")
                    raise ValueError
    return(parfx, parlist)
