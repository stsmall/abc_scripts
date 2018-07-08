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
import re


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
    pattern = re.compile(r'r[aA-zZ]* \d*\.?\d* \d*\.?\d*')
    with open(paramsFile, 'r') as par:
        for line in par:
            if line.startswith("#"):
                pass
            elif line.startswith("tbi"):
                x = line.split()
                parlist.append("{}{}".format(x[1], x[2]))
                rDist = re.findall(pattern, line)
                parPart = []
                for r in rDist:
                    y = r.split()
                    if "U" in y[0]:
                        if "int" in y[0]:
                            # divergence, Ne
                            low = int(y[1])
                            high = int(y[2])
                            parPart.append(partial(unifint, low, high+1))
                        elif "flt" in y[0]:
                            # inheritance
                            low = float(y[1])
                            high = float(y[2])
                            parPart.append(partial(unif, low, high))
                    elif "N" in y[0]:
                        if "int" in y[0]:
                            # more confidence in prior for divergence, Ne
                            mu = int(y[1])
                            sigma = int(y[2])
                            parPart.append(partial(normint, mu, sigma))
                        elif "log" in y[0]:
                            # draws values from normal then log tansforms (no 0s)
                            mu = int(y[1])
                            sigma = int(y[2])
                            parPart.append(partial(lognormint(mu, sigma)))
                    elif "B" in y[0]:
                        # more confidence on inheritance
                        a = float(y[1])
                        b = float(y[2])
                        parPart.append(partial(beta, a, b))
                    elif "C" in y[0]:
                        if "." in y[1]:
                            c = float(y[1])
                        else:
                            c = int(y[1])
                        parPart.append(partial(constant, c))
                    else:
                        print("not a recognized distribution")
                        raise ValueError
                if len(rDist) > 1:
                    parfx.append(parPart)
                else:
                    parfx.append(parPart[0])
            else:
                x = line.split()
                demodict[int(x[0])].append(["{}{}".format(x[1], x[2]), x[3], x[4]])
    return(parfx, parlist, demodict)
