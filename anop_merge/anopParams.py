#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:31:12 2020
@author: Scott T. Small

module for generating draws from priors for filet_sims.py
"""

import numpy as np
from collections import defaultdict
from functools import partial
from math import log
import re


def unif(low, high):
    """draws a random value from a uniform distribution, float

    Parameters
    ----------
    low: float
        lower bound
    high: float
        upper bound

    """
    return(np.random.uniform(low, high))


def unifint_L(low, high):
    """draws a random value from a uniform distribution, float

    Parameters
    ----------
    low: float
        lower bound
    high: float
        upper bound

    """
    log_low = log(low, 10)
    log_high = log(high, 10)
    rn = np.random.uniform(log_low, log_high)
    return(int(10**(rn)))


def unifint(low, high):
    """draws a random value from a uniform distribution, int

    Parameters
    ----------
    low: int
        lower bound
    high: int
        upper bound

    """
    return(np.random.randint(low, high+1))


def normint(mu, sigma):
    """draws a random value from a normal distribution, int

    Parameters
    ----------
    mu: int
        mean
    sigma: int
        variance

    """
    return(np.round(sigma*np.random.randn()+mu))


def lognormint(mu, sigma):
    """draws a random value from a log normal distribution

    Parameters
    ----------
    mu: float
        log mean
    sigma: float
        variance

    """
    # TODO: log normal
    pass
    return(None)


def beta_plot(alpha, beta):
    """draws a random value from a log normal distribution

    Parameters
    ----------
    mu: float
        log mean
    sigma: float
        variance

    """
    import matplotlib.pyplot as plt
    points = np.random.beta(alpha, beta, size=100000)
    count, bins, ignored = plt.hist(points, 100, normed=True, align='mid')


def beta(alpha, beta, plot=False):
    """draws a random value from a log normal distribution

    Parameters
    ----------
    a: float
        beta
    b: float
        beta

    """
#    alpha = ((((1-mu)/var) - (1/mu))*mu)**2
#    beta = alpha*((1/mu) - 1)
    mu = alpha / (alpha + beta)
    var = alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))
    mode = (alpha - 1) / (alpha + beta - 2)
    assert 0 < mu < 1
    assert 0 < var < .25
    if plot:
        beta_plot(alpha, beta)
        print(f"mode:{mode}, mu:{mu}, var: {var}")
    return(np.random.beta(alpha, beta))


def constant(c):
    """The number is a constant, just returns the number

    Parameters
    ----------
    c: float, int

    """
    return(c)


def get_dist(rDist):
    """
    """
    parPart = []
    for rD in rDist:
        y = rD.split()
        if "U" in y[0]:
            if "int" in y[0]:
                # divergence, Ne
                if "L" in y[0]:
                    low = int(y[1])
                    high = int(y[2])
                    partial(unifint_L, low, high)
                else:
                    low = int(y[1])
                    high = int(y[2])
                    partial(unifint, low, high)
            elif "flt" in y[0]:
                # p, 1-p is admixed proportion
                low = float(y[1])
                high = float(y[2])
                partial(unif, low, high)
    return(parPart)


def drawParams(params_file):
    """Creates a generator for each parameter based on the distribution

    Parameters
    ----------
    params_file: str
        file

    Returns
    -------
    par_gen: gen
        list of distributions
    par_list: list
        par
    demo_dict: default dict
        dict

    """
    par_gen = []
    par_list = []
    demo_dict = defaultdict(list)
    prior_cond = []
    pattern = re.compile(r'(r[aA-zZ]+) (tbi\d|.?\d*) (tbi\d|.?\d*)')
    with open(params_file, 'r') as par:
        for line in par:
            if line.startswith("#"):
                # indicates a comment line
                pass
            elif line.startswith("set"):
                prior_cond.append(line.split()[1:])
            elif line.startswith("tbi"):
                parms = line.split()
                tbd = parms[0]
                event = parms[1]
                pops = parms[2]
                par_list.append(f"{event}{pops}")
                rDist = re.findall(pattern, line)
                parPart = []
                for y in rDist:
                    breakpoint()
                    if "U" in y[0]:
                        if "int" in y[0]:
                            if "L" in y[0]:
                                low = int(y[1])
                                high = int(y[2])
                                parPart.append(partial(unifint_L, low, high+1))
                            else:
                                low = int(y[1])
                                high = int(y[2])
                                parPart.append(partial(unifint, low, high+1))
                        elif "flt" in y[0]:
                            # p, 1-p is admixed proportion
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
                    par_gen.append(parPart)
                else:
                    par_gen.append(parPart[0])
            else:
                # demography lists
                parms = line.split()
                time = parms[0]
                event = parms[1]
                pop = parms[2]
                Ne = parms[3]
                growth = parms[4]
                demo_dict[int(time)].append([f"{event}_{pop}_{Ne}_{growth}"])
    return(prior_cond, par_gen, par_list, demo_dict)
