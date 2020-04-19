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
    """Draws a random value from a uniform distribution, float.

    Parameters
    ----------
    low: float
        lower bound
    high: float
        upper bound

    """
    return(np.random.uniform(low, high))


def unifint_L(low, high):
    """Draws a random value from a uniform distribution, float.

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
    """Draws a random value from a uniform distribution, int.

    Parameters
    ----------
    low: int
        lower bound
    high: int
        upper bound

    """
    return(np.random.randint(low, high+1))


def normint(mu, sigma):
    """Draws a random value from a normal distribution, int.

    Parameters
    ----------
    mu: int
        mean
    sigma: int
        variance

    """
    return(np.round(sigma*np.random.randn()+mu))


def lognormint(mu, sigma):
    """Draws a random value from a log normal distribution.

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
    """Draws a random value from a log normal distribution.

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
    """Draws a random value from a log normal distribution.

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


def constant(low, high):
    """Return a constant.

    Parameters
    ----------
    c: float, int

    """
    if low == high:
        c = low
    return(c)


def get_dist(par_gen):
    """Create parameter list from function list."""
    # TODO: tbi can only be in first distribution
    params_list = []
    tmp_list = []
    i = -1
    while abs(i) <= len(par_gen):
        px = par_gen[i]
        if type(px) is list:
            dist, low, high = px[0]
            try:
                rr = dist(low, high)
            except TypeError:
                try:
                    rr = dist(low, tmp_list[-1])
                except TypeError:
                    l_ix = int(low[-1])+1
                    h_ix = int(high[-1])+1
                    rr = dist(tmp_list[-l_ix], tmp_list[-h_ix])
            tmp_list.append(rr)
            dist, low, high = px[1]
            rrf = dist(low, high)
            if len(px) > 2:
                dist, low, high = px[2]
                rrfc = dist(low, high)
                params_list.append([rr, rrf, rrfc])
            else:
                params_list.append([rr, rrf])
        else:
            dist, low, high = px
            try:
                rr = dist(low, high)
            except TypeError:
                try:
                    rr = dist(low, tmp_list[-1])
                except TypeError:
                    l_ix = int(low[-1]) + 1
                    h_ix = int(high[-1]) + 1
                    rr = dist(tmp_list[-l_ix], tmp_list[-h_ix])
            params_list.append(rr)
            tmp_list.append(rr)
        i -= 1
    return(params_list)


def drawParams(params_file):
    """Create a generator for each parameter based on the distribution.

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
    pattern = re.compile(r'(r[aA-zZ]+) (tbi\d|0?.?\d*) (tbi\d|0?.?\d*)')
    with open(params_file, 'r') as par:
        for line in par:
            if line.startswith("#"):
                # indicates a comment line
                pass
            elif line.startswith("tbi"):
                parms = line.split()
                tbd = parms[0]
                event = parms[1]
                pops = parms[2]
                par_list.append(f"{event}{pops}")
                rDist = re.findall(pattern, line)
                parPart = []
                for y in rDist:
                    dist = y[0]
                    low = y[1]
                    high = y[2]
                    if low.isdigit():
                        low = int(y[1])
                    elif "." in low:
                        low = float(y[1])
                    if high.isdigit():
                        high = int(y[2])
                    elif "." in high:
                        high = float(y[2])
                    if "U" in dist:
                        if "int" in dist:
                            if "L" in dist:
                                parPart.append((unifint_L, low, high))
                            else:
                                parPart.append((unifint, low, high))
                        elif "flt" in dist:
                            parPart.append((unif, low, high))
                    elif "N" in dist:
                        if "int" in dist:
                            # more confidence in prior for divergence, Ne
                            mu = low
                            sigma = high
                            parPart.append((normint, mu, sigma))
                        elif "log" in dist:
                            # draws values from normal then log tansforms (no 0s)
                            mu = low
                            sigma = high
                            parPart.append((lognormint, mu, sigma))
                    elif "B" in dist:
                        # more confidence on inheritance
                        a = low
                        b = high
                        parPart.append((beta, a, b))
                    elif "C" in dist:
                        if low == high:
                            parPart.append((constant, low, high))
                        else:
                            print("constant is not constant")
                            raise ValueError
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
    return(par_gen, par_list, demo_dict)
