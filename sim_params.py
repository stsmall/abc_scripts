#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:31:12 2020
@author: Scott T. Small

module for generating draws from priors for filet_sims.py
"""

import numpy as np
from collections import defaultdict
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


def get_dist(par_dict):
    """Create parameter list from function list."""
    event_dict = {}
    key_count = list(par_dict.keys()).copy()
    while len(key_count) > 0:
        keys = key_count
        for tbi in keys:
            tmp_list = []
            for event in par_dict[tbi]:
                dist, low, high = event
                if any(type(t) == str for t in event):
                    try:
                        if type(low) == str:
                            low = event_dict[low][0]
                        if type(high) == str:
                            high = event_dict[high][0]
                        tmp_list.append(dist(low, high))
                    except KeyError:
                        break
                else:
                    tmp_list.append(dist(low, high))
            if tmp_list:
                event_dict[tbi] = tmp_list
                key_count.remove(tbi)
    return(event_dict)


def drawParams(params_file):
    """Create a generator for each parameter based on the distribution.

    Parameters
    ----------
    params_file: str
        file

    Returns
    -------
    par_dict: Dict
        list of distributions
    event_list: List
        par
    demo_dict: default dict
        dict
    cond_list: List
        list of parameter conditions like lt, gt

    """
    demo_dict = defaultdict(list)
    par_dict = {}
    cond_list = []
    event_dict = defaultdict(list)
    pattern = re.compile(r'(r[aA-zZ]+) (tbi\d|0?.?\d*) (tbi\d|0?.?\d*)')
    with open(params_file, 'r') as par:
        for line in par:
            if line.strip():
                if line.startswith("#"):
                    pass
                elif line.startswith("set"):
                    cond_list.append(line.split()[1:])
                elif line.startswith("tbi"):
                    tbi, event, pops, *_ = line.split()
                    event_str = f"{event}_{pops}"
                    event_dict[event_str].append(tbi)
                    r_dist = re.findall(pattern, line)
                    parPart = []
                    for y in r_dist:
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
                                    dist = unifint_L
                                else:
                                    dist = unifint
                            elif "flt" in dist:
                                dist = unif
                            parPart.append((dist, low, high))
                        elif "N" in dist:
                            mu = low
                            sigma = high
                            if "int" in dist:
                                dist = normint
                            elif "log" in dist:
                                # draws values from normal then log tansforms (no 0s)
                                dist = lognormint
                            parPart.append((dist, mu, sigma))
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
                    par_dict[tbi] = parPart
                else:
                    # demography lists
                    time, event, pop, Ne, growth = line.split()
                    demo_dict[int(time)].append(f"{event}_{pop}_{Ne}_{growth}")
    return(cond_list, par_dict, event_dict, demo_dict)
