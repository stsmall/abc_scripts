#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 17:21:53 2018

anopSims_main.py --cfg FILE.cfg --model 1

@author: stsmall
"""
import sys
import numpy as np
import msprime as msp
import allel
from collections import defaultdict
from anopParams import drawParams
from anopModels import Model
import anopStats
import argparse


def checkDemo(popconfig, migmat, demoevents):
    """Print out demography for debugging
    """
    #f = open("demographyHistory.txt", 'a')
    mig = list(migmat)
    dd = msp.DemographyDebugger(
        population_configurations=popconfig,
        migration_matrix=mig,
        demographic_events=demoevents)
    dd.print_history()
    #f.write("{}".format(dd.print_history()))
    #f.close()
    return(None)

def simulate(model, demodict, ix, parfx, parlist, thetaarray, rhoarray):
    """Runs msprime with model
    """
    treelist = []
    params = [p() for p in parfx]
    # avoid exactly equal times by small perturbation
    while len(params) != len(set(params)):
        params = [i + np.random.randint(-10, 10) for i in params]
    # simulate
    for loc in range(model["loci"]):
        Ne = int(np.round((np.random.choice(thetaarray)/(4*model["mutation_rate"]))))
        m = Model()
        dem_events = m.genDem(ix, params, demodict, parlist, Ne, model["mMax"], model["mIso"])
        recomb = np.random.choice(rhoarray)/(4*Ne)
        checkDemo(model["pop_config"], model["migMat"], dem_events)
        ts = msp.simulate(
                          population_configurations=model["pop_config"],
                          demographic_events=dem_events,
                          Ne=model["Ne"],
                          migration_matrix=model["migMat"],
                          length=model["contig_length"],
                          mutation_rate=model["mutation_rate"],
                          recombination_rate=recomb,
                          )
        treelist.append(ts)
    # stats_out = anopStats(treearray)
    # writeABC(stats_out, params, ix)
    return(None)


def setInitial(sampleSize, intitialSize, growthRate):
    """Set initial population sizes, growth rates, and samples

    Parameters
    ------
    sampleSize: list, of diploid sample sizes
    initialSize: list, of Ne
    growthRate: list, of growth rates

    Returns
    ------
    pop_config: msprime object

    """
    # inital effective sizes
    N_fun, N_like, N_van, N_lon, N_par = initialSize
#    N_fun, N_like, N_van, N_lon, N_par, N_riv = initialSize
    # sample size as diploid
    S_fun, S_like, S_van, S_lon, S_par = sampleSize
#    S_fun, S_like, S_van, S_lon, S_par, S_riv = sampleSize
    # intial growth rates
    G_fun, G_like, G_van, G_lon, G_par = growthRate
 #   G_fun, G_like, G_van, G_lon, G_par, G_riv = growthRate
    pop_config = [
                  msp.PopulationConfiguration(sample_size=S_fun,
                                              initial_size=N_fun,
                                              growth_rate=G_fun),
                  msp.PopulationConfiguration(sample_size=S_like,
                                              initial_size=N_like,
                                              growth_rate=G_like),
                  msp.PopulationConfiguration(sample_size=S_van,
                                              initial_size=N_van,
                                              growth_rate=G_van),
                  msp.PopulationConfiguration(sample_size=S_lon,
                                              initial_size=N_lon,
                                              growth_rate=G_lon),
                  msp.PopulationConfiguration(sample_size=S_par,
                                              initial_size=N_par,
                                              growth_rate=G_par)
                  ]
#                msp.PopulationConfiguration(sample_size=S_riv,
#                                              initial_size=N_riv,
#                                              growth_rate=G_riv)
#                  ]

    return(pop_config)


