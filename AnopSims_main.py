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
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from anopParams import drawParams
from anopModels import Model
#import anopStats
import argparse
# check versions
assert sys.version_info[0] >= 3
assert msp.__version__ == "0.6.0"
assert allel.__version__ == "1.1.10"

parser = argparse.ArgumentParser()
parser.add_argument('-cfg', "--configFile", required=True,
                    help="config file")
parser.add_argument('-m', "--model_ix", type=int, required=True,
                    help="model index")
parser.add_argument('-i', "--iterations", type=int, required=True,
                    help="number of iterations")
args = parser.parse_args()


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


def writeABC(stats_file, stats_out, params_out):
    """Prints results of simulations and stats to text file

    Parameters
    ------
    stats_file: openFile
    stats_out: list of calculated stats
    params_out: list of input parameters

    Returns
    ------
    stats_file: openFile

    """

    return(stats_file)


def simulate(model, demodict, ix, parfx, parlist, thetaarray, rhoarray):
    """Runs msprime with model
    """
    treelist = []
    # simulate
    for loc in range(model["loci"]):
        Ne = np.random.choice(thetaarray)/(4*model["mutation_rate"])
        params = [p() for p in parfx]
        # avoid exactly equal times by small perturbation
        while len(params) != len(set(params)):
            params = [i + np.random.randint(0, 10) for i in params]
        m = Model()
        dem_events = m.genDem(ix, params, demodict, parlist, Ne)
        recomb = np.random.choice(rhoarray)/(4*Ne)
        checkDemo(model["pop_config"], model["migMat"], dem_events)
        import ipdb;ipdb.set_trace()
        ts = msp.simulate(
                          population_configurations=model["pop_config"],
                          demographic_events=dem_events,
                          Ne=Ne,
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
    N_fun, N_like, N_van, N_par, N_long = initialSize
    # sample size as diploid
    S_fun, S_like, S_van, S_par, S_long = sampleSize
    # intial growth rates
    G_fun, G_like, G_van, G_par, G_long = growthRate
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
                  msp.PopulationConfiguration(sample_size=S_par,
                                              initial_size=N_par,
                                              growth_rate=G_par),
                  msp.PopulationConfiguration(sample_size=S_long,
                                              initial_size=N_long,
                                              growth_rate=G_long)
                  ]
    return(pop_config)


if __name__ == "__main__":
    config = configparser.ConfigParser()
    config.read(args.configFile)
    #
    sh = "simulation"
    contiglen = config.getint(sh, "contiglen")
    loci = config.getint(sh, "loci")
    iterations = config.getint(sh, "iterations")
    recombRate = config.getfloat(sh, "recombination_rate")
    mutationRate = config.getfloat(sh, "mutation_rate")
    effectiveSize = config.getint(sh, "effective_population_size")
    #
    sh = "initialize"
    sampleSize = list(map(int, config.get(sh, "sample_sizes").split(",")))
    initialSize = list(map(int, config.get(sh, "initial_sizes").split(",")))
    growthRate = list(map(float, config.get(sh, "growth_rates").split(",")))
    #
    migFile = config.get(sh, "migration_matrix")
    migration_matrix = np.genfromtxt(migFile, delimiter=",")
    #
    demoFile = config.get(sh, "demographic_file")
    demodict = defaultdict(list)
    with open(demoFile, 'r') as f:
        for line in f:
            key, *val = line.split()
            demodict[int(key)].append([float(val[0]), val[1]])
    #
    sh = "parameters"
    thetaFile = config.get(sh, "theta_distribution")
    thetaarray = np.loadtxt(thetaFile)
    rhoFile = config.get(sh, "rho_distribution")
    rhoarray = np.loadtxt(rhoFile)
    paramFile = config.get(sh, "params")
    parfx, parlist = drawParams(paramFile)
    # start functions
    popcfg = setInitial(sampleSize, initialSize, growthRate)
    model = {"contig_length": contiglen,
             "Ne": effectiveSize,
             "recombination_rate": recombRate,
             "mutation_rate": mutationRate,
             "loci": loci,
             "migMat": migration_matrix.tolist(),
             "pop_config": popcfg}
    # open file
    # statsFile = open("ABCsims.model-{}.txt".format(args.model_ix), 'w')
    for i in range(args.iterations):
        simulate(model,
                 demodict,
                 args.model_ix,
                 parfx,
                 parlist,
                 thetaarray,
                 rhoarray
                 )
##        scrm_base = ("scrm {nhaps} 1 -t {theta} -r {rho} {basepair} "
##                     "-G {exp_growth} -eG {time_growth} 0.0 -SC abs -p"
##                     " {sig_digits} ")
##        mscmd = scrm_base.format(**ms_params)
