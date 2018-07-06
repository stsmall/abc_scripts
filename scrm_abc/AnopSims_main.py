#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 17:21:53 2018

anopSims_main.py --cfg FILE.cfg --model 1

@author: stsmall
"""
from __future__ import print_function
from __future__ import division
import sys
import numpy as np
import allel
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
from anopParams import drawParams
from anopModels import Model
from anopStats import SimStats
from parse_ms import read_msformat
import argparse
import subprocess
# check versions
assert sys.version_info[0] >= 3
assert allel.__version__ == "1.1.10"

parser = argparse.ArgumentParser()
parser.add_argument('-cfg', "--configFile", required=True,
                    help="config file")
parser.add_argument('-m', "--model_ix", type=int, required=True,
                    help="model index")
parser.add_argument('-i', "--iterations", type=int, required=True,
                    help="number of iterations")
parser.add_argument("--filetpath", type=str, help="path to filet exe")
args = parser.parse_args()


def writeABC(stats, seed, scrmline, params, ix, block, filetpath, filet=True):
    """Prints results of simulations and stats to text file

    Parameters
    ------
    stats: Obj of class SimStats
    seed: int
    scrmline: str
    params: list
    modelix: int

    Returns
    ------
    stats_file: openFile

    """
    jsfslist = stats.jsfsStats(fold=False)
    jsfstotal = np.sum(jsfslist, axis=1)
    props = [j/jsfstotal[i] for i, j in enumerate(jsfslist)]
    jsfs = " ".join(map(str, np.concatenate(props).ravel()))
    # jsfs = " ".join(map(str, props))
    tots = " ".join(map(str, jsfstotal))
    asfslist = stats.asfsStats(fold=False)
    asfs = " ".join(map(str, [i for t in asfslist for i in t]))
    if filet:
        filet_list = stats.filetStats(block, filetpath)
        filetstats = " ".join(map(str, np.concatenate(filet_list).ravel()))
    f = open("abc.{}.{}.out".format(ix, seed), 'w')
    x = scrmline.split()
    theta = x[4]
    rho = x[6]
    if ix == 1:
        nalist = 'NA ' * 8
    elif ix == 2:
        nalist = 'NA ' * 8
    elif ix == 3:
        nalist = 'NA ' * 8
    elif ix == 4:
        nalist = 'NA ' * 8
    elif ix == 5:
        nalist = 'NA ' * 8
    elif ix == 6:
        nalist = 'NA ' * 8
    else:
        pass
    f.write("{} {} {} {} {} {} {} {} {} {}".format(seed, theta, rho, ix,
            " ".join(map(str, params)), nalist.rstrip(), asfs, jsfs, tots, filetstats))
    f.close()
    return(None)


def simulate(model, demodict, ix, parfx, parlist, thetaarray, rhoarray):
    """Runs scrm by building 1 file with lines equal to it[eration]s. Each line
    is then executed in the shell and output is printed to a file with
    parameter choices in the file name. The loci parameter controls the size of
    the fragment.
    Example: its = 10,000 will produce a file with 10,000 lines that
    when executed will create 10,000 output files.
    """
    # simulate
    Ne = int(np.round((np.random.choice(thetaarray)/(4*model["mutation_rate"]))))
    theta = 4*Ne*model["mutation_rate"] * model["contig_length"]
    rho = np.random.choice(rhoarray) * model["contig_length"]
    subpops = "-I {} {} 0".format(len(model["sampleSize"]),' '.join(map(str, (model["sampleSize"]))))
    nesubpops = ["-n {} {}".format(i+1, ne/Ne) for i, ne in enumerate(model["initialSize"])]
    gsubpops = ["-g {} {}".format(i+1, g*4*Ne) for i, g in enumerate(model["growthRate"])]
    miglist = model["migMat"].tolist()
    migmat = [val for sublist in miglist for val in sublist]
    params = [p() for p in parfx]
    # avoid exactly equal times by small perturbation
    while len(params) != len(set(params)):
        params = [i + np.random.randint(-10, 10) for i in params]
    m = Model()
    dem_events = m.genDem(ix, params, demodict, parlist, Ne, model["mMax"], model["mIso"])
    # all to dict
    ms_params = {
                'nhaps': sum(model["sampleSize"]),
                'loci': model["loci"],
                'theta': theta,
                'rho': rho,
                'basepairs': model["contig_length"] + 1,
                'subpops': subpops,
                'ne_subpop': " ".join(nesubpops),
                'growth_subpop': " ".join(gsubpops),
                'migmat': " ".join(map(str, migmat)),
                'demo': " ".join(dem_events)
                 }
    # scrm command line
    scrm_base = ("scrm {nhaps} {loci} -t {theta} -r {rho} {basepairs} "
                 "{subpops} {ne_subpop} {growth_subpop} -ma {migmat} {demo} "
                 "-p 12")
    mscmd = scrm_base.format(**ms_params)
    return(mscmd, params)


if __name__ == "__main__":
    # =========================================================================
    # Config parsing
    # =========================================================================
    config = configparser.ConfigParser()
    config.read(args.configFile)
    #
    sh = "simulation"
    contiglen = config.getint(sh, "contiglen")
    loci = config.getint(sh, "loci")
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
    assert len(sampleSize) == migration_matrix.shape[0]
    #
    mMax = config.getfloat(sh, "mMax")
    mIso = config.getint(sh, "mIso")
    #
    sh = "parameters"
    thetaFile = config.get(sh, "theta_distribution")
    thetaarray = np.loadtxt(thetaFile)
    rhoFile = config.get(sh, "rho_distribution")
    rhoarray = np.loadtxt(rhoFile)
    paramFile = config.get(sh, "params")
    parfx, parlist, demodict = drawParams(paramFile)
    # start functions
    ix = args.model_ix
    model = {"contig_length": contiglen,
             "recombination_rate": recombRate,
             "mutation_rate": mutationRate,
             "sampleSize": sampleSize,
             "initialSize": initialSize,
             "growthRate": growthRate,
             "migMat": migration_matrix,
             "loci": loci,
             "mMax": mMax,
             "mIso": mIso
             }
    for i in range(args.iterations):
        # =====================================================================
        # Simulations
        # =====================================================================
        mscmd, params = simulate(model,
                                 demodict,
                                 ix,
                                 parfx,
                                 parlist,
                                 thetaarray,
                                 rhoarray)
        # =====================================================================
        # Parse
        # =====================================================================
        print(mscmd)
        msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
        print("\nsim complete ... reading file")
        gtlist, pos, pops, block, scrmline, seed = read_msformat(msout)
        # gtlist, pos, pops, block, scrmline, seed = read_msformat_file(basename)
        # =====================================================================
        # Stats
        # =====================================================================
        print("calculating stats")
        stats = SimStats(gtlist, pops, pos)
        # =====================================================================
        # Write to File
        # =====================================================================
        print("writing stats")
        writeABC(stats, seed, scrmline, params, ix, block, args.filetpath)
