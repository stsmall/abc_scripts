#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 17:21:53 2018

anopSims_main.py --cfg FILE.cfg --model 1 -i 1

@author: stsmall
"""
from __future__ import print_function
from __future__ import division
import sys
import os
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
import multiprocessing
from itertools import combinations
# check versions
assert sys.version_info[0] >= 3
assert allel.__version__ == "1.1.10"

parser = argparse.ArgumentParser()
parser.add_argument('-cfg', "--configFile", required=True,
                    help="config file")
parser.add_argument('-m', "--model_ix", type=int, required=True,
                    help="model index")
parser.add_argument('-i', "--iterations", type=int, default=1,
                    help="number of iterations")
parser.add_argument("--filet", type=str, help="path to filet exe")
parser.add_argument("--ms", type=str, help="path to ms")
args = parser.parse_args()


def writeABC(stats, seed, scrmline, params, parlist, ix, block, filetpath,
             filet=True, MP=False, nprocs=10):
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
    if nprocs > multiprocessing.cpu_count():
        nprocs = multiprocessing.cpu_count()
    jsfslist = stats.jsfsStats(rand=True)
    jsfstotal = np.sum(jsfslist, axis=1)
    props = [j/jsfstotal[i] for i, j in enumerate(jsfslist)]
    jsfs = " ".join(map(str, np.concatenate(props).ravel()))
    asfslist = stats.asfsStats(rand=True)
    asfs = " ".join(map(str, [i for t in asfslist for i in t]))
    if filet:
        if MP:
            print("filet")
            pool = multiprocessing.Pool(nprocs)
            argslist = []
            for pop1, pop2 in combinations(stats.pops, 2):
                argslist.append([pop1, pop2, block, filetpath])
            filetlist = pool.map(stats.filetStatsMP, argslist)
            pool.close()
            filetstats = " ".join(map(str, np.concatenate(filetlist).ravel()))
        else:
            filet_list = stats.filetStats(block, filetpath)
            filetstats = " ".join(map(str, np.concatenate(filet_list).ravel()))
    x = scrmline.split()
    theta = x[4]
    rho = x[6]
    par = []
    if ix < 7:
        mMax_ix = [i for i, p in enumerate(parlist) if "mMax" in p]
        mIso_ix = [i for i, p in enumerate(parlist) if "mIso" in p]
        if mMax_ix:
            rmv = min(mMax_ix + mIso_ix)
            params = params[0:rmv]
    else:
        pass  # model 7 and 8 are for inferring mig/iso
    for i in params:
        if type(i) is list:
            par.extend(i)
        else:
            par.append(i)
    if ix == 1:
        # 'ej21 ej45 NA ej35 NA ej15 NA NA NA NA NA'
        par.insert(2, 'NA')
        par.insert(4, 'NA')
        nalist = 'NA ' * 5
    elif ix == 2:
        # 'ej21 ej45 NA ej35 NA ej15 NA NA NA NA NA'
        par.insert(2, 'NA')
        par.insert(4, 'NA')
        nalist = 'NA ' * 5
    elif ix == 3:
        # 'ej21 ej45 NA ej35 NA ej15 NA NA NA NA NA'
        par.insert(2, 'NA')
        par.insert(4, 'NA')
        nalist = 'NA ' * 4
    elif ix == 4:
        # 'ej21 ej45 NA es3 esa ej15 NA NA NA NA NA'
        par.insert(2, 'NA')
        nalist = 'NA ' * 5
    elif ix == 5:
        # 'ej21 es4 esa4 es3 esa3 ej15 NA NA NA NA NA'
        nalist = 'NA ' * 5
    elif ix == 6:
        # 'NA es4 esa4 es3 es3a ej15 ej25 em m12 m21 em'
        nalist = None
        par.insert(0, 'NA')
    elif ix == 7:
        # version of model 6 allowing for a comparison with model 8
        # 'NA es4 esa4 es3 es3a ej15 ej25 em m12 m21 em'
        nalist = 'NA ' * 10
        par.insert(0, 'NA')
    elif ix == 8:
        # version of model 6 to infer mMax and mIso
        # 'NA es4 esa4 es3 es3a ej15 ej25 em m12 m21 em mMax35 mIso35 mMax45 mIso45 mMax25 mIso45'
        # cut -d" " -f 4-26
        nalist = None
        par.insert(0, 'NA')
    else:
        pass

    if nalist is None:
        fabc = "{} {} {} {} {} {} {} {}".format(seed, theta, rho, ix,
                                                " ".join(map(str, par)),
                                                asfs, jsfs,
                                                filetstats)
    else:
        fabc = "{} {} {} {} {} {} {} {} {}".format(seed, theta, rho, ix,
                                                   " ".join(map(str, par)),
                                                   nalist.rstrip(), asfs, jsfs,
                                                   filetstats)
    return(fabc)


def simulate(model, demodict, ix, parfx, parlist, thetaarray, rhoarray, ms):
    """Runs simulations by building 1 file with lines equal to it[eration]s.
    Each line is then executed in the shell and output is printed to a file
    with parameter choices in the file name. The loci parameter controls the
    size of the fragment.
    Example: its = 10,000 will produce a file with 10,000 lines that
    when executed will create 10,000 output files.
    """
    # simulate
    Ne = int(np.round((np.random.choice(thetaarray)/(4*model["mutation_rate"]))))
    theta = 4*Ne*model["mutation_rate"] * model["contig_length"]
    rho = np.random.choice(rhoarray) * model["contig_length"]
    npops = len(model["sampleSize"])
    subpops = "-I {} {} 0".format(npops,' '.join(map(str, (model["sampleSize"]))))
    nesubpops = ["-n {} {}".format(i+1, ne/Ne) for i, ne in enumerate(model["initialSize"])]
    gsubpops = ["-g {} {}".format(i+1, g*4*Ne) for i, g in enumerate(model["growthRate"])]
    miglist = model["migMat"].tolist()
    migmat = [val for sublist in miglist for val in sublist]
    params = []
    for px in parfx:
        if type(px) is list:
            params.append([p() for p in px])
        else:
            params.append(px())
    m = Model()
    dem_events = m.genDem(ix,
                          npops,
                          params,
                          demodict,
                          parlist,
                          Ne)
    # all to dict
    ms_params = {
                'ms': ms,
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
    scrm_base = ("{ms} {nhaps} {loci} -t {theta} -r {rho} {basepairs} "
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
    mdir = os.path.abspath(args.configFile)
    mfile = os.path.split(mdir)[0]
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
    migration_matrix = np.genfromtxt(os.path.join(mfile, migFile), delimiter=",")
    assert len(sampleSize) == migration_matrix.shape[0]
    #
    sh = "parameters"
    thetaFile = config.get(sh, "theta_distribution")
    thetaarray = np.loadtxt(os.path.join(mfile, thetaFile))
    rhoFile = config.get(sh, "rho_distribution")
    rhoarray = np.loadtxt(os.path.join(mfile, rhoFile))
    paramFile = config.get(sh, "params")
    parfx, parlist, demodict = drawParams(os.path.join(mfile, paramFile))
    # start functions
    ix = args.model_ix
    model = {"contig_length": contiglen,
             "recombination_rate": recombRate,
             "mutation_rate": mutationRate,
             "sampleSize": sampleSize,
             "initialSize": initialSize,
             "growthRate": growthRate,
             "migMat": migration_matrix,
             "loci": loci
             }
    its = args.iterations
    abcfile = os.path.join(mfile, "abc.{}.{}.out".format(ix, np.random.rand()))
    with open(abcfile, 'w') as f:
        for i in range(its):
            # =================================================================
            # Simulations
            # =================================================================
            mscmd, params = simulate(model,
                                     demodict,
                                     ix,
                                     parfx,
                                     parlist,
                                     thetaarray,
                                     rhoarray,
                                     args.ms)
            # =================================================================
            # Parse
            # =================================================================
            msout = subprocess.Popen(mscmd, shell=True, stdout=subprocess.PIPE)
            print(mscmd)
            print("\nsim complete ... reading file")
            gtlist, pos, pops, block, msline, seed = read_msformat(msout)
            # gtlist,pos,pops,block,msline,seed = read_msformat_file(base)
            # =================================================================
            # Stats
            # =================================================================
            stats = SimStats(gtlist, pops, pos)
            # =================================================================
            # Write to File
            # =================================================================
            fabc = writeABC(stats,
                            seed,
                            msline,
                            params,
                            parlist,
                            ix,
                            block,
                            args.filet)
            f.write("{}\n".format(fabc))
