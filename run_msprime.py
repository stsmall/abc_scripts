#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""ABC analysis using LD, jSFS, AF-IBS
requires:
    dadi
    msprime
    scikit-allel
Input:
    ped file
    dadi in file (vcf2dadi.py)
    vcf file (1 for each population)
    chrlen file (chrName length)
Output:
    pop.sfs
    pop1-pop2.jsfs
    pop.ld
    *pop.afibs (Theunert MBE 2012)
Future Additions:
    muage scikit-allel  # Mathieson and McVean PLosGen 2014, quantiles?
    rarecoal: Schiffels NatCom 2016  # distance to doubletons

"""
from __future__ import print_function
from __future__ import division
import allel
import numpy as np
import msprime as msp
from collections import defaultdict
from abc_obsstats import SFSstats
from abc_obsstats import Popsizeabc
from abs_simstats import Simstats
import abc_models
import argparse
from gwas import summary_stat as ss
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
assert msp.__version__ == "0.4.0"
assert allel.__version__ == "1.1.9"

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", type=str, required=True,
                    help='path to vcf file')
parser.add_argument('-p', "--pedFile", type=str, required=True,
                    help="path to ped file")
parser.add_argument('-c', "--chrlen", type=str, required=True,
                    help="path to file of chromosome lengths and IDs")
parser.add_argument('-d', "--dadiFile", type=str, required=True,
                    help="path to dadi infile")
parser.add_argument('-m', "--model", type=str, required=True,
                    help="model name")
parser.add_argument('-i', "--configFile", type=str, required=True,
                    help="path to config infile")
parser.add_argument('-o', "--outfile", type=str, required=True,
                    help="outfile name")
args = parser.parse_args()


def vcf2obsstats(vcfFile, chrlist, dadiFile, pedFile):
    """Calculates summary statistics from VCF file
    """
    # sfs and jsfs
    sfs = SFSstats()
    sfs.sfs_obs(dadiFile, pedFile, mask=True, fold=True)
    jsfsdict = sfs.jSFS_summary()
    sfs.summarywrite(jsfsdict)
    # ld and AF-IBD
    peddict = defaultdict(list)
    with open(pedFile, 'r') as ped:
        for line in ped:
            if line.strip():
                x = line.strip().split()
                peddict[x[0]].append(x[1])
            else:
                continue
    for pix, pop in enumerate(peddict.keys()):
        p = Popsizeabc(pop)
        p.ldstats(vcfFile, chrlist, pedFile, configdict, pix)
        p.printstats()
    return(p.ldints)


def sims2stats(tree_sequence, ldints, configdict):
    """Calculates summary statistics from msprime simulation object
    """
    posdict = defaultdict(list)
    hapdict = defaultdict(list)
    countdict = defaultdict(list)
    npop = len(configdict["sample_size"])
    for ts in tree_sequence:
        muts = ts.get_num_mutations()
        if muts > 0:
            statout = Simstats()
            statlist = []
            # LD feed to PopSizeABC
            hapdict, posdict, countdict = statout.ld_stats(ts, ldints, configdict, hapdict, posdict)
            # JSFS feed to allel
            jsfs = statout.jsfs(ts, npop)
            statlist.append(jsfs)
            # AF-IBD, custom script
            afibs = statout.afibs(ts, configdict["contig_length"], npop)
            statlist.append(afibs)
    return(statlist, hapdict, posdict, countdict)


def countstats(sizelist):
    """
    """
    stats = 0
    pops = len(sizelist)
    stats += sum(sizelist) / 2  # sfs
    stats += ((pops*(pops-1)) / 2) * 23  # jsfs
    stats += sum(sizelist) / 2  # AFIBS
    return(stats)


def runsims(model, configdict, iterations, outfile, ldints):
    """Runs simulations using model file and config file
    """
    # set population parameters
    pop_config = []
    for i, pop in enumerate(configdict["sample_size"]):
        size = pop
        init = configdict["initial_size"][i]
        grow = 0
        pop_config.append(msp.PopulationConfiguration(sample_size=size,
                                                      initial_size=init,
                                                      growth_rate=grow))
    # set model
    model = getattr(__import__('demographic_models'), model)
    # print demography
    dp = msp.DemographyDebugger(population_configurations=pop_config,
                                demographic_events=model)
    dp.print_history()
    # get prior and stat lens
    dem_list, priors = abc_models.model(configdict)
    stats = countstats(configdict["sample_size"])
    # run simulate
    param_arr = np.zeros(iterations, len(priors))
    stats_arr = np.zeros(iterations, len(stats))
    rmin = configdict["recombination_rate"][0]
    rmax = configdict["recombination_rate"][2]
    for i in range(iterations):
        dem_list, priors, N_ANC = abc_models.model(configdict)
        rrate = np.random.uniform(rmin, rmax)
        priors[0] = rrate
#        Ne=configdict["Ne"]
#        Ne=N_ANC
        stats = []
        tree_sequence = msp.simulation(
                                       population_configurations=pop_config,
                                       demographic_events=dem_list,
                                       length=configdict["contig_length"],
                                       mutation_rate=configdict["muation_rate"],
                                       recombination_rate=rrate,
                                       random_seed=np.random.randint(1, 1000001),
                                       num_replicates=configdict["nb_seg"]
                                       )

        statlist, hapdict, posdict, countdict = sims2stats(tree_sequence, ldints, configdict)
        stats.append(statlist)
        ldlist = []
        aflist = []
        for pop, n in range(len(configdict["sample_size"])):
            res_afs = ss.histo(countdict[pop], n/2)
            geno_list = ss.hap_to_geno(hapdict[pop])
            res_ld_zyg = ss.distrib_zyg_r2(posdict[pop], geno_list, ldints)
            ldlist.append(res_ld_zyg)
            aflist.append(res_afs)
        afarr = np.array(aflist)
        ldarr = np.array(ldlist)
        stats_arr[i] = np.concatenate(np.mean(np.array(stats), axis=0), afarr, ldarr)
        param_arr[i] = priors
    # write stats and params
    np.savetxt("{}.params".format(outfile), param_arr)
    np.savetxt("{}.stats".format(outfile), stats_arr)
    return(None)


if __name__ == "__main__":
    # stats from VCF
    vcfFile = args.vcfFile
    # config for simulations
    config = configparser.ConfigParser()
    config.read(args.configfile)
    sh = "simulation"
    pops = config.getint(sh, "num_pops")
    sampleSize = list(map(int, config.get(sh, "sample_size").split(",")))
    initialSize = list(map(int, config.get(sh, "initial_size").split(",")))
    try:
        assert len(sampleSize) == pops
        assert len(initialSize) == pops
    except AssertionError:
        print("each popualtion must have intial values")
    iterations = config.getint(sh, "iterations")
    contiglen = config.getint(sh, "contiglen")
    effectiveSize = config.getint(sh, "effective_size")
    recombRate = config.getfloat(sh, "recombination_rate")
    recombMin = config.getfloat(sh, "recombinmin")
    recombMax = config.getfloat(sh, "recombinmax")
    mutationRate = config.getfloat(sh, "mutation_rate")
    nb_times = config.getint(sh, "nb_times")
    tmax = config.getint(sh, "tmax")
    a = config.getfloat(sh, "a")
    pererr = config.getint(sh, "pererr")
    nb_seg = config.getint(sh, "nb_seg")
    mac_ld = config.getint(sh, "mac_ld")
    mac = config.getint(sh, "mac")
    configdict = {"sample_size": sampleSize,
                  "initial_size": initialSize,
                  "contig_length": contiglen,
                  "Ne": effectiveSize,
                  "recombination_rate": (recombMin, recombRate, recombMax),
                  "mutation_rate": mutationRate,
                  "nb_times": nb_times,
                  "tmax": tmax,
                  "a": a,
                  "pererr": pererr,
                  "nb_seg": nb_seg,
                  "mac_ld": mac_ld,
                  "mac": mac}
    ldints = vcf2obsstats(vcfFile, args.pedFile, args.chrlen, args.dadiFile, configdict)
    runsims(args.model, configdict, iterations, args.outfile, ldints)

####test
dem_list=[msp.MassMigration(time=10, source=0, destination=1, proportion=1.0), msp.MassMigration(time=100, source=1, destination=2, proportion=1.0)]
popcfg = [msp.PopulationConfiguration(sample_size=10, initial_size=100), msp.PopulationConfiguration(sample_size=10, initial_size=1000),msp.PopulationConfiguration(sample_size=5, initial_size=500)]
tree_sequence = msp.simulate(population_configurations=popcfg, Ne=1000, length=1e4, recombination_rate=2e-8, mutation_rate=2e-8, demographic_events=dem_list)
