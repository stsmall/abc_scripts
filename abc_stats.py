#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 13:05:00 2020
@author: Scott T. Small

Main script for generating statistics for ABC and ML training.
depends: parse_sims.py, sims_stats.py, obs_stats.py

Example
-------

abc_stats.py sim --infile ms.out --outfile ms.stats --pops 4 --pairs 0-1 0-2 0-3
    --stats sfs jsfs --mask --gff foo.gff --mode split-run

    --mask for FILET, expects names as 0-1.mask.txt 0-2.mask.txt ...
    --gff to avoid genes

abc_stats.py obs --infile vcf/h5 --fasta foo.mask.fa --gff foo.gff --pops 4
    --pairs 0-1 0-2 0-3 --stats filet --anc_fasta anc.fasta

    --fasta masked fasta for FILET masking
    --anc_fasta for FILET polarizing and unfolded SFS

Notes
-----

Creates summary stats from coalescent simulations: ms, msmove, discoal, msprime
There are two main modes: sims and obs

 Mode 'sims' is for ms-style formats but will also use msprime treeseqs (TODO)

 Mode 'obs' is for generating the same stats but with a starting vcf. Take a look
    at generating mask files from FILET and diploshic


Example
-------

    $ python abc_stats.py --infile --pairs --stats --out_file --split

"""
import sys
import argparse
from itertools import combinations
import multiprocessing
import glob
import os
import numpy as np
from timeit import default_timer as timer
from sim_parse import ms_parse
from sim_parse import split2pairs
from sim_stats import SumStats


def pair_split(msdict, out_file, pairs):
    """Generate pairs for splitting.

    Parameters
    ----------
    infile : TYPE
        DESCRIPTION.
    outfile : TYPE
        DESCRIPTION.
    pairs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    for p in pairs:
        p1, p2 = p.split("-")
        p1, p2 = int(p1), int(p2)
        split2pairs(msdict, out_file, p1, p2)
    return None


def write_stats_out(out_file, pairs, asfsdict, jsfsdict, afibsdict, filetdict):
    """Write summary stats to a file.

    Parameters
    ----------
    out_file : TYPE
        DESCRIPTION.
    asfsdict : TYPE
        DESCRIPTION.
    jsfsdict : TYPE
        DESCRIPTION.
    afibsdict : TYPE
        DESCRIPTION.
    filetdict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # asfs_default = "afS1 afS2_2"
    # jsfs_default = "jfS1 jfS2 jfS3 jfS4 jfS5 jfS6 jfS7 " \
    #                "jfS8 jfS9 jfS10 jfS11 jfS12 jfS13 "   \
    #                "jfS14 jfS15 jfS16 jfS17 jfS18 jfS19 " \
    #                "jfS20 jfS21 jfS22 jfS23"
    filet_default = "pi1 hetVar1 ss1 private1 tajd1 ZnS1 pi2 hetVar2 ss2 "    \
                    "private2 tajd2 ZnS2 Fst snn dxy_mean dxy_min gmin zx "  \
                    "dd1 dd2 ddRank1 ddRank2"
    # note that this FILET does not include haplotype or IBS stats

    # prep headers
    afibs_header = []
    asfs_header = []
    jsfs_header = []
    filet_header = []

    sub_pops = set([int(item) for x in [x.split("-") for x in pairs] for item in x])

    if afibsdict:
        for p in sub_pops:
            pass
            #afibs_header.append()
    if asfsdict:
        for p in sub_pops:
            asfs_header.append(f"afS1_{p}")
        for p in sub_pops:
            asfs_header.append(f"afS2_{p}")
    if jsfsdict:
        for p in pairs:
            for j in range(1, 24):
                jsfs_header.append(f"jfS{j}_{p}")
    if filetdict:
        for p in pairs:
            for f in filet_default.split():
                filet_header.append(f"{f}_{p}")

    header = afibs_header + asfs_header + jsfs_header + filet_header

    # write
    for d in [afibsdict, asfsdict, jsfsdict, filetdict]:
        if d:
            reps = list(d.keys())
            break
    if os.path.exists(out_file):
        out_file = f"{out_file}_{round(np.random.rand(), 2)}"
    with open(out_file, "w") as out:
        out.write(f"{' '.join(header)}\n")
        for rep in reps:
            line = []
            for d in [afibsdict, asfsdict, jsfsdict]:
                if d:
                    line.append(d[rep])
            if filetdict:
                line.append(" ".join(filetdict[rep]))
            out.write(f"{' '.join(line)}\n")
    return None


def calc_simstats(msdict, pairs, stats, filetpath, nprocs, window=10000):
    """Calculate summary stats from ms.

    Parameters
    ----------
    msdict : TYPE
        DESCRIPTION.
    pairs : TYPE
        DESCRIPTION.
    stats : TYPE
        DESCRIPTION.
    filetpath : TYPE
        DESCRIPTION.
    nprocs : TYPE
        DESCRIPTION.
    window : TYPE, optional
        DESCRIPTION. The default is 10000.

    Returns
    -------
    asfsdict : TYPE
        DESCRIPTION.
    jsfsdict : TYPE
        DESCRIPTION.
    afibsdict : TYPE
        DESCRIPTION.
    filetdict : TYPE
        DESCRIPTION.

    """
    popconfig = msdict["pops"]
    basepairs = msdict["basepairs"]
    pos_list = msdict["pos"]
    hap_list = msdict["haps"]
    # format for allel
    hap_arrT = []
    for hap in hap_list:
        if len(hap) > 1:
            hap_tmp = []
            for h in hap:
                hap_tmp.append(h.transpose())
            hap_arrT.append(hap_tmp)
        else:
            hap_arrT.append(hap[0].transpose())
    # stats class
    sub_pops = set([int(item) for x in [x.split("-") for x in pairs] for item in x])
    pops = [popconfig[i] for i in sub_pops]
    sum_stats = SumStats(hap_arrT, pos_list, pops)
    # calc stats
    if "sfs" in stats:
        asfsdict = sum_stats.asfsStats(fold=False, rand=True, randn=100000)
        # asfs = " ".join(map(str, [i for t in aSFS for i in t]))
    else:
        asfsdict = {}
    if "jsfs" in stats:
        jsfsdict = sum_stats.jsfsStats(pairs, rand=True)
    else:
        jsfsdict = {}
    if "afibs" in stats:
        afibsdict = sum_stats.afibs(basepairs, fold=False)
    else:
        afibsdict = {}
    if "filet" in stats:
        if nprocs == 1:
            filetdict = sum_stats.filetStats(basepairs, filetpath, window)
        else:
            filet_list = []
            argslist = []
            filetdict = {}

            # check that there are not more requested than available
            if nprocs > multiprocessing.cpu_count():
                nprocs = multiprocessing.cpu_count()
            # set pool and map
            pool = multiprocessing.Pool(nprocs)
            for pop1, pop2 in combinations(sum_stats.pops, 2):
                argslist.append([pop1, pop2, basepairs, filetpath, window])
            filet_list.append(pool.map(sum_stats.filetStatsMP, argslist))
            pool.close()
            # resize and zip
            filet_zip = list(zip(*filet_list[0]))
            for r in range(len(filet_zip)):
                filetdict[r] = filet_zip[r]
    else:
        filetdict = {}
    return asfsdict, jsfsdict, afibsdict, filetdict


# def calc_observedstats(sample_names, pairs, vcf, gvcf, anc_fa, frac_miss, qual, repeat_gff):
#     """Calculate observed statistics.

#     Parameters
#     ----------
#     vcf : TYPE
#         DESCRIPTION.
#     gvcf : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     """
#     from obs_stats import *
#     import normalizePgStats as norm_stats

#     # sfs, jsfs, afibs
#     # TODO: use functions from diploshic for importing VCF to allel
#     # use sim_stats functions once in allel format

#     # FILET
#     # this could likely be done with SnakeMake since it is just passing files around
#     for arm in chrom:
#         windows = make_coords(chrom, chr_len, start, end, window_size, step)
#     vcfs, masks = build_mask(pairs, sample_names, vcf, gvcf, repeat_gff, ref_fa, prcnt_miss, qual)
#     fastas = vcf_to_fasta(vcf_path, vcfs, mask_path, masks)
#     for pair in fastas:
#         p1, p2, arm = pair.split("-")
#         f"pgStatsBedSubpop_forML {p1}.{p1}-{p2}.{arm}.masked.fasta {p2}.{p1}-{p2}.{arm}.masked.fasta"
#         f"{anc_fa} windows.{arm}.txt {frac_miss} | norm_stats {window_size}"


def parse_args(args_in):
    """Parse args.

    Parameters
    ----------
    args_in : TYPE
        DESCRIPTION.

    Returns
    -------
    argsDict : TYPE
        DESCRIPTION.

    """
    parser = argparse.ArgumentParser(description="calculate summary statistics"
                                     " from simulated data or observed data")
    parser._positionals.title = f"enter 'python {sys.argv[0]} modeName -h' for modeName's help message"
    subparsers = parser.add_subparsers(help='sub-command help')
    parser_a = subparsers.add_parser('sim', help="Generate stats from sim data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_b = subparsers.add_parser('obs', help="Generate stats from data in a VCF",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_a.set_defaults(mode='sim')
    parser_a._positionals.title = "required arguments"
    parser_a.add_argument('ms_file', help="path to simulation output file"
                          "must be same format used by Hudson\'s ms")
    parser_a.add_argument('out_file', help="path to file where stats or outfile"
                          "will be written")
    parser_a.add_argument("--split", action="store_true",
                          help="split is preprocessing")
    parser_a.add_argument("--pairs", nargs='+', default="all",
                          help="list of pairs separate by hyphen, 0 indexed")
    parser_a.add_argument("--stats", nargs='+', default="filet",
                          choices=["sfs", "jsfs", "filet", "afibs"],
                          help="which stats to calculate")
    parser_a.add_argument("--nprocs", type=int, default=1, help="try to run"
                          " with parallel. If > 1 will parallel on current machine"
                          " else will run through each file(s) in order")
    parser_a.add_argument('--filet_path', type=str, default="current_dir",
                          help="path to FILET dir w/ programs")

    # calculate stats from a VCF file
    parser_b.set_defaults(mode='obs')
    parser_b._positionals.title = "required arguments"
    parser_b.add_argument('chr_arm_vcf', help="VCF format file containing data"
                          "for our chromosome arm (other arms will be ignored)")
    parser_b.add_argument('chr_arm', help="Exact name of the chromosome arm for"
                          "which feature vectors will be calculated")
    parser_b.add_argument('chr_len', type=int, help="Length of the chromosome arm")
    parser_b.add_argument('out_file', help="path to file where feature vectors "
                          "will be written")
    parser_b.add_argument("--pairs", nargs='+', default="all",
                          help="list of pairs separate by hyphen, 0 indexed")
    parser_b.add_argument("--stats", nargs='+', default="filet",
                          choices=["sfs", "jsfs", "filet", "afibs", "all"],
                          help="which stats to calculate")
    parser_b.add_argument('--win_size', type=int, default=10000,
                          help="Length of the window")
    parser_b.add_argument('--win_slide', type=int, default=0,
                          help="overlap/slide between windows")
    parser_b.add_argument('--mask_file', default=None,
                          help="Path to a fasta-formatted file that contains"
                          "masking information (marked by \'N\'). If specified,"
                          "simulations will be masked in a manner mirroring"
                          "windows drawn from this file.")
    parser_b.add_argument('--unmasked_frac_cutoff', type=float, default=0.25,
                          help="Minimum fraction of unmasked sites, if masking simulated data")
    parser_b.add_argument('--ancestral_fasta', default=None,
                          help="Path to a fasta-formatted file that contains"
                          " inferred ancestral states")
    parser_b.add_argument('--segment_start', default=None,
                          help="Left boundary of region in which feature vectors"
                          " are calculated (whole arm if omitted)")
    parser_b.add_argument('--segment_end', default=None,
                          help="Right boundary of region in which feature vectors"
                          "are calculated (whole arm if omitted)")
    parser_b.add_argument('--filet_path', type=str, default="current_dir",
                          help="path to FILET dir w/ programs")
    args = parser.parse_args(args_in)
    argsDict = vars(args)

    return argsDict


if __name__ == "__main__":
    argsDict = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    if argsDict["mode"] == "sim":
        ms_file = argsDict["ms_file"]
        out_file = argsDict["out_file"]
        split = argsDict["split"]
        pairs = argsDict["pairs"]
        stats = argsDict["stats"]
        processors = argsDict["nprocs"]
        filet_path = os.path.abspath(argsDict["filet_path"])
    else:
        vcf = argsDict["chr_arm_vcf"]
        chrom = argsDict["chr_arm"]
        chrom_len = argsDict["chr_len"]
        out_file = argsDict["out_file"]
        win_size = argsDict["win_size"]
        win_slide = argsDict["win_slide"]
        mask_file = argsDict["mask_file"]
        unmasked_frac = argsDict["unmasked_frac_cutoff"]
        unmasked_geno = argsDict["unmasked_geno_cutoff"]
        anc_fasta = argsDict["ancestral_fasta"]
        seg_start = argsDict["segment_start"]
        seg_end = argsDict["segment_end"]
        pairs = argsDict["pairs"]
        stats = argsDict["stats"]
        filet_path = argsDict["filet_path"]
    # =========================================================================
    #  Main executions
    # =========================================================================
    if argsDict["mode"] == "sim":
        ms_file = os.path.abspath(ms_file)
        if os.path.isdir(ms_file):
            msfiles = glob.glob("*.ms")
        else:
            msfiles = [ms_file]
        for ms in msfiles:
            msdict = ms_parse(ms)
            if pairs == "all":
                npairs = []
                p = list(range(len(msdict["pops"])))
                for p1, p2 in combinations(p, 2):
                    npairs.append(f"{p1}-{p2}")
            elif "-" not in pairs[0]:
                npairs = []
                for p1, p2 in combinations(pairs, 2):
                    npairs.append(f"{p1}-{p2}")
            else:
                npairs = pairs
            if split:
                pair_split(msdict, out_file, npairs)
            else:
                #start = timer()
                asfsdict, jsfsdict, afibsdict, filetdict = calc_simstats(msdict, npairs, stats, filet_path, processors)
                #end = timer()
                #print(end - start)
    elif argsDict["mode"] == "obs":
        pass
        # calc_observedstats(sample_names, pairs, vcf, gvcf, anc_fa, frac_miss, qual, repeat_gff)
    # write stats to a flile
    write_stats_out(out_file, npairs, asfsdict, jsfsdict, afibsdict, filetdict)
