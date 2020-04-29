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
import numpy as np
import allel
from itertools import combinations
import multiprocessing
import glob
import os
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
    if type(pairs) == list:
        for p in pairs:
            p1, p2 = p.split("-")
            p1, p2 = int(p1), int(p2)
            split2pairs(msdict, out_file, p1, p2)
    else:
        p = list(range(pairs))
        for p1, p2 in combinations(p):
            split2pairs(msdict, out_file, p1, p2)
    return None


def calc_simstats(msdict, out_file, pairs, stats, procs, mask, filetpath):
    """Calculate summary stats from ms.

    Parameters
    ----------
    msdict : TYPE
        DESCRIPTION.
    outfile : TYPE
        DESCRIPTION.
    pairs : TYPE
        DESCRIPTION.
    stat_list : TYPE
        DESCRIPTION.
    nprocs : TYPE
        DESCRIPTION.
    mask : TYPE
        DESCRIPTION.
    filet_path : TYPE
        DESCRIPTION.

    Returns
    -------
    asfs : TYPE
        DESCRIPTION.
    jsfs : TYPE
        DESCRIPTION.
    afibs : TYPE
        DESCRIPTION.
    filet : TYPE
        DESCRIPTION.

    """
    pops = msdict["popconfig"]
    basepairs = msdict["basepairs"]
    pos_list = msdict["pos"]
    hap_list = msdict["haps"]
    # TODO: pairs use it to remove from popconfig/pops
    hap_arr = [h.transpose() for h in hap_list]
    gt = allel.HaplotypeArray(np.vstack(hap_arr), dtype='i1')
    sum_stats = SumStats(gt, pos_list, pops)

    if any("sfs" in t for t in stats):
        asfslist = sum_stats.asfsStats(rand=True)
        asfs = " ".join(map(str, [i for t in asfslist for i in t]))
    if any("jsfs" in t for t in stats):
        jsfslist = sum_stats.jsfsStats(rand=True)
        jsfstotal = np.sum(jsfslist, axis=1)
        props = [j/jsfstotal[i] for i, j in enumerate(jsfslist)]
        jsfs = " ".join(map(str, np.concatenate(props).ravel()))
    if any("afibs" in t for t in stats):
        afibs = sum_stats.afibs()
    if any("filet" in t for t in stats):
        if procs == 1:
            filet_list = sum_stats.filetStats(basepairs, mask, filetpath)
            filet = " ".join(map(str, np.concatenate(filet_list).ravel()))
        else:
            # check that there are not more requested than available
            if procs > multiprocessing.cpu_count():
                nprocs = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(nprocs)
            argslist = []
            for pop1, pop2 in combinations(sum_stats.pops, 2):
                argslist.append([pop1, pop2, basepairs, mask, filetpath])
            filet_list = pool.map(sum_stats.filetStatsMP, argslist)
            pool.close()
            filet = " ".join(map(str, np.concatenate(filet_list).ravel()))
    return asfs, jsfs, afibs, filet


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
    parser_a = subparsers.add_parser('sim', help="Generate stats from sim data")
    parser_b = subparsers.add_parser('obs', help="Generate stats from data in a VCF")

    parser_a.set_defaults(mode='sim')
    parser_a._positionals.title = "required arguments"
    parser_a.add_argument('ms_file', help="path to simulation output file"
                          "must be same format used by Hudson\'s ms")
    parser_a.add_argument('out_file', help="path to file where stats or outfile will be written")
    parser_a.add_argument("--split", action="store_true",
                          help="split is preprocessing")
    parser_a.add_argument("--pairs", nargs='+', action='append', default="all",
                          help="list of pairs separate by hyphen, 0 indexed")
    parser_a.add_argument("--stats", nargs='+', action='append',
                          choices=["sfs", "jsfs", "filet", "afibs", "all"],
                          required=True, help="which stats to calculate")
    parser_a.add_argument("--processors", type=int, default=1, help="try to run"
                          " with parallel. If > 1 will parallel on current machine"
                          " else will run through each file(s) in order")
    parser_a.add_argument('--mask_file', default=None,
                          help="Path to a fasta-formatted file that contains"
                          "masking information (marked by \'N\'). If specified,"
                          "simulations will be masked in a manner mirroring"
                          "windows drawn from this file.")
    parser_a.add_argument('--filet_path', type=str,
                          help="path to FILET dir w/ programs")


    # calculate stats from a VCF file
    parser_b.set_defaults(mode='obs')
    parser_b._positionals.title = "required arguments"
    parser_b.add_argument('chr_arm_vcf', help="VCF format file containing data"
                          "for our chromosome arm (other arms will be ignored)")
    parser_b.add_argument('chr_arm', help="Exact name of the chromosome arm for"
                          "which feature vectors will be calculated")
    parser_b.add_argument('chr_len', type=int, help="Length of the chromosome arm")
    parser_b.add_argument('out_file', help="path to file where feature vectors will be written")
    parser_b.add_argument('--win_size', type=int, default=10000,
                          help="Length of the window")
    parser_b.add_argument('--win_slide', type=int, default=0,
                          help="overlap/slide between windows")
    parser_a.add_argument('--mask_file', default=None,
                          help="Path to a fasta-formatted file that contains"
                          "masking information (marked by \'N\'). If specified,"
                          "simulations will be masked in a manner mirroring"
                          "windows drawn from this file.")
    parser_a.add_argument('--unmasked_frac_cutoff', type=float, default=0.25,
                          help="Minimum fraction of unmasked sites, if masking simulated data")
    parser_b.add_argument('--ancestral_fasta', default=None,
                          help="Path to a fasta-formatted file that containsvinferred ancestral states")
    parser_b.add_argument('--segment_start', default=None,
                          help="Left boundary of region in which feature vectors"
                          " are calculated (whole arm if omitted)")
    parser_b.add_argument('--segment_end', default=None,
                          help="Right boundary of region in which feature vectors"
                          "are calculated (whole arm if omitted)")
    parser_b.add_argument("--pairs", nargs='+', action='append', default="all",
                          help="list of pairs separate by hyphen, 0 indexed")
    parser_b.add_argument("--stats", nargs='+', action='append',
                          choices=["sfs", "jsfs", "filet", "afibs", "all"],
                          default="filet", help="which stats to calculate")
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
        processors = argsDict["processors"]
        mask_file = argsDict["mask_file"]
        filet_path = argsDict["filet"]
    else:
        vcf = argsDict["chr_arm_vcf"]
        chrom = argsDict["chr_arm"]
        chrom_len = argsDict["chr_len"]
        out_file = argsDict["out_file"]
        target_pop = argsDict["target_pop"]
        sample_names = argsDict["sample_names"]
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

    # =========================================================================
    #  Main executions
    # =========================================================================
    if argsDict["mode"] == "sim":
        if os.path.isdir(ms_file):
            msfiles = glob.glob("*.ms")
        else:
            msfiles = [ms_file]
        for ms in msfiles:
            msdict = ms_parse(ms)
            if split:
                pair_split(msdict, out_file, pairs)
            else:
                calc_simstats(msdict, out_file, pairs, stats, processors, mask_file, filet_path)
    elif argsDict["mode"] == "obs":
        pass
        # calc_observedstats(sample_names, pairs, vcf, gvcf, anc_fa, frac_miss, qual, repeat_gff)
