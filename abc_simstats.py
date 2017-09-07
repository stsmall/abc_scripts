#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:49:19 2017
This python code directly calculates statistics from ms input. Preferrably this
should be used with msprime
ms2stats.py
@author: scott
"""
import argparse
import numpy as np
import msprime as msp
import dadi
import allel
from libsequence.msprime import makeSimData

parser = argparse.ArgumentParser()
parser.add_argument('-s', "--stats", type=str, required=True,
                    help="name of ingroup")


def msprime_stats(ts):
    """default stats from msprime

    Parameters
    ------
    ts: treesequence, msprime.Simulate()

    Returns
    ------
    f2
    f3
    f4
    pi
    doubletons
    r2: array, r2 statistic of ld in bins from 500-1.5Mb

    """
    # ts_ld = ms.LdCalculator(ts)


def msp2allel(ts):
    """Outputs msprime into format for scikit-allel
    """
    tree_sequence = ts
    shape = tree_sequence.get_num_mutations(), tree_sequence.get_sample_size()
    A = np.empty(shape, dtype="u1")
    for variant in tree_sequence.variants():
        A[variant.index] = variant.genotypes
    return(A)


def allel_stats(ts):
    """Use scikit-allel to calculate stats from numpy arrays
    """
    # RoH
    gt = msp2allel(ts)


def pylib_stats(ts):
    """Use functions from pylibseq to calculate stats, will accept msprime as
    input
    """
    s = make_SimData(ts)


# I would like to avoid dadi if I can, or rewrite to accept numpy arrays from msprime similar to scikit
def dadi_stats():
    """Use dadi functions to calculate SFS and jSFS for input data. Dadi will
    accept ms-style data. msprime > msstype, or use same seed
    """

    fs = dadi.Spectrum ([0, 100, 20, 10, 1, 0])
    # need to mask corners
    theta_w_pop0 = fs.marginalize([1, 2])
    pi_pop0 = fs.marginalize([1, 2])
    tajD_pop0 = fs.marginalize([1, 2])
    segSite_pop0_1 = fs.marginalize([2])
    fst_pop0_1 = fs.marginalize([2])
    # SFS
    fs_pop0 = fs.marginalize([1, 2])
    fs_pop0.mask[1, :] = True
    # jSFS
    jfs_pop0_1 = fs.marginalize([2])