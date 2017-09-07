#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:58:30 2017
stats2abc.py
@author: scott
"""


def obs_stats(vcfin):
    """Take vcfs and calculates observed stats for ABC analysis

    Parameters
    ------
    vcfin
    pops
    downsample

    Returns
    ------
    sfs-fold: dadi
    jsfs: dadi
    3dsfs: dadi
    4dsfs: dadi
    theta_w: dadi
    segsites: dadi
    fst: dadi
    pi: dadi
    r2: popsizeabc
    f2: allel
    f3: allel
    f4: allel
    roh: allel

    """
    vcf2dadi > masking > sfs, jsfs, 3d, 4d, theta, S, Fst, pi
    popsize.vcf2stats > ld bins, maybe IBS
    scikit: f2, f3, f4, RoH

