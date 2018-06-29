#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:05:14 2018
Class for calculating summary statistics from simulations
@author: stsmall
"""

import numpy as np
import allel
import msprime as msp
import bisect
from itertools import combinations


class Simstats:
    """Calculates stats from pos/gt and msprime object
    """
    def __init__(self, sfs=None, jsfs=None, afibs=None):
        self.sfs = sfs
        self.jsfs = jsfs
        self.afibs = afibs

    def jsfsStats(self, treelist, npop, fold=True):
        """Joint site frequency spectrum with scikit-allel
        """
        gtmat = [ts.genotype_matrix() for ts in treelist]
        v = np.vstack(gtmat)
        gt = allel.HaplotypeArray(v)
        pix = []
        jsfslist = []
        pix = [treelist[0].get_samples(pop) for pop in range(npop)]
        for i, j in combinations(pix, 2):
            gtpop1 = gt[:, i]
            gtpop2 = gt[:, j]
            ac1 = gtpop1.count_alleles()
            ac2 = gtpop2.count_alleles()
            if fold:
                # pad for allel as well
                popsizeA, popsizeB = len(i[0])/2, len(i[1])/2
                fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
                jsfs = allel.joint_sfs(ac1[:, 1], ac2[:, 1])
                fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
            else:
                # pad for allel as well
                popsizeA, popsizeB = len(i[0]), len(i[1])
                fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
                jsfs = allel.stats.joint_sfs_folded(ac1, ac2)
                fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
            jsfsarray = np.zeros(23)
            jsfsarray[0] = np.sum(fs[0, 1:3])
            jsfsarray[1] = np.sum(fs[1:3, 0])
            jsfsarray[2] = np.sum(fs[0, 3:-3])
            jsfsarray[3] = np.sum(fs[3:-3, 0])
            jsfsarray[4] = np.sum(fs[0, -3:-1])
            jsfsarray[5] = np.sum(fs[-3:-1, 0])
            jsfsarray[6] = np.sum(fs[1:3, 1:3])
            jsfsarray[7] = np.sum(fs[1:3, 3:-3])
            jsfsarray[8] = np.sum(fs[3:-3, 1:3])
            jsfsarray[9] = np.sum(fs[-3:-1, 3:-3])
            jsfsarray[10] = np.sum(fs[3:-3, -3:-1])
            jsfsarray[11] = np.sum(fs[1:3, -3:-1])
            jsfsarray[12] = np.sum(fs[-3:-1, 1:3])
            jsfsarray[13] = np.sum(fs[3:-3, 3:-3])
            jsfsarray[14] = np.sum(fs[-3:-1, -3:-1])
            jsfsarray[15] = np.sum(fs[0, -1])
            jsfsarray[16] = np.sum(fs[-1, 0])
            jsfsarray[17] = np.sum(fs[-1, 1:3])
            jsfsarray[18] = np.sum(fs[1:3, -1])
            jsfsarray[19] = np.sum(fs[-1, 3:-3])
            jsfsarray[20] = np.sum(fs[3:-3, -1])
            jsfsarray[21] = np.sum(fs[-1, -3:-1])
            jsfsarray[22] = np.sum(fs[-3:-1, -1])
            jsfslist.append(jsfsarray)
        return(jsfslist)

    def treeStats(self, treelist, npop):
        """Tree stats with msprime
        """
        for ts in treelist:
            x = msp.GeneralStatCalculator(ts)
            y = msp.BranchLengthStatCalculator(ts)
            z = msp.SiteStatCalculator(ts)

    def filetStats(self, treelist, npop):
        """use filet for stats array
        """




















