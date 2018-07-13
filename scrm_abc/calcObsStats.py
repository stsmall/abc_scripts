#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 13:43:26 2018

@author: scott
"""

from __future__ import division
from __future__ import print_function
import allel
import numpy as np
import pandas as pd
from allel_class import Chr
import autil as autil
from itertools import combinations
import argparse

import apca as apca

parser = argparse.ArgumentParser()
parser.add_arguement('-v', "--vcfFile", help="path to vcf")
parser.add_argument('--h5', action="store_true", help="h5 exists")
parser.add_arguement('-m', "--meta", required=True, help="path to meta data")
args = parser.parse_args()


def makeh5fromvcf(vcfin, altnum, hf5):
    """
    """
    h5out = "{}.h5".format(vcfin.split(".")[:-2])
    if hf5:
        pass
    else:
        fieldsfromvcf = ['samples', 'calldata/GQ', 'variants/ALT',
                         'variants/REF', 'variants/QUAL', 'variants/CHROM',
                         'variants/POS', 'variants/AF', 'variants/AB',
                         'variants/MQM', 'variants/DP', 'calldata/DP',
                         'calldata/AD', 'calldata/GT']
        allel.vcf_to_hdf5(vcfin, h5out, fields=fieldsfromvcf,
                          types={'calldata/GQ': 'float32'}, alt_number=2)
    # callset = h5py.File(h5out, mode='r')
    return(None)


#    def jsfsStats(gt, pops, fold=False):
#        """Joint site frequency spectrum with scikit-allel
#        """
#        print("jsfs")
#        gtT = [g.transpose() for g in self.gtlist]
#        gt = allel.HaplotypeArray(np.vstack(gtT))
#        jsfslist = []
#        for i, j in combinations(pops, 2):
#            gtpops = gt.take(i+j, axis=1)
#            acpops = gtpops.count_alleles()
#            seg = acpops.is_segregating()
#            gtseg = gt.compress(seg)
#            # random snps
#            n = 100000  # number of SNPs to choose randomly
#            try:
#                vidx = np.random.choice(gtseg.shape[0], n, replace=False)
#            except ValueError:
#                vidx = np.random.choice(gtseg.shape[0], gtseg.shape[0], replace=False)
#            vidx.sort()
#            gtr = gtseg.take(vidx, axis=0)
#            gtpop1 = gtr.take(i, axis=1)
#            gtpop2 = gtr.take(j, axis=1)
#            ac1 = gtpop1.count_alleles()
#            ac2 = gtpop2.count_alleles()
#            if fold:
#                # pad for allel as well
#                popsizeA, popsizeB = len(i)/2, len(j)/2
#                fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
#                jsfs = allel.joint_sfs_folded(ac1, ac2)
#                fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
#            else:
#                # pad for allel as well
#                popsizeA, popsizeB = len(i), len(j)
#                fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
#                jsfs = allel.joint_sfs(ac1[:, 1], ac2[:, 1])
#                fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
#            jsfsarray = np.zeros(23)
#            jsfsarray[0] = np.sum(fs[0, 1:3])
#            jsfsarray[1] = np.sum(fs[1:3, 0])
#            jsfsarray[2] = np.sum(fs[0, 3:-3])
#            jsfsarray[3] = np.sum(fs[3:-3, 0])
#            jsfsarray[4] = np.sum(fs[0, -3:-1])
#            jsfsarray[5] = np.sum(fs[-3:-1, 0])
#            jsfsarray[6] = np.sum(fs[1:3, 1:3])
#            jsfsarray[7] = np.sum(fs[1:3, 3:-3])
#            jsfsarray[8] = np.sum(fs[3:-3, 1:3])
#            jsfsarray[9] = np.sum(fs[-3:-1, 3:-3])
#            jsfsarray[10] = np.sum(fs[3:-3, -3:-1])
#            jsfsarray[11] = np.sum(fs[1:3, -3:-1])
#            jsfsarray[12] = np.sum(fs[-3:-1, 1:3])
#            jsfsarray[13] = np.sum(fs[3:-3, 3:-3])
#            jsfsarray[14] = np.sum(fs[-3:-1, -3:-1])
#            jsfsarray[15] = np.sum(fs[0, -1])
#            jsfsarray[16] = np.sum(fs[-1, 0])
#            jsfsarray[17] = np.sum(fs[-1, 1:3])
#            jsfsarray[18] = np.sum(fs[1:3, -1])
#            jsfsarray[19] = np.sum(fs[-1, 3:-3])
#            jsfsarray[20] = np.sum(fs[3:-3, -1])
#            jsfsarray[21] = np.sum(fs[-1, -3:-1])
#            jsfsarray[22] = np.sum(fs[-3:-1, -1])
#            jsfslist.append(jsfsarray)
#        return(jsfslist)
#


def jsfsStats(gt, pops, fold=False):
    """Joint site frequency spectrum with scikit-allel
    """
    print("jsfs")
    n = 100000  # number of SNPs to choose randomly
    try:
        vidx = np.random.choice(gt.shape[0], n, replace=False)
    except ValueError:
        vidx = np.random.choice(gt.shape[0], gt.shape[0], replace=False)
    vidx.sort()
    gtr = gt.take(vidx, axis=0)
    jsfslist = []
    for i, j in combinations(pops, 2):
        gtpop1 = gtr.take(i, axis=1)
        gtpop2 = gtr.take(j, axis=1)
        ac1 = gtpop1.count_alleles()
        ac2 = gtpop2.count_alleles()
        if fold:
            # pad for allel as well
            popsizeA, popsizeB = len(i)/2, len(j)/2
            fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
            jsfs = allel.joint_sfs_folded(ac1, ac2)
            fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
        else:
            # pad for allel as well
            popsizeA, popsizeB = len(i), len(j)
            fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
            jsfs = allel.joint_sfs(ac1[:, 1], ac2[:, 1])
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


def asfsStats(gt, pops, fold=False):
    """Aggregate SFS, singletons and doubletons
    """
    print("asfs")
    n = 100000  # number of SNPs to choose randomly
    try:
        vidx = np.random.choice(gt.shape[0], n, replace=False)
    except ValueError:
        vidx = np.random.choice(gt.shape[0], gt.shape[0], replace=False)
    vidx.sort()
    gtr = gt.take(vidx, axis=0)
    aSFS1 = []
    aSFS2 = []
    for p in pops:
        gtp = gtr.take(p, axis=1)
        sfsp = (allel.sfs(gtp.count_alleles()[:, 1]))
        tots = np.sum(sfsp)
        aSFS1.append(sfsp[1]/tots)
        aSFS2.append(sfsp[2]/tots)
    return(aSFS1, aSFS2)


#    def asfsStats(gt, pops, fold=False):
#        """Aggregate SFS, singletons and doubletons
#        """
#        print("asfs")
#        aSFS1 = []
#        aSFS2 = []
#        for p in pops:
#            gtpop = gt.take(p, axis=1)
#            acpop = gtpop.count_alleles()
#            seg = acpop.is_segregating()
#            gtseg = gtpop.compress(seg)
#            # random snps
#            n = 100000  # number of SNPs to choose randomly
#            try:
#                vidx = np.random.choice(gtseg.shape[0], n, replace=False)
#            except ValueError:
#                vidx = np.random.choice(gtseg.shape[0], gtseg.shape[0], replace=False)
#            vidx.sort()
#            gtp = gtseg.take(vidx, axis=0)
#            sfsp = (allel.sfs(gtp.count_alleles()[:, 1]))
#            tots = np.sum(sfsp)
#            aSFS1.append(sfsp[1]/tots)
#            aSFS2.append(sfsp[2]/tots)
#        return(aSFS1, aSFS2)


if __name__ == "__main__":
    makeh5fromvcf(args.vcfFile, 1)
    meta = args.meta
    meta = pd.read_csv(meta, delimiter=",")
    var = Chr('All', "{}.h5".format(args.vcfFile))
    popdict = autil.subpops(var, meta, bypop=True, bykary=False)
    pop2color = autil.popcols(popdict)
    chrlist = np.unique(var.chrm[:])
    pops = popdict.keys()
    for c in chrlist:
        var.geno(c, meta)
        asfsStats(var.gt, pops)
        jsfsStats(var.gt, pops)


    # PCA and LD thin
    thinpos = {}
    pcadict = {}
    gnudict = {}
    for c in chrlist:
        # LD thin
        print(c)
        var.geno(c, meta)
        var.miss(var.gt, var.pos, 0)  # use only sites without missing data
        gn, thinp = autil.ldthin(var.gt, var.pos, "thin", iters=5)
        gnudict[c] = gn
        thinpos[c] = thinp
    for c in gnudict.keys():
        # PCA
        gn = gnudict[c]
        coords, model = apca.pca_fx(gn, meta, c, pop2color, False, var.pop,
                                    bykary=True)
        pcadict[c] = (coords, model)
