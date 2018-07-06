#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:05:14 2018
Class for calculating summary statistics from simulations
@author: stsmall
"""
from __future__ import print_function
from __future__ import division
from subprocess import run, PIPE
import numpy as np
import allel
from itertools import combinations


class SimStats:
    """Calculates stats from pos/gt and msprime object
    """
    def __init__(self, gtlist=None, pops=None, pos=None):
        self.gtlist = gtlist
        self.pops = pops
        self.pos = pos

    def jsfsStats(self, fold=False):
        """Joint site frequency spectrum with scikit-allel
        """
        print("jsfs")
        gtT = [g.transpose() for g in self.gtlist]
        gt = allel.HaplotypeArray(np.vstack(gtT))
        n = 100000  # number of SNPs to choose randomly
        try:
            vidx = np.random.choice(gt.shape[0], n, replace=False)
        except ValueError:
            vidx = np.random.choice(gt.shape[0], gt.shape[0], replace=False)
        vidx.sort()
        gtr = gt.take(vidx, axis=0)
        jsfslist = []
        for i, j in combinations(self.pops, 2):
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

    def asfsStats(self, fold=False):
        """Aggregate SFS, singletons and doubletons
        """
        print("asfs")
        gtT = [g.transpose() for g in self.gtlist]
        gt = allel.HaplotypeArray(np.vstack(gtT))
        n = 100000  # number of SNPs to choose randomly
        try:
            vidx = np.random.choice(gt.shape[0], n, replace=False)
        except ValueError:
            vidx = np.random.choice(gt.shape[0], gt.shape[0], replace=False)
        vidx.sort()
        gtr = gt.take(vidx, axis=0)
        aSFS1 = []
        aSFS2 = []
        for p in self.pops:
            gtp = gtr.take(p, axis=1)
            sfsp = (allel.sfs(gtp.count_alleles()[:, 1]))
            tots = np.sum(sfsp)
            aSFS1.append(sfsp[1]/tots)
            aSFS2.append(sfsp[2]/tots)
        return(aSFS1, aSFS2)

    def filetStats(self, block, filetpath):
        """use filet for stats array
        """
        print("filet")
        filet = []
        loci = len(self.gtlist)
        norm = np.array([block, block**2, block, block, block, 1, block, 1, 1,
                         block, block**2, block, block, block, 1, block, 1, 1,
                         1, 1, block, block, block, 1, 1, 1, 1, 1, 1, 1, 1])
        for pop1, pop2 in combinations(self.pops, 2):
            n1 = len(pop1)
            n2 = len(pop2)
            fakems = []
            fakems.append("ms {} {} -t tbs -r tbs {} -I 2 {} {}\n1234\n".format(n1+n2, loci, block, n1, n2))
            for i, g in enumerate(self.gtlist):
                gt = allel.HaplotypeArray(g.transpose())
                gt12 = gt.take(pop1+pop2, axis=1)
                ac = gt12.count_alleles()
                seg_pos = ac.is_segregating()
                seg = np.count_nonzero(seg_pos)
                gt_seg = gt.compress(seg_pos, axis=0)
                posit = self.pos[i][seg_pos]
                fakems.append("\n//\nsegsites: {}\npositions: {}\n".format(seg, " ".join(map(str, posit))))
                # this needs to be retransposed using g12_seg and pop1 pop2
                gt1 = gt_seg.take(pop1, axis=1)
                for hap in gt1.transpose():
                    fakems.append("{}\n".format("".join(map(str, hap))))
                gt2 = gt_seg.take(pop2, axis=1)
                for hap in gt2.transpose():
                    fakems.append("{}\n".format("".join(map(str, hap))))
            msinput = "".join(fakems)
            cmd = ["{}/twoPopnStats_forML".format(filetpath), str(n1), str(n2)]
            proc = run(cmd, stdout=PIPE, input=msinput, encoding='ascii')
            lstats = proc.stdout.rstrip().split('\n')[1:]
            filetlist = [list(map(float, l.split())) for l in lstats]
            filetmean = np.mean(np.vstack(filetlist), axis=0)
            filetnorm = filetmean / norm
            filet.append(filetnorm)
        return(filet)
