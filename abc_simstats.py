#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:49:19 2017
This python code directly calculates statistics from ms input. Preferrably this
should be used with msprime
ms2stats.py
@author: scott
"""
import numpy as np
import allel
from gwas import hapdata
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

    def ld_stats(self, tree_sequence, interval_list, configdict, hapdict, posdict, countdict, prob_err=0):
        """Linkage Disequilibrium decay using popsizeABC code
        """
        # TODO: set this to work with msprime.LdCalculator()
        # with popsizeabc
        pops = len(configdict["sample_size"])
        haps = [i*2 for i in configdict["sample_size"]]
        muts = tree_sequence.get_num_mutations()
        for pop in range(pops):
            p_ix = tree_sequence.get_samples(pop)
            mydata = hapdata.HapDataset('ms_rep', nsnp=muts, nhaplo=haps[pop])
            pos = []
            for mut in tree_sequence.mutations():
                pos.append(mut[0])
            mydata.snp_pos = np.array(pos, dtype='float')
            # read haplotypes
            j = 0
            pophaps = np.array(list(tree_sequence.haplotypes()))[p_ix]
            for h in pophaps:  # 1 pop at a time
                mydata.UpdatePop(j, 'pop1')
                mydata.Data[j, :] = np.array(list(h)[:muts], dtype='int16')
                j += 1
            # updates lists
            if prob_err > 0:
                mydata.introduce_errors(prob_err)
            u = hapdata.Haplos_and_counts_fast(mydata, 'pop1', mac=configdict["mac"])
            countdict[pop].append(u[1])
            u = hapdata.Haplos_and_counts_fast(mydata, 'pop1', mac=configdict["mac_ld"])
            posdict[pop].append(u[0])
            hapdict[pop].append(u[1])
        # with msprime
#        ld_calc = msp.LdCalculator(tree_sequence)
#        A = ld_calc.get_r2_array()
#        distance between tree_sequence.mutations()
        return(hapdict, posdict, countdict)

    def tree2gtarray(self, tree_sequence):
        """
        """
        V = np.zeros((tree_sequence.get_num_mutations(),
                      tree_sequence.get_sample_size()), dtype=np.int8)
        for variant in tree_sequence.variants():
            V[variant.index] = variant.genotypes
        gt = allel.HaplotypeArray(V)
        pos = allel.SortedIndex([int(v.position) for v in tree_sequence.variants()])
        return(pos, gt)

    def sfs_stats(self, tree_sequence, npop):
        """
        """
        pos, gt = self.tree2gtarray(tree_sequence)
        pix = []
        sfslist = []
        pix = [tree_sequence.get_samples(pop) for pop in range(npop)]
        for p in pix:
            gtpop = gt[:, p]
            ac = gtpop.count_alleles()
            # TODO: this needs to be padded so the length is consistent; check the msp2dadi.py script
            sfslist.append(allel.sfs_folded(ac))
        return(sfslist)

    def jsfs_stats(self, tree_sequence, npop):
        """Joint site frequency spectrum
        """
        pos, gt = self.tree2gtarray(tree_sequence)
        pix = []
        jsfslist = []
        pix = [tree_sequence.get_samples(pop) for pop in range(npop)]
        for i, j in combinations(pix, 2):
            gtpop1 = gt[:, i]
            gtpop2 = gt[:, j]
            ac1 = gtpop1.count_alleles()
            ac2 = gtpop2.count_alleles()
            fs = allel.stats.joint_sfs_folded(ac1, ac2)
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

    def afibs_stats(self, tree_sequence, length, npop, fold=True):
        """Allele Frequency spectrum of ibs lengths
        """
        pos, gt = self.tree2gtarray(tree_sequence)
        pix = [tree_sequence.get_samples(pop) for pop in range(npop)]
        afibsdict = {}
        for i, p in enumerate(pix):  # single pop
            ibsarray = np.zeros((len(p), len(p)-1))
            gtpop = gt[:, p]
            # segregating positions only
            het = np.sum(gtpop, axis=1) > 0  # mask
            # subsample for seg only
            poslist = pos[het]
            gtpop_seg = gtpop[het]
            freqlist = np.sum(gtpop_seg, axis=1)
            for ind in range(len(p)):
                indpos = poslist[gtpop_seg[:, ind] > 0]
                indpos = np.insert(indpos, 0, 0)
                indpos = np.insert(indpos, len(indpos), length)
                for freq in range(1, len(p)):
                    ibs = 0
                    mut_ix = np.where(freqlist == freq)[0]
                    for m in mut_ix:
                        start = bisect.bisect_left(indpos, poslist[m])
                        end = bisect.bisect_right(indpos, poslist[m])
                        ibs += indpos[end] - indpos[start - 1]
                    try:
                        ibsarray[ind, freq-1] = (ibs / len(mut_ix))
                    except ZeroDivisionError:
                        # nothing in that freq class
                        ibsarray[ind, freq-1] = 0
            afibsdict[i] = np.mean(ibsarray, axis=0)
        # fold and export to list
        ibsstats = []
        if fold:
            for i, p in enumerate(pix):
                ibs_flip = np.flip(afibsdict[i], axis=0)
                ibs_fold = (ibs_flip + afibsdict[i]) / 2
                haps = len(p)/2
                ibs = ibs_fold[0:haps]
                ibsstats.append(ibs)
        else:
            for i in range(len(pix)):
                ibsstats.append(afibsdict[i])
        return(ibsstats)
