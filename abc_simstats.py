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
            gtpop = gt[:, p]
            for ind in range(len(pix)):
                ibslist = []
                for freq in range(1, len(p)):
                    ibs = 0
                    # snp position only
                    het = np.sum(gtpop, axis=1) > 0  # mask
                    # positions of snps
                    poslist = pos[het]  # genome positions not index
                    # freq of snp positions
                    freqlist = np.sum(gtpop, axis=1)[het]
                    mut_ix = np.where(freqlist == freq)[0]
                    for m in mut_ix:
                        if m == 0:
                            ibs += poslist[m + 1]
                        else:
                            try:
                                ibs += poslist[m + 1] - poslist[m - 1]
                            except IndexError:
                                ibs += length - poslist[m - 1]
                    try:
                        ibslist.append(ibs / len(mut_ix))
                    except ZeroDivisionError:
                        # nothing in that freq class
                        ibslist.append(0)
            afibsdict[i] = ibslist
        # get mean
        ibsstats = []
        if fold:
            for i, p in enumerate(pix):
                ibs_mean = np.mean(afibsdict[i], axis=0)
                ibs_meanr = np.flip(ibs_mean, axis=0)
                ibs_fold = (ibs_mean + ibs_meanr) / 2
                ibs = ibs_fold[0:len(p)]
                ibsstats.append(ibs)
        else:
            for i in range(len(pix)):
                ibs = np.mean(afibsdict[i], axis=0)
                ibsstats.append(ibs)
        return(ibsstats)
