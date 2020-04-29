#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 18:05:14 2018
Class for calculating summary statistics from simulations
@author: stsmall
"""
from subprocess import run, PIPE
import numpy as np
import allel
from itertools import combinations
import bisect
import dadi


class SumStats:
    """Calculates stats from pos/gt and msprime object."""

    def __init__(self, gtlist=None, pops=None, pos=None):
        self.gtlist = gtlist
        self.pops = pops
        self.pos = pos

    def treeseq2gtarray(self, tree_sequence):
        """Calculate stats from msprime tree sequence."""
        V = np.zeros((tree_sequence.get_num_mutations(),
                      tree_sequence.get_sample_size()), dtype=np.int8)
        for variant in tree_sequence.variants():
            V[variant.index] = variant.genotypes
        gt = allel.HaplotypeArray(V)
        pos = allel.SortedIndex([int(v.position) for v in tree_sequence.variants()])
        return(pos, gt)

    def asfsStats(self, rand=True, fold=False, randn=100000):
        """Calculate the aggregate SFS, singletons and doubletons.

        Parameters
        ----------
        rand : TYPE, optional
            DESCRIPTION. The default is True.
        fold : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        aSFS1 : TYPE
            DESCRIPTION.
        aSFS2 : TYPE
            DESCRIPTION.
        aSFS : TYPE
            DESCRIPTION.

        """
        gt = self.gtlist
        aSFS1 = []
        aSFS2 = []
        aSFS = []
        for p in self.pops:
            gtpop = gt.take(p, axis=1)
            acpop = gtpop.count_alleles()
            seg = acpop.is_segregating()
            gtseg = gtpop.compress(seg)
            # random snps
            if rand:
                n = randn  # number of SNPs to choose randomly
                try:
                    vidx = np.random.choice(gtseg.shape[0], n, replace=False)
                except ValueError:
                    vidx = np.random.choice(gtseg.shape[0], gtseg.shape[0], replace=False)
            else:
                vidx = np.random.choice(gtseg.shape[0], gtseg.shape[0], replace=False)
            vidx.sort()
            gtp = gtseg.take(vidx, axis=0)
            # sfs
            if fold:
                # TODO: check this
                sfsp = (allel.sfs_fold(gtp.count_alleles()[:, 1]))
            else:
                sfsp = (allel.sfs(gtp.count_alleles()[:, 1]))
            tots = np.sum(sfsp)
            aSFS1.append(sfsp[1]/tots)
            aSFS2.append(sfsp[2]/tots)
            aSFS.append(sfsp)
        return aSFS1, aSFS2, aSFS

    def sfs_dadi(self, npop, seg=True, fold=False, jsfs=True):
        """Calculate the site frequency spectrum with dadi.

        Parameters
        ----------
        npop : TYPE
            DESCRIPTION.
        seg : TYPE, optional
            DESCRIPTION. The default is True.
        fold : TYPE, optional
            DESCRIPTION. The default is False.
        jsfs : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        sfs_list :
            asdf
        jsfs_list :
            asdf

        """
        gt = self.gtlist
        sfs_list = []
        jsfs_list = []
        for i, p in enumerate(self.pops):
            gtpop = gt.take(p, axis=1)
            acpop = gtpop.count_alleles()
            seg = acpop.is_segregating()
            gtseg = gtpop.compress(seg)
            ac = gtseg.count_alleles()
            sfs = allel.sfs(ac[:, 1])
            sfs.resize(len(p), refcheck=False)
            fs = dadi.Spectrum(sfs, mask_corners=False, pop_ids=["pop{}".format(i)])
            if fold:
                fs = fs.fold()
            sfs_list.append(list(fs))
            jsfs.append(fs)
        if jsfs:
            jsfs_list = self.jsfs_dadi(self, jsfs)
            return sfs_list, jsfs_list
        else:
            return sfs_list

    def jsfsStats(self, fold=False, seg=True, rand=True,randn=100000):
        """Calculate joint site frequency spectrum (jsfs) with scikit-allel.

        Parameters
        ----------
        fold : TYPE, optional
            DESCRIPTION. The default is False.
        seg : TYPE, optional
            DESCRIPTION. The default is True.
        rand : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        jsfs_list : TYPE
            DESCRIPTION.

        """
        gt = self.gtlist
        jsfs_list = []
        for i, j in combinations(self.pops, 2):
            if seg:
                gtpops = gt.take(i+j, axis=1)
                acpops = gtpops.count_alleles()
                segpops = acpops.is_segregating()
                gt = gt.compress(segpops)
            # random snps
            if rand:
                n = randn  # number of SNPs to choose randomly
                try:
                    vidx = np.random.choice(gt.shape[0], n, replace=False)
                except ValueError:
                    vidx = np.random.choice(gt.shape[0], gt.shape[0], replace=False)
            else:
                vidx = np.random.choice(gt.shape[0], gt.shape[0], replace=False)
            vidx.sort()
            gtr = gt.take(vidx, axis=0)
            gtpop1 = gtr.take(i, axis=1)
            gtpop2 = gtr.take(j, axis=1)
            ac1 = gtpop1.count_alleles()
            ac2 = gtpop2.count_alleles()
            # jsfs
            if fold:
                # pad for allel as well
                popsizeA, popsizeB = len(i)/2, len(j)/2
                fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
                jsfs = allel.joint_sfs_folded(ac1, ac2)
                fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
            else:
                # pad for allel as well
                popsizeA, popsizeB = len(i)*2, len(j)*2
                fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
                jsfs = allel.joint_sfs(ac1[:, 1], ac2[:, 1])
                fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
            # summarize
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
            jsfs_list.append(jsfsarray)
        return jsfs_list

    def jsfs_dadi(self, jsfs, fold=True):
        """Calculate the joint site freq spectrum with dadi.

        Parameters
        ----------
        jsfs : TYPE
            DESCRIPTION.
        fold : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        jsfs_list

        """
        jsfs_list = []
        for i, j in combinations(jsfs, 2):
            fs = dadi.Spectrum([i, j])
            if fold:
                fs = fs.fold()
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
            jsfs_list.append(jsfsarray)
        return jsfs_list

    def afibs_stats(self, fold=True):
        """Allele Frequency spectrum of ibs lengths."""
        afibs_dict = {}
        gt = self.gtlist
        for i, p in enumerate(self.pops):
            ibsarray = np.zeros((len(p), len(p)-1))
            gtpop = gt.take(p, axis=1)
            acpop = gtpop.count_alleles()
            seg = acpop.is_segregating()
            gtseg = gtpop.compress(seg)
            # ac = gtseg.count_alleles()
            poslist = self.pos[seg]
            freqlist = np.sum(gtseg, axis=1)
            for ind in range(len(p)):
                indpos = poslist[gtseg[:, ind] > 0]
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
            afibs_dict[i] = np.mean(ibsarray, axis=0)
        # fold and export to list
        afibs_list = []
        if fold:
            for i, p in enumerate(self.pops):
                ibs_flip = np.flip(afibs_dict[i], axis=0)
                ibs_fold = (ibs_flip + afibs_dict[i]) / 2
                haps = len(p)/2
                ibs = ibs_fold[0:haps]
                afibs_list.append(ibs)
        else:
            for i in range(len(self.pops)):
                afibs_list.append(afibs_dict[i])
        return afibs_list

    def filetStats(self, block, mask, filetpath):
        """Calculate stats using FILET."""
        filet = []
        loci = len(self.gtlist)
        norm = np.array([block, block**2, block, block, block, 1, block, 1, 1,
                         block, block**2, block, block, block, 1, block, 1, 1,
                         1, 1, block, block, block, 1, 1, 1, 1, 1, 1, 1, 1])
        for pop1, pop2 in combinations(self.pops, 2):
            n1 = len(pop1)
            n2 = len(pop2)
            fakems = []
            fakems.append(f"ms {n1+n2} {loci} -t tbs -r tbs {block} -I 2 {n1} {n2}\n1234\n")
            for i, g in enumerate(self.gtlist):
                gt = g[pop1+pop2]
                seg_pos = np.sum(gt, axis=0)
                seg_mask = (seg_pos > 0) & (seg_pos < (n1+n2))
                seg = np.count_nonzero(seg_mask)
                posit = self.pos[i][seg_mask]
                gt_seg = gt[:, seg_mask]
                fakems.append(f"\n//\nsegsites: {seg}\npositions: {' '.join(map(str, posit))}\n")
                for a in gt_seg:
                    fakems.append(f"{''.join(map(str, a))}\n")
            msinput = "".join(fakems)

            # TODO: fix with mask
            cmd = [f"{filetpath}msMaskAllRows {mask}"]
            proc = run(cmd, stdout=PIPE, input=msinput, encoding='ascii')
            # {filetpath}removeNedOutColumnsFromMsFileKeepStars.py stdin"]
            # TODO: fix with masked frac
            cmd = [f"{filetpath}twoPopnStats_forML {n1} {n2}"]
            proc = run(cmd, stdout=PIPE, input=msinput, encoding='ascii')
            # normalizeTwoPopnStats.py {mask}-unmaskedFrac $window

            lstats = proc.stdout.rstrip().split('\n')[1:]
            filetlist = [list(map(float, l.split())) for l in lstats]
            filetnan = np.vstack(filetlist)
            filetnan[np.isinf(filetnan)] = 'nan'
            filetmean = np.nanmean(filetnan, axis=0)
            filetnorm = filetmean / norm
            filet.append(filetnorm)
        return filet

    # def filetStatsMP(self, args):
    #     """Calculate stats using FILET w/ multiprocessors."""
    #     pop1, pop2, block, mask, filet = args
    #     norm = np.array([block, block**2, block, block, block, 1, block, 1, 1,
    #                      block, block**2, block, block, block, 1, block, 1, 1,
    #                      1, 1, block, block, block, 1, 1, 1, 1, 1, 1, 1, 1])
    #     loci = len(self.gtlist)
    #     n1 = len(pop1)
    #     n2 = len(pop2)
    #     fakems = []
    #     fakems.append(f"ms {n1+n2} {loci} -t tbs -r tbs {block} -I 2 {n1} {n2}\n1234\n")
    #     for i, g in enumerate(self.gtlist):
    #         gt = g[pop1+pop2]
    #         seg_pos = np.sum(gt, axis=0)
    #         seg_mask = (seg_pos > 0) & (seg_pos < (n1+n2))
    #         seg = np.count_nonzero(seg_mask)
    #         posit = self.pos[i][seg_mask]
    #         gt_seg = gt[:, seg_mask]
    #         fakems.append(f"\n//\nsegsites: {seg}\npositions: {' '.join(map(str, posit))}\n")
    #         for a in gt_seg:
    #             fakems.append(f"{''.join(map(str, a))}\n")
    #     msinput = "".join(fakems)
    #     cmd = [f"{filet}twoPopnStats_forML {n1} {n2}"]
    #     proc = run(cmd, stdout=PIPE, input=msinput, encoding='ascii')
    #     lstats = proc.stdout.rstrip().split('\n')[1:]
    #     filetlist = [list(map(float, l.split())) for l in lstats]
    #     filetnan = np.vstack(filetlist)
    #     filetnan[np.isinf(filetnan)] = 'nan'
    #     filetmean = np.nanmean(filetnan, axis=0)
    #     filetnorm = filetmean / norm
    #     return filetnorm
