#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 13:43:26 2018

@author: scott

SFS/JSFS/AFIBS OBSERVED
--------------
vcf to allel ~ examples from diploshic
+ gff to find sites far from genes
+ mask?
+ mispolarization
+ phased for afibs
+ can make masking files directly
https://github.com/kern-lab/diploSHIC/blob/master/makeFeatureVecsForChrArmFromVcfDiploid.py
https://github.com/kern-lab/diploSHIC/blob/master/makeFeatureVecsForChrArmFromVcf_ogSHIC.py


FILET OBSERVED
-----
https://github.com/kern-lab/FILET/tree/master/maskingScripts
1) starting with: filtered and phased vcf
    maskedRef.fasta  # mask a fasta for repeat, missing data etc ...
    phasedVcfsToFastas.py vcf maskedRef.fasta fastaout.1 fastaout.2
    pgStatsBedSubpop_forML pop1.fa pop2.fa anc.fa coords.bed 0.5 | python normalizePgStats.py windowSize > out.fv
2) starting with: None
    makeMaskArmFilesFromQualAndRM.py
    makeFilteredSnpVcfsForEachArm.py
    PHASING -- try shapeit
    phasedVcfsToFastas.py
    pgStatsBedSubpop_forML pop1.fa pop2.fa anc.fa coords.bed 0.5 | python normalizePgStats.py windowSize > out.fv
"""

import argparse
import allel
import numpy as np
import pandas as pd
from itertools import combinations

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')

from allel_class import Chr
import autil as autil

from sim_stats import *

parser = argparse.ArgumentParser()
parser.add_argument('-v', "--vcfFile", help="path to vcf")
parser.add_argument('--h5', action="store_true", help="h5 exists")
parser.add_argument('-m', "--meta", required=True, help="path to meta data")
args = parser.parse_args()


# TO DO: use discoal functions to clean this up. additional input of gff to avoid genes
def build_mask(pairs, sample_names, vcf, gvcf, repeat_gff, ref_fa, prcnt_miss, qual):
    """Make mask file.

    Parameters
    ----------
    pairs : TYPE
        DESCRIPTION.
    sample_names : TYPE
        DESCRIPTION.
    gvcf : TYPE
        DESCRIPTION.
    repeat_gff : TYPE
        DESCRIPTION.
    ref_fa : TYPE
        DESCRIPTION.
    prcnt_miss : TYPE
        DESCRIPTION.
    qual_thresh : TYPE
        DESCRIPTION.

    Returns
    -------
    maskout : TYPE
        DESCRIPTION.

    """
    mask_files = []
    vcf_files = []
    for p1, p2 in combinations(pairs):
        # get these gvcfs or split from larger using sample names
        mask_out = f"{p1}-{p2}.mask.fa"
        import makeMaskArmsFilesFromQualAndRM as make_mask
        make_mask(gvcf1, gvcf2, reapeat_gff, ref_fa, qual, prcnt_miss, mask_out)
        mask_files.append(mask_out)
        import makeFilteredSnpVcfsFroEachArm as filt_vcf
        filt_vcf(vcf1, vcf2, gvcf1, gvcf2, mask_out, qual, vcf_out)
        vcf_files.append(mask_out)
    return mask_files, vcf_files


def vcf_to_fasta(vcf_path, vcf_files, mask_path, mask_files):
    """Create fasta from masked and phased VCF.

    Parameters
    ----------
    vcf : TYPE
        DESCRIPTION.
    mask_out : TYPE
        DESCRIPTION.
    fasta_out1 : TYPE
        DESCRIPTION.
    fasta_out2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None

    """
    # vcf_list = ['Fun_Par_3R.vcf', 'Fun_Van_3R.vcf']
    for vcf in vcf_list:
        import phasedVcfsToFastas as vcf2fasta
        pop1, pop2, arm = vcf.strip("vcf").split("_")
        mask = f"{pop1}-{pop2}.mask.fa"
        assert mask in mask_list
        vcf2fasta(PATHvcf, PATHmask, f"{pop1}.{pop1}-{pop2}.{arm}.masked.fasta", f"{pop2}.{pop1}-{pop2}.{arm}.masked.fasta")
        fasta_files.append(f"{pop1}-{pop2}-{arm}")
    return fasta_files


def make_coords(chrom, chrom_len, start, stop, size, step):
    """Make windows for FILET.

    Parameters
    ----------
    chrom : TYPE
        DESCRIPTION.
    chrom_len : TYPE
        DESCRIPTION.
    start : TYPE
        DESCRIPTION.
    stop : TYPE
        DESCRIPTION.
    size : TYPE
        DESCRIPTION.
    step : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if start:
        window = f"{chrom}.{size}-{step}.{start}-{stop}.bed"
    else:
        window = f"{chrom}.{size}-{step}.bed"


def makeh5fromvcf(vcfin, altnum, hf5):
    """
    """
    h5out = "{}.h5".format(vcfin)
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


def asfsStats(gt, pops, chrm, rand=True, plot=False):
    """Aggregate SFS, singletons and doubletons
    """
    print("asfs")
    aSFS1 = []
    aSFS2 = []
    for p in pops:
        gtpop = gt.take(p, axis=1)
        acpop = gtpop.count_alleles()
        seg = acpop.is_segregating()
        gtseg = gtpop.compress(seg)
        # random snps
        if rand:
            n = 100000  # number of SNPs to choose randomly
            try:
                vidx = np.random.choice(gtseg.shape[0], n, replace=False)
            except ValueError:
                vidx = np.random.choice(gtseg.shape[0], gtseg.shape[0], replace=False)
        else:
            vidx = np.random.choice(gtseg.shape[0], gtseg.shape[0], replace=False)
        vidx.sort()
        gtp = gtseg.take(vidx, axis=0)
        sfsp = (allel.sfs(gtp.count_alleles()[:, 1]))
        print(sfsp)
        if plot:
            fig, ax = plt.subplots(figsize=(6, 6))
            allel.stats.plot_sfs(sfsp, ax=ax)
        tots = np.sum(sfsp)
        aSFS1.append(sfsp[1]/tots)
        aSFS2.append(sfsp[2]/tots)
    return(aSFS1, aSFS2)


def jsfsStats(self, fold=False, seg=True, rand=True, plot=True):
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
   jsfs_list : None

    """
    # gtT = [g.transpose() for g in self.gtlist]
    # gt = allel.HaplotypeArray(np.vstack(gtT))
    gt = 'x'
    jsfs_list = []
    for i, j in combinations(self.pops, 2):
        if seg:
            gtpops = gt.take(i+j, axis=1)
            acpops = gtpops.count_alleles()
            segpops = acpops.is_segregating()
            gt = gt.compress(segpops)
        # random snps
        if rand:
            n = 100000  # number of SNPs to choose randomly
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
        if plot:
            fig, ax = plt.subplots(figsize=(6, 6))
            allel.stats.plot_joint_sfs(fs, ax=ax)
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


if __name__ == "__main__":
    makeh5fromvcf(args.vcfFile, 1)
    meta = args.meta
    meta = pd.read_csv(meta, delimiter=",")
    var = Chr('All', "{}.h5".format(args.vcfFile))
    popdict = autil.subpops(var, meta, bypop=True, bykary=False)
    pop2color = autil.popcols(popdict)
    chrlist = np.unique(var.chrm[:])
    pops = list(popdict.values())
    sfsdict = {}
    jsfsdict = {}
    for c in chrlist:
        var.geno(c, meta)
        #sfsdict[c] = asfsStatsSeg(var.gt, pops, c, rand=False, plot=False)
        sfsdict[c] = asfsStats(var.gt, pops, c, rand=False, plot=False)
        #jsfsdict[c] = jsfsStatsSeg(var.gt, pops, c, fold=False, rand=False, plot=True)
        #jsfsdict[c] = jsfsStats(var.gt, pops, c)

    # asfs
    s1 = []
    s2 = []
    for chrm in sfsdict.keys():
        s1.append(sfsdict[chrm][0])
        s2.append(sfsdict[chrm][1])
    s1array = np.mean(np.vstack(s1), axis=0)
    s2array = np.mean(np.vstack(s2), axis=0)

    # jsfs
    props = []
    for chrm in jsfsdict.keys():
        jsfslist = jsfsdict[chrm]
        jsfstotal = np.sum(jsfslist, axis=1)
        props.append([j/jsfstotal[i] for i, j in enumerate(jsfslist)])
    jsfs = []
    for pairs in range(len(props[0])):
        p = []
        for chrm in props:
            p.append(chrm[pairs])
        jsfs.append(np.mean(np.vstack(p), axis=0))
    # write out
    s1 = " ".join(map(str, list(s1array)))
    s2 = " ".join(map(str, list(s2array)))
    j23 = " ".join(map(str, np.concatenate(jsfs).ravel()))
    f = open("Observed_summStats.out", 'w')
    f.write("{} {} {}\n".format(s1, s2, j23))
    f.close()
