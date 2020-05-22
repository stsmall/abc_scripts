#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: stsmall
"""

import allel
import numpy as np
import gzip
import bisect


def readSampleToPopFile(sampleToPopFileName):
    table = {}
    with open(sampleToPopFileName) as sampleToPopFile:
        for line in sampleToPopFile:
            sample, pop = line.strip().split()
            table[sample] = pop
    return table


def readMaskDataForScan(maskFileName, chrArm):
    isAccessible = []
    readingMasks = False
    if maskFileName.endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open
    with fopen(maskFileName, 'rt') as maskFile:
        for line in maskFile:
            if line.startswith(">"):
                currChr = line[1:].strip()
                if currChr == chrArm:
                    readingMasks = True
                elif readingMasks:
                    break
            else:
                if readingMasks:
                    for char in line.strip().upper():
                        if char == 'N':
                            isAccessible.append(False)
                        else:
                            isAccessible.append(True)
    return isAccessible


def asfsStats(args, fold=False, rand=True, randn=100000):
    """

    Parameters
    ----------
    pos : TYPE
        DESCRIPTION.
    hap : TYPE
        DESCRIPTION.
    pops : TYPE
        DESCRIPTION.
    fold : TYPE, optional
        DESCRIPTION. The default is False.
    rand : TYPE, optional
        DESCRIPTION. The default is True.
    randn : TYPE, optional
        DESCRIPTION. The default is 100000.

    Returns
    -------
    asfs : TYPE
        DESCRIPTION.

    """
    pos, gt, pops = args
    aSFS12 = []
    aSFS = []
    for pop in pops:
        gtpop = gt.take(pop, axis=1)
        miss_count = gtpop.count_missing(axis=1)
        miss_arr = miss_count > 0
        gtpop = gtpop.compress(miss_arr, axis=0)
        acpop = gtpop.count_alleles()
        seg = acpop.is_segregating()
        try:
            gtseg = gtpop.compress(seg)
            # random snps
            if rand:
                n = randn  # number of SNPs to choose randomly
                try:
                    vidx = np.random.choice(gtseg.shape[0], n, replace=False)
                    vidx.sort()
                    gtp = gtseg.take(vidx, axis=0)
                except ValueError:
                    gtp = gtseg
            else:
                gtp = gtseg
            # sfs
            if fold:
                sfsp = (allel.sfs_folded(gtp.count_alleles()))
            else:
                sfsp = (allel.sfs(gtp.count_alleles()[:, 1]))
            tots = np.sum(sfsp)
        except ValueError:
            # no segregating snps
            tots = 1
            sfsp = [0]*len(pop)
        try:
            aSFS12.append(sfsp[1]/tots)
        except IndexError:
            aSFS12.append(0)
        try:
            aSFS12.append(sfsp[2]/tots)
        except IndexError:
            aSFS12.append(0)
        aSFS.append(sfsp)
    asfs = " ".join(map(str, aSFS12))
    return f"{asfs}\n"


def summarizejsfs(fs):
    """Create summary jsfs.

    Parameters
    ----------
    fs : TYPE
        DESCRIPTION.

    Returns
    -------
    props : TYPE
        DESCRIPTION.

    """
    # # summarize
    jsfsarray = np.zeros(23, dtype=int)
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
    jsfstotal = np.sum(jsfsarray)
    props = jsfsarray/jsfstotal
    return props


def jsfsStats(args, fold=False, rand=True, randn=100000):
    """

    Parameters
    ----------
    pos : TYPE
        DESCRIPTION.
    hap : TYPE
        DESCRIPTION.
    pops : TYPE
        DESCRIPTION.
    pairs : TYPE
        DESCRIPTION.
    fold : TYPE, optional
        DESCRIPTION. The default is False.
    rand : TYPE, optional
        DESCRIPTION. The default is True.
    randn : TYPE, optional
        DESCRIPTION. The default is 100000.

    Returns
    -------
    jsfs : TYPE
        DESCRIPTION.

    """
    pos, gt, pops, pairs = args
    jsfs_list = []
    for pair in pairs:
        i, j = pair.split("-")
        p1 = pops[int(i)]
        p2 = pops[int(j)]
        gtpops = gt.take(p1+p2, axis=1)
        miss_count = gtpops.count_missing(axis=1)
        miss_arr = miss_count > 0
        gtpops = gtpops.compress(miss_arr, axis=0)
        acpops = gtpops.count_alleles()
        segpops = acpops.is_segregating()
        gtseg = gtpops.compress(segpops)
        # random snps
        if rand:
            n = randn  # number of SNPs to choose randomly
            try:
                vidx = np.random.choice(gtseg.shape[0], n, replace=False)
                vidx.sort()
                gtr = gtpops.take(vidx, axis=0)
            except ValueError:
                gtr = gtseg
        else:
            gtr = gtseg
        gtpop1 = gtr.take(range(len(p1)), axis=1)
        gtpop2 = gtr.take(range(len(p1), gtr.shape[1]), axis=1)
        ac1 = gtpop1.count_alleles()
        ac2 = gtpop2.count_alleles()
        #jsfs = allel.joint_sfs(ac1[:, 1], ac2[:, 1])
        # jsfs
        if fold:
            # pad for allel as well
            popsizeA, popsizeB = len(p1)/2, len(p2)/2
            fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
            jsfs = allel.joint_sfs_folded(ac1, ac2)
            fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
        else:
            # pad for allel as well
            popsizeA, popsizeB = len(p1), len(p2)
            fs = np.zeros((popsizeA + 1, popsizeB + 1), dtype=int)
            jsfs = allel.joint_sfs(ac1[:, 1], ac2[:, 1])
            fs[:jsfs.shape[0], :jsfs.shape[1]] = jsfs
        props = summarizejsfs(fs)
        jsfs_list.append(props)
    jsfs = " ".join(map(str, np.concatenate(jsfs_list).ravel()))
    return f"{jsfs}\n"


def calc_afibs(gt, pos, pops, basepairs, fold):
    """Calculate afibs.

    Parameters
    ----------
    gt : TYPE
        DESCRIPTION.
    pos : TYPE
        DESCRIPTION.
    basepairs : TYPE
        DESCRIPTION.
    fold : TYPE
        DESCRIPTION.

    Returns
    -------
    afibs_list : TYPE
        DESCRIPTION.

    """
    subsample = 'all'
    afibs = {}
    for i, pop in enumerate(pops):
        afibs_gtlist = []
        ibsarray = np.zeros((len(pop), len(pop)-1))
        gtpop = gt.take(pop, axis=1)
        miss_count = gtpop.count_missing(axis=1)
        miss_arr = miss_count > 0
        gtpop = gtpop.compress(miss_arr, axis=0)
        pos = pos[miss_arr]
        acpop = gtpop.count_alleles()
        seg = acpop.is_segregating()
        gtseg = gtpop.compress(seg)
        poslist = pos[seg]
        freqlist = np.sum(gtseg, axis=1)
        if subsample == 'all':
            inds_list = range(len(pop))
        else:
            inds_list = np.random.choice(len(pop), subsample, replace=False)
        for ind in inds_list:
            indpos = poslist[gtseg[:, ind] > 0]
            indpos = np.insert(indpos, 0, 0)
            indpos = np.insert(indpos, len(indpos), basepairs)
            for freq in range(1, len(pop)):
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
            afibs_gtlist.append(np.mean(ibsarray, axis=0))
        afibs[i] = np.mean(afibs_gtlist, axis=0)
    # fold and export to list
    afibs_list = []
    if fold:
        for i, pop in enumerate(pops):
            ibs_flip = np.flip(afibs[i], axis=0)
            ibs_fold = (ibs_flip + afibs[i]) / 2
            haps = int(len(pop)/2)
            ibs = ibs_fold[0:haps]
            afibs_list.append(ibs)
    else:
        for i in range(len(pops)):
            afibs_list.append(afibs[i])
    return afibs_list


def afibsStats(args, fold=False):
    """
    For each individual within a population calculate up/down distance to
    nearest SNP for each frequency class of alleles. Average among individuals
    within a population and return a vector of length pop-2 (no fixed classes)
    of sizes for each freq class.

    Parameters
    ----------
    pos : TYPE
        DESCRIPTION.
    hap : TYPE
        DESCRIPTION.
    pops : TYPE
        DESCRIPTION.
    basepairs : TYPE
        DESCRIPTION.
    fold : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    afibs : TYPE
        DESCRIPTION.

    """
    pos, gt, pops, basepairs = args
    afibsmean = calc_afibs(gt, pos, pops, basepairs, fold)
    afibs = " ".join(map(str, [i for t in afibsmean for i in t]))
    return f"{afibs}\n"


# TODO: filet ... make bed, vcf2fasta, calc stats

# def vcf_to_fasta(vcf_path, vcf_files, mask_path, mask_files):
#     """Create fasta from masked and phased VCF.

#     Parameters
#     ----------
#     vcf : TYPE
#         DESCRIPTION.
#     mask_out : TYPE
#         DESCRIPTION.
#     fasta_out1 : TYPE
#         DESCRIPTION.
#     fasta_out2 : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None

#     """
#     # vcf_list = ['Fun_Par_3R.vcf', 'Fun_Van_3R.vcf']
#     for vcf in vcf_list:
#         import phasedVcfsToFastas as vcf2fasta
#         pop1, pop2, arm = vcf.strip("vcf").split("_")
#         mask = f"{pop1}-{pop2}.mask.fa"
#         assert mask in mask_list
#         vcf2fasta(PATHvcf, PATHmask, f"{pop1}.{pop1}-{pop2}.{arm}.masked.fasta", f"{pop2}.{pop1}-{pop2}.{arm}.masked.fasta")
#         fasta_files.append(f"{pop1}-{pop2}-{arm}")
#     return fasta_files


# def make_coords(chrom, chrom_len, start, stop, size, step):
#     """Make windows for FILET.

#     Parameters
#     ----------
#     chrom : TYPE
#         DESCRIPTION.
#     chrom_len : TYPE
#         DESCRIPTION.
#     start : TYPE
#         DESCRIPTION.
#     stop : TYPE
#         DESCRIPTION.
#     size : TYPE
#         DESCRIPTION.
#     step : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     """
#     if start:
#         window = f"{chrom}.{size}-{step}.{start}-{stop}.bed"
#     else:
#         window = f"{chrom}.{size}-{step}.bed"
# def filetStats(args):
#     """

#     Parameters
#     ----------
#     pos : TYPE
#         DESCRIPTION.
#     hap : TYPE
#         DESCRIPTION.
#     pops : TYPE
#         DESCRIPTION.
#     basepairs : TYPE
#         DESCRIPTION.
#     filet_path : TYPE
#         DESCRIPTION.
#     block : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     filet_list : TYPE
#         DESCRIPTION.

#     """
#     pos, hap, pops, basepairs, filet_path, block = args

#     if basepairs > 100000:
#         block = 100000
#     if block == 0:
#         block = basepairs

#     keep_stats = np.array([True, True, True, True, False, True, False, False,
#                            True, True, True, True, True, False, True, False,
#                            False, True, True, True, True, True, True, True,
#                            True, True, True, True, False, False, False])
#     norm = np.array([block, block**2, block, block, 1, 1, block, block**2,
#                      block, block, 1, 1, 1, 1, block, block, block, 1, 1, 1, 1, 1])
#     filet_list = []
#     for pop1, pop2 in combinations(pops, 2):
#         fakems_haps = []
#         n1 = len(pop1)
#         n2 = len(pop2)
#         if type(hap) is list:
#             loci_r = 0
#             for sub_rep in list(zip(pos, hap)):
#                 posr, gtarr = sub_rep
#                 gt = allel.HaplotypeArray(gtarr)
#                 gtpops = gt.take(pop1+pop2, axis=1)
#                 acpops = gtpops.count_alleles()
#                 segpops = acpops.is_segregating()
#                 gtseg = gtpops.compress(segpops)
#                 posit = posr[segpops]
#                 #
#                 if basepairs > block:
#                     start = 0
#                     step = block
#                     end = start + step
#                     while end < basepairs:
#                         loci_r += 1
#                         s_ix = bisect.bisect_left(posit, start)
#                         e_ix = bisect.bisect_right(posit, end) - 1
#                         posit_block = posit[s_ix:e_ix] / basepairs
#                         gtseg_block = gtseg[s_ix:e_ix]
#                         seg = gtseg_block.shape[0]
#                         fakems_haps.append(f"\n//\nsegsites: {seg}\npositions: {' '.join(map(str, posit_block))}\n")
#                         for geno in gtseg_block.transpose():
#                             fakems_haps.append(f"{''.join(map(str, geno))}\n")
#                         start += step
#                         end += step
#                 #
#                 else:
#                     loci_r = len(sub_rep)
#                     posit = posit / block
#                     seg = np.count_nonzero(segpops)
#                     fakems_haps.append(f"\n//\nsegsites: {seg}\npositions: {' '.join(map(str, posit))}\n")
#                     for geno in gtseg.transpose():
#                         fakems_haps.append(f"{''.join(map(str, geno))}\n")
#         else:
#             gt = allel.HaplotypeArray(hap)
#             posr = pos[0]
#             gtpops = gt.take(pop1+pop2, axis=1)
#             acpops = gtpops.count_alleles()
#             segpops = acpops.is_segregating()
#             gtseg = gtpops.compress(segpops)
#             posit = posr[segpops]
#             #
#             if basepairs > block:
#                 loci_r = 0
#                 start = 0
#                 step = block
#                 end = start + step
#                 while end < basepairs:
#                     loci_r += 1
#                     s_ix = bisect.bisect_left(posit, start)
#                     e_ix = bisect.bisect_right(posit, end) - 1
#                     posit_block = posit[s_ix:e_ix] / basepairs
#                     gtseg_block = gtseg[s_ix:e_ix]
#                     seg = gtseg_block.shape[0]
#                     fakems_haps.append(f"\n//\nsegsites: {seg}\npositions: {' '.join(map(str, posit_block))}\n")
#                     for geno in gtseg_block.transpose():
#                         fakems_haps.append(f"{''.join(map(str, geno))}\n")
#                     start += step
#                     end += step
#             #
#             else:
#                 loci_r = 1
#                 posit = posit / block
#                 seg = np.count_nonzero(segpops)
#                 fakems_haps.append(f"\n//\nsegsites: {seg}\npositions: {' '.join(map(str, posit))}\n")
#                 for geno in gtseg.transpose():
#                     fakems_haps.append(f"{''.join(map(str, geno))}\n")
#         fakems_head = f"ms {n1+n2} {loci_r} -t tbs -r tbs {block} -I 2 {n1} {n2}\n1234\n"
#         fakems = "".join(fakems_haps)
#         msinput = fakems_head + fakems
#         filet_prog = os.path.join(filet_path, "twoPopnStats_forML")
#         cmd = [filet_prog, str(n1), str(n2)]
#         proc = run(cmd, stdout=PIPE, input=msinput, encoding='ascii', check=True)
#         # collect stats
#         lstats = proc.stdout.rstrip().split('\n')[1:]
#         stat_vec = [list(map(float, l.split())) for l in lstats]
#         if len(stat_vec) > 1:
#             stat_arr = np.vstack(stat_vec)
#             stat_arr[np.isinf(stat_arr)] = 'nan'
#             filetmean = np.nanmean(stat_arr, axis=0)
#             filet_norm = filetmean[keep_stats] / norm
#         else:
#             stat_arr = np.array(stat_vec[0])[keep_stats]
#             stat_arr[np.isinf(stat_arr)] = 'nan'
#             filet_norm = stat_arr / norm
#         filet_list.append(" ".join(map(str, filet_norm)))
#     return f"{' '.join(filet_list)}\n"

