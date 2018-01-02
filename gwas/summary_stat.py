# -*- coding: utf-8 -*-
"""
This file includes functions to compute summary statistics from haplotype
or genotype data, or to manipulate such data (for instance converting
haplotype data into genotype data).

4 types of data are generally used as input for the functions here :
    - polymorphic positions : these are lists of 1D arrays, where each 1D
        array provides physical positions for one chromosome
    - minor allele counts : these are lists of 1D arrays, where each 1D array
        provides the minor allele counts for the polymorphic positions found
        on one chromosome
    - haplotypes : these are lists of 2D arrays, where each 2D array provides
        the haplotypes for the polymorphic positions found on one chromosome.
        Each line corresponds to an haplotype and each column to a polymorphic
        site. Alleles are coded 0 and 1
    - genotypes : these are lists of 2D arrays, where each 2D array provides
        the genotypes for the polymorphic positions found on one chromosome.
        Each line corresponds to an individual and each column to a polymorphic
        site. Genotypes are coded 0, 1 and 2
"""

import data
import numpy as np
import scipy as sp
import warnings
from scipy import stats
from gwas import missing, complete_cases, is_na


def spatial_histo(pos_list, count_list, M, dmax=np.inf):
    '''

    Returns
    ------
    -the histogram of values in counts, from 1 to M
    -the variance of the distance between sites with count i, for all i.
        This variance needs to be multiplied by the overall proportion of SNPs.
        positions of vector pos are assumed to be SORTED within each
        chromosome
    '''
    histo = np.zeros(shape=M, dtype='float')
    dist = -np.ones(shape=M, dtype='float')
    nb_snp = 0
    for i in range(1, M + 1):
        d = np.zeros(shape=1, dtype='int32')
        for chro in range(0, len(pos_list)):
            sel = (count_list[chro] == i)
            histo[i-1] += sum(sel)
            pos_temp = pos_list[chro]
            pos_temp = pos_temp[sel]
            d_temp = pos_temp[1:] - pos_temp[:len(pos_temp)-1]
            d = np.concatenate((d, d_temp[(d_temp <= dmax)]))
        if len(d) > 2:
            dist[i-1] = np.std(d[1:])
    return(histo/np.sum(histo), dist)


def histo(count_list, M):
    '''
    Returns
    ------
        - the histogram of values in counts, from 1 to M
    '''
    histo = np.zeros(shape=M, dtype='float')
    for i in range(1, M+1):
        for chro in range(0, len(count_list)):
            sel = (count_list[chro] == i)
            histo[i-1] += sum(sel)
    return(histo/np.sum(histo))


def r2(u, v):
    '''
    Returns
    ------
    the r2 value for two haplotype vectors (numpy arrays with alleles coded 1
    and 0)
    '''
    fcross = np.mean(u*v)
    fu = np.mean(u)
    fv = np.mean(v)
    return((fcross-fu*fv)**2/(fu*(1-fu)*fv*(1-fv)))


def distrib_r2(pos_list, hap_list, interval_list):
    '''
    Returns
    ------
    the mean and the variance of r2 for a list of distance intervals.
    - pos_list is a list of 1 dim arrays
    - hap_list is a list of 2 dim arrays
    - interval_list is a list of ordered pairs
    a subset of non overlapping pairs is used for each interval
    '''
    p = len(interval_list)
    moy = -np.ones(shape=p, dtype='float')
    var = -np.ones(shape=p, dtype='float')
    for i in range(0, p):
        r2_list = []
        dmin = interval_list[i][0]
        dmax = interval_list[i][1]
        # looks for snp pairs with the good distance
        for chro in range(0, len(pos_list)):
            nb_snp = len(pos_list[chro])
            if nb_snp > 0:
                i_deb = 0
                i_fin = 1
                while i_fin < nb_snp:
                    while i_fin < nb_snp and (pos_list[chro][i_fin]
                                              - pos_list[chro][i_deb]) < dmin:
                        i_fin += 1
                    if i_fin < nb_snp and (pos_list[chro][i_fin]
                                           - pos_list[chro][i_deb]) <= dmax:
                        # compute r2
                        u_deb = hap_list[chro][:, i_deb]
                        u_fin = hap_list[chro][:, i_fin]
                        r2_list.append(r2(u_deb, u_fin))
                    i_deb = i_fin+1
                    i_fin = i_deb+1
        if len(r2_list) < 2:
            # try a more exhaustive screening of SNP pairs
            r2_list = []
            dmin = interval_list[i][0]
            dmax = interval_list[i][1]
            for chro in range(0, len(pos_list)):
                nb_snp = len(pos_list[chro])
                if nb_snp > 0:
                    i_deb = 0
                    i_fin = 1
                    while i_fin < nb_snp:
                        while i_fin < nb_snp and (pos_list[chro][i_fin]
                                                  - pos_list[chro][i_deb]) < dmin:
                            i_fin += 1
                        if i_fin < nb_snp and (pos_list[chro][i_fin]
                                               - pos_list[chro][i_deb]) <= dmax:
                            # compute r2
                            u_deb = hap_list[chro][:, i_deb]
                            u_fin = hap_list[chro][:, i_fin]
                            r2_list.append(r2(u_deb, u_fin))
                        i_deb += 1
                        i_fin = i_deb+1
        # computes the stat
        if len(r2_list) >= 2:
            moy[i] = np.mean(np.array(r2_list, dtype='float'))
            var[i] = np.std(np.array(r2_list, dtype='float'))
    return(moy, var)


def zyg_r2(u, v):
    '''
    Returns
    ------
    -the zygotic r2 value for two genotype vectors
        (numpy arrays with genotypes coded 0, 1 and 2)
    '''
    return((np.corrcoef(u, v)[0, 1])**2)


def distrib_zyg_r2(pos_list, geno_list, interval_list):
    '''
    Returns
    ------
    the mean and the variance of zygotic r2 for a list of distance intervals.
    pos_list is a list of 1 dim arrays
    geno_list is a list of 2 dim arrays
    interval_list is a list of ordered pairs
    a subset of non overlapping pairs is used for each interval
    SNP where all individuals are heterozygotes are dealt by catching a
        RuntimeWarning, BUT this may not work with old numpy versions
    '''
    warnings.simplefilter("error", RuntimeWarning)
    p = len(interval_list)
    moy = -np.ones(shape=p, dtype='float')
    var = -np.ones(shape=p, dtype='float')
    for i in range(0, p):
        r2_list = []
        dmin = interval_list[i][0]
        dmax = interval_list[i][1]
        for chro in range(0, len(pos_list)):
            nb_snp = len(pos_list[chro])
            if nb_snp > 0:
                i_deb = 0
                i_fin = 1
                while i_fin < nb_snp:
                    while i_fin < nb_snp and (pos_list[chro][i_fin]
                                              - pos_list[chro][i_deb]) < dmin:
                        i_fin += 1
                    if i_fin < nb_snp and (pos_list[chro][i_fin]
                                           - pos_list[chro][i_deb]) <= dmax:
                        # compute r2
                        u_deb = geno_list[chro][:, i_deb]
                        u_fin = geno_list[chro][:, i_fin]
                        try:
                            r2_list.append(zyg_r2(u_deb, u_fin))
                        except RuntimeWarning:
                            pass
                    i_deb = i_fin + 1
                    i_fin = i_deb + 1
        if len(r2_list) < 2:
            # try a more exhaustive screening of SNP pairs
            r2_list = []
            dmin = interval_list[i][0]
            dmax = interval_list[i][1]
            for chro in range(0, len(pos_list)):
                nb_snp = len(pos_list[chro])
                if nb_snp > 0:
                    i_deb = 0
                    i_fin = 1
                    while i_fin < nb_snp:
                        while i_fin < nb_snp and (pos_list[chro][i_fin]
                                                  - pos_list[chro][i_deb]) < dmin:
                            i_fin += 1
                        if i_fin < nb_snp and (pos_list[chro][i_fin]
                                               - pos_list[chro][i_deb]) <= dmax:
                            # compute r2
                            u_deb = geno_list[chro][:, i_deb]
                            u_fin = geno_list[chro][:, i_fin]
                            try:
                                r2_list.append(zyg_r2(u_deb, u_fin))
                            except RuntimeWarning:
                                pass
                        i_deb += 1
                        i_fin = i_deb + 1
        if len(r2_list) >= 2:
            moy[i] = np.mean(np.array(r2_list, dtype='float'))
            var[i] = np.std(np.array(r2_list, dtype='float'))
    return(moy, var)


def break_chr(pos_list, hap_list, dmax=2000000):
    '''Breaks a list of long chromosomes into an equivalent list of chromosomes
    with length lower than dmax. To be used before ibs_quantiles in the case of
    real data sets with unequal chromosomes legnths
    '''
    pos_list_new = []
    hap_list_new = []
    for chro in range(0, len(pos_list)):
        pos_temp = pos_list[chro]
        hap_temp = hap_list[chro]
        outlier_ind = (pos_temp > dmax)
        while np.sum(outlier_ind) > 0:
            if np.prod(outlier_ind) == 0:
                pos_list_new.append(pos_temp[np.logical_not(outlier_ind)])
                hap_list_new.append(hap_temp[:, np.logical_not(outlier_ind)])
            pos_temp = pos_temp[outlier_ind]-dmax
            hap_temp = hap_temp[:, outlier_ind]
            outlier_ind = (pos_temp > dmax)
        pos_list_new.append(pos_temp)
        hap_list_new.append(hap_temp)
    return(pos_list_new, hap_list_new)


def hap_to_geno(hap_list):
    '''Transforms a list of haplotypes into a list of genotypes
    pairs of haplotypes are randomly sampled for each chromosome
    '''
    geno_list = []
    for hap in hap_list:
        n = hap.shape[0]
        p = hap.shape[1]
        permut = np.random.permutation(n)
        geno = -np.ones(shape=(n/2, p), dtype='int32')
        for i in range(n/2):
            geno[i, :] = hap[permut[2*i], :]+hap[permut[2*i+1], :]
        geno_list.append(geno)
    return(geno_list)


def hap_to_geno_sorted(hap_list):
    '''Transforms a list of haplotypes into a list of genotypes
    pairs of haplotypes follow the order of lines
    (lines 1 and 2 = indiv 1, ...)
    '''
    geno_list = []
    for hap in hap_list:
        n = hap.shape[0]
        p = hap.shape[1]
        geno = -np.ones(shape=(n/2, p), dtype='int32')
        for i in range(n/2):
            geno[i, :] = hap[2*i, :] + hap[2*i+1, :]
        geno_list.append(geno)
    return(geno_list)
