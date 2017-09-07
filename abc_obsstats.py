# -*- coding: utf-8 -*-
"""This script calculated summary statistic from a vcf file. These statistics
are observed summary statistics for ABC inference via ABCrf or abc in R.
Contents:
    LD: using popsizeABC, format VCF
    IBS: using popsizeABC, format VCF

    SFS: using dadi, vcf2dadi
    jSFS: using dadi, vcf2dadi
    Tajima's D: using dadi, vcf2dadi
    Fst: using dadi, vcf2dadi

    bSFS: using ABLE, vcf2able
    jbSFS: using ABLE, vcf2able

    WCFST: popstats.py, vcf2tped
    FST: popstats.py, vcf2tped
    F2: popstats.py, vcf2tped
    F3: popstats.py, vcf2tped
    F4: popstats.py, vcf2tped
"""
from gwas import IO
from gwas import data
from gwas import summary_stat as ss
import numpy as np
import dadi
from collections import defaultdict
from itertools import combinations


class Dadistats:
    """stats from dadi package
    """

    def __init__(self, sfs, jsfs, tajd, fst):
        self.sfs = {}
        self.jsfs = {}
        self.tajd = {}
        self.fst = {}

    def dadi_obs(self, infile, pedfile, mask=True, fold=True):
        """
        """
        # sample sizes in diploid
        # Haiti 7
        # Mali 11
        # Kenya 9
        # PNG 20
        dd = dadi.Misc.make_data_dict(infile)
        peddict = defaultdict(list)
        with open(pedfile, 'r') as ped:
            for line in ped:
                if line.strip():
                    x = line.strip().split()
                    peddict[x[0]].append(x[1])
                else:
                    continue
        sfs = {}
        jsfs = {}
        tajd = {}
        fst = {}
        poplist = peddict.keys()
        size = [len(peddict[p]) for p in poplist]
        for p, s in zip(poplist, size):
            if fold:
                fs = dadi.Spectrum.from_data_dict(dd, pop_ids=p,
                                                  projections=s,
                                                  polarized=False)
                if mask:
                    fs.mask[1] = True
            else:
                fs = dadi.Spectrum.from_data_dict(dd, pop_ids=p,
                                                  projections=s,
                                                  polarized=True)
            sfs[p] = fs
            tajd[p] = fs.Tajima_D()

        for p, s in zip(combinations(poplist, 2), combinations(size, 2)):
            if fold:
                fs = dadi.Spectrum.from_data_dict(dd, pop_ids=p,
                                                  projections=s,
                                                  polarized=False)
                if mask:
                    fs.mask[1, 1] = True
            else:
                fs = dadi.Spectrum.from_data_dict(dd, pop_ids=p,
                                                  projections=s,
                                                  polarized=True)
            jsfs["{}_{}".format(p[0], p[1])] = fs
            fst["{}_{}".format(p[0], p[1])] = fs.Fst()
        self.jsfs = jsfs
        self.sfs = sfs
        self.tajd = tajd
        self.fst = fst
        # 23 summary stats from Naduvilezhath 2011
        s1 = np.sum(fs[0, 1:3])
        s2 = np.sum(fs[1:3, 0])
        s3 = np.sum(fs[0, 3:-3])
        s4 = np.sum(fs[3:-3, 0])
        s5 = np.sum(fs[0, -3:-1])
        s6 = np.sum(fs[-3:-1, 0])
        s7 = np.sum(fs[1:3, 1:3])
        s8 = np.sum(fs[1:3, 3:-3])
        s9 = np.sum(fs[3:-3, 1:3])
        s10 = np.sum(fs[-3:-1, 3:-3])
        s11 = np.sum(fs[3:-3, -3:-1])
        s12 = np.sum(fs[1:3, -3:-1])
        s13 = np.sum(fs[-3:-1, 1:3])
        s14 = np.sum(fs[3:-3, 3:-3])
        s15 = np.sum(fs[-3:-1, -3:-1])
        s16 = np.sum(fs[0, -1])
        s17 = np.sum(fs[-1, 0])
        s18 = np.sum(fs[-1, 1:3])
        s19 = np.sum(fs[1:3, -1])
        s20 = np.sum(fs[-1, 3:-3])
        s21 = np.sum(fs[3:-3, -1])
        s22 = np.sum(fs[-1, -3:-1])
        s23 = np.sum(fs[-3:-1, -1])


class Popsizeabc:
    """Calculates from popsizeABC program. Takes a straight up vcf with no
    intermediate conversion needed
        obsstats(haps, interval_list, chromlist, vcf, list_ani, popsped,
             pop, mac, mac_ld, L, args.out)
    """
    def __init__(self, pop, ld, ldint, ibs):
        self.pop = pop  # x.Popsizeabc("PNG")
        self.ld = ld
        self.ldint = ldint
        self.ibs = ibs

    def ldwindow(self, nb_times=21, r=2.9E-9, L=2E6, per_err=5,
                 Tmax=13000, a=0.06):
        """Creation of the bins of physical distance for which the average
            LD will be computed, based on the time windows defined above.

           Parameters:
                nb_times: int,
                times: array,
                r: float,recomb rate per generation per bp
                L: int, size of each segment, in bp.
                per_err: int,
                Tmax: int,

            Returns:
                intervals_list: list

            interval_list = ldstats(nb_times, times, r, L, per_err, Tmax
            """
        times = -np.ones(shape=nb_times, dtype='float')
        for i in range(nb_times):
            times[i] = (np.exp(np.log(1 + a * Tmax) * i/(nb_times - 1)) - 1)/a
        interval_list = []
        for i in range(nb_times - 1):
            t = (times[i + 1] + times[i])/2
            d = 1/(2 * r * t)
            if d <= L:
                interval_list.append([d - per_err * d/100,
                                      d + per_err * d/100])
        t = Tmax + times[nb_times - 1] - times[nb_times - 2]
        d = 10**8/(2 * t)
        # d = 1/(2*r*t)
        interval_list.append([d-per_err * d/100, d + per_err * d/100])
        return(interval_list)

    def parsevcf(self, vcf, chrlist, popsped, mac=0, mac_ld=0):
        """
        Parameters
        ------
        chromlist: list, list of chromosomes; name chrom.vcf.gz
        vcf: file, vcf file
        popsped: file, ped file

        """
        list_ani = []
        with open(popsped, 'r') as ped:
            for line in ped:
                p = line.strip().split()
                if p[0] in self.pop:
                    list_ani.append(p[1])
        chromlist = []
        with open(chrlist, 'r') as chrm:
            for line in chrm:
                chromlist.append(line.strip())
        count_list = []
        pos_list = []  # for ld
        geno_list = []  # for ld
        nb_snp = 0
        Lchr = 0
        for chrom in chromlist:
            # store data and pre-analyse it
            print("Processing chromosome,{}".format(chrom))
            infile_vcf = "{}/{}.vcf.gz".format(vcf, chrom)
            [mydata, mymap] = IO.parseVcfFile(infile_vcf, includeInd=list_ani)
            pedfile = popsped
            IO.parsePedFile_nogeno(pedfile, mydata)
            u = data.Genos_and_counts(mydata, mymap, self.pop, mac=mac)
            count_list.append(u[1][0])
            p = len(u[0][0])
            nb_snp += p
            Lchr += u[0][0][p - 1]
            u = data.Genos_and_counts(mydata, mymap, self.pop, mac=mac_ld)
            pos_list.append(u[0][0])
            geno_list.append(u[2][0])
        return(pos_list, geno_list, Lchr, nb_snp)

    def ldstats(self, vcf, chrlist, popsped, L=2E6):
        """program to calc summary stats from vcf

            Parameters:
                interval_list: list, list of time intervals
                pop: str, population outfile id
                mac: int, minor allele count
                mac_ld: int, minor allele count LD
                L: int, length of regions

            Returns:
                file: file with string of stats
        """
        size_list = [1]
        prob_list = [0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99,
                     0.999, 0.9999]
        interval_list = self.ldstats()
        self.ldint = interval_list
        pos_list, geno_list, Lchr, nb_snp = self.parsevcf(vcf, chrlist,
                                                          popsped)
        # compute summary statistics
        u = ss.break_chr(pos_list, geno_list, L)
        pos_list = u[0]
        geno_list = u[1]
        res_ld_zyg = ss.distrib_zyg_r2(pos_list, geno_list, interval_list)
        self.ld = res_ld_zyg
        res_ibs = ss.ibs_quantiles_from_geno(size_list[0], pos_list, geno_list,
                                             prob_list, dmax=L)
        for m in size_list[1:]:
            res_ibs = np.concatenate((res_ibs, ss.ibs_quantiles_from_geno(m,
                                      pos_list, geno_list, prob_list, dmax=L)))
        self.ibs = res_ibs
        # fnp = np.array([np.float(nb_snp)/np.float(Lchr)], dtype='float')
        # np.savetxt(fname, np.concatenate((fnp, res_afs, res_ld_zyg[0])),
        #           fmt='%.3e')
