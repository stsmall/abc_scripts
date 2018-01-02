# -*- coding: utf-8 -*-
"""This script calculated summary statistic from a vcf file. These statistics
are observed summary statistics for ABC inference via ABCrf or abc in R.
"""
from __future__ import print_function
from gwas import IO
from gwas import data
from gwas import summary_stat as ss
import numpy as np
import gzip
from collections import defaultdict
from itertools import combinations
import dadi


class SFSstats:
    """Stats from dadi package
    """

    def __init__(self, sfs, jsfs, tajd, fst):
        self.sfs = {}
        self.jsfs = {}
        self.tajd = {}
        self.fst = {}

    def sfs_obs(self, infile, pedfile, mask=True, fold=True):
        """Calculates stats from dadi infile
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

    def jSFS_summary(self):
        """23 summary stats from Naduvilezhath 2011
        """
        jsfsdict = {}
        for pair in self.jsfs.keys():
            jsfsarray = np.zeros(24)
            fs = self.jsfs[pair]
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
            # jsfsarray[23] = self.fst[pair]
            jsfsdict[pair] = jsfsarray
        return(jsfsdict)

    def summarywrite(self, jsfsdict):
        """Write summary stats to output files
        """
        for pair in jsfsdict.keys():
            np.savetxt("{}.jsfs-fst".format(pair), jsfsdict[pair])
#        for pop in self.sfs.keys():
#            np.savetxt("{}.sfs".format(pop), self.sfs[pop])


class Popsizeabc:
    """Calculates from popsizeABC program. Takes a straight up vcf with no
    intermediate conversion needed
        obsstats(haps, interval_list, chromlist, vcf, list_ani, popsped,
             pop, mac, mac_ld, L, args.out)
    """
    def __init__(self, pop, afs=None, ld=None, ldints=None, ibs=None):
        self.pop = pop  # x.Popsizeabc("PNG")
        self.afs = afs
        self.ld = ld
        self.ldints = ldints
        self.ibs = ibs

    def ldwindow(self, configdict):
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
        nb_times = configdict["nb_times"]
        r = configdict["recombination_rate"]
        L = configdict["contig_length"]
        per_err = configdict["pererr"]
        Tmax = configdict["tmax"]
        a = configdict["a"]
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
        # d = 10**8/(2 * t)
        d = 1/(2*r*t)
        interval_list.append([d-per_err * d/100, d + per_err * d/100])
        self.ldints = interval_list
        return(interval_list)

    def parsevcf(self, vcfFile, chrlist, popsped, mac, mac_ld):
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
        chrlen = []
        with open(chrlist, 'r') as chrm:
            for line in chrm:
                chromlist.append(line.strip().split()[0])
                chrlen.append(int(line.strip().split()[1]))
        count_list = []
        pos_list = []  # for ld
        geno_list = []  # for ld
        nb_snp = 0
        Lchr = 0
        for chrom in chromlist:
            # store data and pre-analyse it
            f = gzip.open("{}.{}.vcf.gz".format(vcfFile, chrom), 'w')
            print("Processing chromosome,{}".format(chrom))
            print(chrom)
            with gzip.open("{}.vcf.gz".format(vcfFile), 'r') as vcf:
                for line in vcf:
                    if line.startswith("#"):
                        f.write(line)
                    else:
                        x = line.strip().split()
                        if x[0] == chrom:
                            f.write(line)
            f.close()
            infile_vcf = "{}.{}.vcf.gz".format(vcfFile, chrom)
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
        return(pos_list, geno_list, count_list, Lchr, nb_snp, chrlen)

    def ldstats(self, vcf, chrlist, popsped, configdict, pix):
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
        L = configdict["contig_length"]
        interval_list = self.ldwindow(configdict)
        pos_list, geno_list, count_list, Lchr, nb_snp, chrlen = self.parsevcf(vcf, chrlist,
                                                                              popsped, configdict["mac"],
                                                                              configdict["mac_ld"])
        haps = configdict["sample_size"][pix]
        res_afs = ss.histo(count_list, int(haps/2))
        # afs = np.insert(res_afs, 0, np.float(nb_snp)/np.float(Lchr))
        u = ss.break_chr(pos_list, geno_list, L)
        pos_list = u[0]
        geno_list = u[1]
        res_ld_zyg = ss.distrib_zyg_r2(pos_list, geno_list, interval_list)
        self.ld = res_ld_zyg
        self.afs = res_afs
        ibs = self.ibsStats(self, pos_list, geno_list, chrlen)
        self.ibs = ibs
        return(None)

    def ibsStats(self, pos_list, geno_list, chrlen, fold=True):
        """program to calc summary stats from vcf

            Parameters:
                pos_list: array, array of position of each snp
                geno_list: array, array of genotypes
                Lchr: int, total bases
                nb_snp: int, total snps

            Returns:
                file: file with string of stats
        """
        inds = geno_list[0].shape[0]
        ibssfslist = []
        for i in range(0, len(geno_list)):
            ibslist = []
            for freq in range(1, inds * 2):
                ibs = 0
                # snp position only
                het = np.sum(geno_list[i], axis=0) > 0  # mask
                # positions of snps
                poslist = pos_list[i][het]  # genome positions not index
                # freq of snp positions
                freqlist = np.sum(geno_list[i], axis=0)[het]
                mut_ix = np.where(freqlist == freq)[0]
                for m in mut_ix:
                    if poslist[m] == 0:
                        ibs += poslist[m + 1]
                    else:
                        try:
                            ibs += poslist[m + 1] - poslist[m - 1]
                        except IndexError:
                            ibs += chrlen[i] - poslist[m - 1]
                try:
                    ibslist.append(ibs / len(mut_ix))
                except ZeroDivisionError:
                    ibslist.append(0)
            ibssfslist.append(ibslist)
        if fold:
            ibs_mean = np.mean(ibssfslist, axis=0)
            ibs_meanr = np.flip(ibs_mean, axis=0)
            ibs_fold = (ibs_mean + ibs_meanr) / 2
            ibs = ibs_fold[0:inds]
        else:
            ibs = np.mean(ibssfslist, axis=0)
        return(ibs)

    def printstats(self, pop):
        """Prints stats per popualtion from vcf
        """
        np.savetxt("pop.afs-ld.obs", np.concatenate((self.afs, self.ld[0])),
                   fmt='%.3e')
        print("{}".format(self.ldints))
        return(None)
