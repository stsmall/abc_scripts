#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 12:12:54 2017

@author: scott
"""
import sys
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--vcf", type=str, required=True,
                    help="name of infile vcf")
parser.add_argument('-o', "--out", type=str, required=True,
                    help="name of outfile")
parser.add_argument('-p', "--ped", type=str, required=True,
                    help="name of pedfile")
args = parser.parse_args()

count = 0
progress_every = 10000
pedf = open(args.out + '.tped', 'w')

with open(args.vcf, 'r') as vcfin:
    for line in vcfin:
        if line.startswith("#CHROM"):
            line = vcfin.next()
            chrom0 = line.strip().split()[0]
            break
chrom_num = 1
with open(args.vcf, 'r') as vcfin:
    for line in vcfin:
        line = line.strip()
        if not line.startswith("##"):
            if line.startswith('#CHROM'):
                subject_ids = line.strip().split()[9:]
            else:
                if (count + 1) % progress_every == 0:
                    sys.stderr.write(str(count+1) + '\t')
                    sys.stderr.flush()
                count += 1
                words = line.strip().split()
                chrom = words[0]
                pos = words[1]
                _id = words[2]
                ref = words[3]  # ref allele
                alt = words[4].split(',')[0]  # alt allele
                # if _id is '.' give sequential id
                if _id == '.':
                    _id = "{}:{}".format(chrom, pos)
                if type(chrom) is str:
                    if chrom0 == chrom:
                        chrom = str(chrom_num)
                    else:
                        chrom_num += 1
                        chrom = str(chrom_num)
                    chrom0 = words[0]
                out = [chrom, _id, '0', pos]
                alleles = [ref, alt]
                for i in xrange(9, len(words)):
                    gt = words[i].split(':')[0]
                    if '.' in gt[0]:
                        gt = "{}\t{}".format(0, 0)
                    else:
                        gt0 = re.split('/|\|', gt)
                        gt = alleles[int(gt0[0])] + '\t' + alleles[int(gt0[1])]
                    out.append(gt)
                pedf.write('\t'.join(out) + '\n')
pedf.close()
if count > progress_every:
    sys.stderr.write('\n')

# write fam file
inds = []
famf = open(args.out + '.tfam', 'w')
with open(args.ped, 'r') as ped:
    for line in ped:
        if line.strip().split()[1] in subject_ids:
            famf.write("{}\n".format("\t".join(line.strip().split())))
            inds.append(line.strip().split()[1])
    temp3 = [x for x in subject_ids if x not in set(inds)]
    temp4 = [x for x in inds if x not in set(subject_ids)]
    if temp3 or temp4:
        raise AssertionError("Ind {} in vcf, but not in ped".format(temp3))
famf.close()
