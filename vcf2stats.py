#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 16:05:25 2017

@author: scott
"""

import numpy as np
import argparse
from abc_obsstats import Popsizeabc
from abc_obsstats import Dadistats
from vcf2dadi import parse_vcf

parser = argparse.ArgumentParser()
parser.add_argument('-vcf', "--vcffile", type=str, required=True,
                    help="path to vcf directory, vcfs should be named as"
                    "chr.vcf.gz")
parser.add_argument('-ped', "--pedfile", type=str, required=True,
                    help="path to pedfile")
parser.add_argument('-chr', "--chrlist", type=str, required=True,
                    help="path to chromfile")
args = parser.parse_args()

# LD, IBS :  png.ldint, png.ld, png.ibs
png = Popsizeabc("PNG")
png.ldstats(args.vcf, args.chrlist, args.pedfile)

# SFS, jSFS, tajd, fst : fs.sfs, fs.jsfs, fs.tajd, fs.fst
parse_vcf(args.vcf, args.pedfile, "Wb", "Baml")
fs = Dadistats()
fs.dadi_obs("dadi.notriplets.in", args.pedfile)

# bSFS, bjSFS : bfs.sfs, bfs.jsfs
vcf2able()

# td, f2, f3, f4 : pop.Tdiv, pop.F2, pop.F3, pop.F4
vcf2tped()
