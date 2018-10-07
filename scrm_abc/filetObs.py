#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 17:46:23 2018

@author: scott
"""
from __future__ import print_function
from __future__ import division
import glob
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', "--pair", help="pair")
args = parser.parse_args()

chrmfiles = glob.glob("{}.*.obsstats.out".format(args.pair))
v = []
for file in chrmfiles:
    print(file)
    with open(file, 'r') as stats:
        next(stats)
        try:
            for line in stats:
                x = line.split()
                xx = np.array(map(float, x[4:]))
                if any(np.isinf(xx)) or any(np.isnan(xx)):
                    continue
                v.append(xx)
        except StopIteration:
            break
p = np.mean(np.vstack(v), axis=0)
r = np.std(np.vstack(v), axis=0)
print("{}".format(" ".join(map(str, p))))
print("{}".format(" ".join(map(str, r))))
