#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 16:24:18 2018

@author: stsmall
"""
import numpy as np
import allel
import argparse
import glob
from anopStats import SimStats
assert allel.__version__ == "1.1.10"

parser = argparse.ArgumentParser()
parser.add_argument('-base', "--baseName", required=True,
                    help="Base Name of input files for glob")
parser.add_argument('-m', "--model_ix", type=int, required=True,
                    help="model index")
parser.add_argument("--filet", action="store_true", help="use filet to calc"
                    " summary statistics for all pairs")
args = parser.parse_args()


def read_msformat_file(msFile):
    """
    """
    pos_count = 0
    pos_list = []
    gt_list = []
    with open(f, 'r') as ms:
        for line in ms:
            if line.startswith("scrm"):
                x = line.split()
                nind = int(x[1])
                block = int(x[7])
                pops = int(x[9])
                popsize = map(int, x[10:10+pops])
                ind = 0
                popconfig = []
                for p in popsize:
                    popconfig.append(list(range(ind, p+ind)))
                    ind = p
            elif line.startswith("positions"):
                # collisions can result here when theta is high
                pos = np.round(np.array(line.strip().split()[1:], dtype=np.float64))
                prev = 0
                # :TODO double check if this works correctly
                for idx, item in enumerate(pos, start=0):
                    while prev >= item:
                        item += 1
                    pos[idx] = item
                    prev = pos[idx]
                pos_list.append(pos.astype(np.int64) + pos_count)  # append
                pos_count += block + 10000
                line = next(ms)
                gt = np.zeros((nind, pos.shape[0]), dtype=np.uint8)
                cix = 0
                try:
                    while line:
                        line = list(line.strip())
                        try:
                            gt[cix, :] = np.array(line, dtype=np.uint8)
                        except IndexError:
                            break
                        cix += 1
                        line = next(ms)
                except StopIteration:
                    gt_list.append(gt)
                    break
                gt_list.append(gt)
        return(gt_list, np.concatenate(pos_list, axis=0), popconfig, block)


def read_msformat(msFile):
    """
    """
    pos_count = 0
    pos_list = []
    gt_list = []
    with open(f, 'r') as ms:
        for line in ms:
            if line.startswith("scrm"):
                x = line.split()
                nind = int(x[1])
                block = int(x[7])
                pops = int(x[9])
                popsize = map(int, x[10:10+pops])
                ind = 0
                popconfig = []
                for p in popsize:
                    popconfig.append(list(range(ind, p+ind)))
                    ind = p
            elif line.startswith("positions"):
                # collisions can result here when theta is high
                pos = np.round(np.array(line.strip().split()[1:], dtype=np.float64))
                prev = 0
                # :TODO double check if this works correctly
                for idx, item in enumerate(pos, start=0):
                    while prev >= item:
                        item += 1
                    pos[idx] = item
                    prev = pos[idx]
                pos_list.append(pos.astype(np.int64) + pos_count)  # append
                pos_count += block + 10000
                line = next(ms)
                gt = np.zeros((nind, pos.shape[0]), dtype=np.uint8)
                cix = 0
                try:
                    while line:
                        line = list(line.strip())
                        try:
                            gt[cix, :] = np.array(line, dtype=np.uint8)
                        except IndexError:
                            break
                        cix += 1
                        line = next(ms)
                except StopIteration:
                    gt_list.append(gt)
                    break
                gt_list.append(gt)
        return(gt_list, np.concatenate(pos_list, axis=0), popconfig, block)


if __name__ == "__main__":
    input_file = glob.glob("{}*".format(args.baseName))
    for f in input_file:
        gtlist, pos, pops, block = read_msformat(f)
        # allel
        s = SimStats(gtlist, pops, pos)
        jsfs = s.jsfsStats(fold=False)
        asfs = s.asfsStats(fold=False)
        if args.filet:
            filet = s.filetStats(block, '/home/scott/programs_that_work/FILET')
