#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 12:30:50 2018

@author: scott
"""
from __future__ import print_function
from __future__ import division
import numpy as np
import glob

# TODO: Allow reading of msfiles from simulate_var
def read_msformat(msout):
    """Read and parse file of ms coalescent simulations from stdout
    """
    pos_count = 0
    pos_list = []
    gt_list = []
    ms = iter(msout.stdout.readline, '')
    for line in ms:
        if line != b'':
            line = line.decode('utf-8')
            if line.startswith("scrm"):
                x = line.split()
                nind = int(x[1])
                block = int(x[7])
                pops = int(x[9])
                popsize = map(int, x[10:10+pops])
                ind = 0
                popconfig = []
                scrmline = line
                for p in popsize:
                    if p == 0:
                        pass
                    else:
                        popconfig.append(list(range(ind, p+ind)))
                        ind += p
                line = next(ms)
                seed = line.split()[0].decode('utf-8')
            elif line.startswith("positions"):
                pos = np.array(line.strip().split()[1:], dtype=np.float64)
                pos_list.append(pos)
                # format as integer for scikit-allel
                # need to re-add -SC abs
                # need to enable posit = np.concatenate(pos_list, axis=0)
#                pos = np.round(np.array(line.strip().split()[1:], dtype=np.float64))
#                prev = 0
#                # collisions can result here when theta is high
#                for idx, item in enumerate(pos, start=0):
#                    while prev >= item:
#                        item += 1
#                    pos[idx] = item
#                    prev = pos[idx]
#                pos_list.append(pos.astype(np.int64) * pos_count)  # append
#                pos_count += block + 1
                line = next(ms).decode('utf-8')
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
                        line = next(ms).decode('utf-8')
                except StopIteration:
                    gt_list.append(gt)
                    break
                gt_list.append(gt)
        else:
            break
#    posit = np.concatenate(pos_list, axis=0)
    posit = pos_list
    return(gt_list, posit, popconfig, block, scrmline, seed)


def read_msformat_file(baseName):
    """Read and parse simulations from ms formatted file from folder of files
    with basename
    """
    input_file = glob.glob("{}*".format(baseName))
    for f in input_file:
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
                    scrmline = line
                    seed = int(next(ms.split()[0]))
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
        return(gt_list, np.concatenate(pos_list, axis=0), popconfig, block, scrmline, seed)
