#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:04:02 2018

@author: scott
"""
import sys
f = open("{}.format".format(sys.argv[1]), 'w')
f.write("theta,t,p\n")
with open(sys.argv[1], 'r') as abc:
    for line in abc:
        if line.startswith("-t"):
            theta = line.split()[1]
            line = abc.next()
            if line.startswith("-ev"):
                ev = line.split()
                t = ev[1]
                p = ev[4]
                f.write("{},{},{}\n".format(theta, t, p))
            else:
                f.write("{},{},{}\n".format(theta, "NA", "NA"))
                theta = line.split()[1]
                line = abc.next()
                ev = line.split()
                t = ev[1]
                p = ev[4]
                f.write("{},{},{}\n".format(theta, t, p))
f.close()
