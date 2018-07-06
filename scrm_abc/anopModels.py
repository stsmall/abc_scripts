#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 14:43:07 2018
Class for model building
@author: stsmall
"""
from __future__ import print_function
from __future__ import division
import numpy as np

from collections import OrderedDict
from collections import defaultdict


class Model(object):
    def __init__(self):
        """
        """
        return(None)

    def genDem(self, ix, params, demodict, parlist, Ne, mMax, mIso):
        """
        """
        if ix == 1:
            dem_list = self.model1_scrm(params, demodict, parlist, Ne)
        elif ix == 2:
            dem_list = self.model2_scrm(params, demodict, parlist, Ne, mMax, mIso)
        elif ix == 3:
            dem_list = self.model3_scrm(params, demodict, parlist, Ne, mMax, mIso)
        elif ix == 4:
            dem_list = self.model3_scrm(params, demodict, parlist, Ne, mMax, mIso)
        elif ix == 5:
            dem_list = self.model3_scrm(params, demodict, parlist, Ne, mMax, mIso)
        elif ix == 6:
            dem_list = self.model3_scrm(params, demodict, parlist, Ne, mMax, mIso)

        return(dem_list)

    def model1_scrm(self, params, demodict, parlist, Ne):
        """Simulate speciation followed by isolation

        Parameters
        ------
        params: list, list of parameters
            times, ancestral effective sizes
        demodict: dict, times and growth rate changes for each species
            i.e., time: rate

        Returns
        ------
        dem_list: list, demography list for msprime simulation

        """
        dem_list = []
        timedict = defaultdict(list)
        for i, event in enumerate(parlist):
            if "ej" in event:
                timedict[params[i]].append([event])
        dd = defaultdict(list)
        for d in (timedict, demodict):
            for key, value in d.items():
                dd[key].extend(value)
        od = OrderedDict(sorted(dd.items()))
        sourcelist = []
        for k, event in od.items():
            for v in event:
                if "ej" in v[0]:
                    dem_list.append("-ej {} {} {}".format(k/(4*Ne), v[0][2], v[0][3]))
                    dem_list.append("-eg {} {} {}".format(k/(4*Ne), v[0][2], 0))
                    sourcelist.append(v[0][2])
                elif "Ne" in v[0]:
                    if v[0][2] not in sourcelist:
                        if v[1] == "None":
                            dem_list.append("-eg {} {} {}".format(k/(4*Ne), v[0][2], 4*Ne*float(v[2])))
                        else:
                            dem_list.append("-en {} {} {}".format(k/(4*Ne), v[0][2], int(v[1])/Ne))
                            dem_list.append("-eg {} {} {}".format(k/(4*Ne), v[0][2], 4*Ne*float(v[2])))
        return(dem_list)

    def model2_scrm(self, params, demodict, parlist, Ne, mMax, mIso, t_ints=10):
        """Simulate speciation followed by migration

        Parameters
        ------
        params: list, list of parameters
            times, ancestral effective sizes, time_isolation, migrationMax
        demodict: dict, times and growth rate changes for each species
            i.e., time: rate

        Returns
        ------
        dem_list: list, demography list for msprime simulation

        """

        dem_list = []
        timedict = defaultdict(list)
        for i, event in enumerate(parlist):
            if "ej" in event:
                timedict[params[i]].append([event])
        # calculate these in coal time, then transform into gens
        NemMax = mMax * 4*Ne
        migdict = defaultdict(list)
        tsc = mIso/(4*Ne)
        for td in list(timedict.keys()):
            # ts = mIso_list  # list of different speciation/isolation times
            tdc = td/(4*Ne)
            if tdc-tsc < 0:
                tlin = np.linspace(0, tdc, t_ints)
                mlist = [(NemMax/(tdc+tsc))*(t+tsc) for t in tlin]
            else:
                tlin = np.linspace(tdc-tsc, tdc, t_ints)
                #tlist = [np.round(t) for t in tlin]
                mlist = [(NemMax/(tdc-(tdc-tsc)))*(t-(tdc-tsc)) for t in tlin]
            for t, m in zip(tlin, mlist):
                migdict[np.round(t*4*Ne)].append(['mg' + timedict[td][0][0][2:], m/(4*Ne)])
            # migdict[mIso] = ['mg' + timedict[td][0][2:], 0]
        # merge dicts
        dd = defaultdict(list)
        for d in (timedict, demodict, migdict):
            for key, value in d.items():
                dd[key].extend(value)
        od = OrderedDict(sorted(dd.items()))
        sourcelist = []
        for k, event in od.items():
            # event = [y for x in event for y in x]
            for v in event:
                if "ej" in v[0]:
                    dem_list.append("-ej {} {} {}".format(k/(4*Ne), v[0][2], v[0][3]))
                    dem_list.append("-eg {} {} {}".format(k/(4*Ne), v[0][2], 0))
                    dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], 0))
                    dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], 0))
                    sourcelist.append(v[0][2])
                elif "Ne" in v[0]:
                    if v[0][2] not in sourcelist:
                        if v[1] == "None":
                            dem_list.append("-eg {} {} {}".format(k/(4*Ne), v[0][2], 4*Ne*float(v[2])))
                        else:
                            dem_list.append("-en {} {} {}".format(k/(4*Ne), v[0][2], int(v[1])/Ne))
                            dem_list.append("-eg {} {} {}".format(k/(4*Ne), v[0][2], 4*Ne*float(v[2])))
                elif "mg" in v[0]:
                    if v[0][2] in sourcelist or v[0][3] in sourcelist:
                        pass
                    else:
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], v[1]*4*Ne))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], v[1]*4*Ne))
        return(dem_list)

    def model3_scrm(self, params, demodict, parlist, Ne):
        """
        """
        # migration from demodict or params is listed as m, where as gradual speciation is listed as mg
        return(None)
