#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 14:43:07 2018
Class for model building
@author: stsmall
"""
import msprime as msp
import numpy as np
from math import log
from collections import OrderedDict
from collections import defaultdict


class Model(object):
    def __init__(self):
        """
        """
        return(None)

    def genDem(self, ix, params, demodict, parlist, Ne):
        """
        """
        if ix == 1:
            dem_list = self.model1(params, demodict, parlist)
        elif ix == 2:
            dem_list = self.model2(params, demodict, parlist, Ne)
        elif ix == 3:
            dem_list = self.model3(params, demodict, parlist, Ne)

        return(dem_list)

    def model1(self, params, demodict, parlist):
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
        assert len(params) == 8
        dem_list = []
        t01, N01, t43, N43, t23, N23, t13, N13 = params
        timedict = defaultdict(list)
        for i, event in enumerate(parlist):
            if "t" in event:
                p = "N" + event[1:]
                timedict[params[i]].append([params[parlist.index(p)], event])
        dd = defaultdict(list)
        for d in (timedict, demodict):
            for key, value in d.items():
                dd[key].extend(value)
        od = OrderedDict(sorted(dd.items()))
        for k, event in od.items():
            # event = [y for x in event for y in x]
            for v in event:
                if "t" in v[1]:
                    dem_list.append(msp.MassMigration(time=k,
                                                      source=int(v[1][1]),
                                                      destination=int(v[1][2]),
                                                      proportion=1.0))
                    dem_list.append(msp.PopulationParametersChange(time=k,
                                                                   initial_size=v[0],
                                                                   population_id=int(v[1][2])))
                elif "g" in v[1]:
                    dem_list.append(msp.PopulationParametersChange(time=k,
                                                                   growth_rate=v[0],
                                                                   population_id=int(v[1][1])))
        return(dem_list)

    def model2(self, params, demodict, parlist, Ne):
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

        assert len(params) == 10
        isotime = lambda x1, x2: 0 if (x1 - x2) < 0 else x1-x2
        dem_list = []
        t01, N01, t43, N43, t23, N23, t13, N13, mMax, mIso = params
        timedict = defaultdict(list)
        for i, event in enumerate(parlist):
            if "t" in event:
                p = "N" + event[1:]
                timedict[params[i]].append([params[parlist.index(p)], event])

        # calculate these in coal time, then transform into gens
        maxMig = mMax * 4*Ne
        migdict = {}
        for ts in list(timedict.keys()):
            tsp = ts / (4*Ne)
            cur = (ts - 100) / (4*Ne)  # add to iso rather than subtract from speciation
            tis = isotime(ts, mIso) / (4*Ne)
            slope = maxMig / (tsp - tis)
            intercept = (slope*tis)
            # cur = tis + (tis * 0.0005)
            while cur < tsp and cur > tis:
                while True:
                    u = np.random.uniform()
                    v = np.random.uniform()
                    t = cur - (1/maxMig) * log(1-u)
                    if v > (t - tis)/(tsp-tis):
                        break
                migdict[np.round(cur*4*Ne)] = [(slope*cur - intercept)/(4*Ne),
                                               'm' + timedict[ts][0][1][1:]]
                cur -= t - cur
            if ts > mIso:
                migdict[mIso] = [0, 'm' + timedict[ts][0][1][1:]]
        # merge dicts
        dd = defaultdict(list)
        for d in (timedict, demodict, migdict):
            for key, value in d.items():
                dd[key].extend(value)
        od = OrderedDict(sorted(dd.items()))
        for k, event in od.items():
            # event = [y for x in event for y in x]
            for v in event:
                if "t" in v[1]:
                    dem_list.append(msp.MassMigration(time=k,
                                                      source=int(v[1][1]),
                                                      destination=int(v[1][2]),
                                                      proportion=1.0))
                    dem_list.append(msp.PopulationParametersChange(time=k,
                                                                   initial_size=v[0],
                                                                   population_id=int(v[1][2])))
                elif "g" in v[1]:
                    dem_list.append(msp.PopulationParametersChange(time=k,
                                                                   growth_rate=v[0],
                                                                   population_id=int(v[1][1])))
                elif "m" in v[1]:
                    dem_list.append(msp.MigrationRateChange(time=k,
                                                            rate=v[0],
                                                            matrix_index=(int(v[1][1]),
                                                                          int(v[1][2]))))
                    dem_list.append(msp.MigrationRateChange(time=k,
                                                            rate=v[0],
                                                            matrix_index=(int(v[1][2]),
                                                                          int(v[1][1]))))
        return(dem_list)

    def model3(self, params, demodict, parlist, Ne):
        """
        """
        return(None)
