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
import collections


def model1(params, growthdict):
    """Simulate speciation followed by isolation

    Parameters
    ------
    params: list, list of parameters
        times, ancestral effective sizes
    growthdict: dict, times and growth rate changes for each species
        i.e., time: rate

    Returns
    ------
    dem_list: list, demography list for msprime simulation

    """

    assert len(params) == 8
    dem_list = []
    t_01, N_01, t_43, N_43, t_23, N_23, t_13, N_13 = params
    timedict = {t_01: [N_01, "t01"],
                t_43: [N_43, "t43"],
                t_23: [N_23, "t23"],
                t_13: [N_13, "t13"]
                }
    d = {**growthdict, **timedict}
    od = collections.OrderedDict(sorted(d.items()))
    for k, v in od.items():
        if "t" in v[1]:
            dem_list.append(msp.MassMigration(time=k,
                                              source=v[1][1],
                                              destination=v[1][2],
                                              proportion=1.0))
            dem_list.append(msp.PopulationParametersChange(time=k,
                                                           initial_size=v[0],
                                                           population_id=v[1][2]))
        elif "g" in v[1]:
            dem_list.append(msp.PopulationParametersChange(time=k,
                                                           growth_rate=v[0],
                                                           population_id=v[1][1]))
    return(dem_list)


def model2(params, growthdict):
    """Simulate speciation followed by migration

    Parameters
    ------
    params: list, list of parameters
        times, ancestral effective sizes, time_isolation, migrationMax
    growthdict: dict, times and growth rate changes for each species
        i.e., time: rate

    Returns
    ------
    dem_list: list, demography list for msprime simulation

    """

    assert len(params) == 10
    dem_list = []
    t_01, N_01, t_43, N_43, t_23, N_23, t_13, N_13, mMax, mIso = params
    timedict = {t_01: [N_01, "t01"],
                t_43: [N_43, "t43"],
                t_23: [N_23, "t23"],
                t_13: [N_13, "t13"]
                }
    # calculate these in coal time, then transform into gens
    tis = mIso
    maxMig = mMax
    migdict = {}
    for tsp in list(timedict.keys()):
        slope = maxMig / (tsp - tis)
        intercept = (slope*tis)
        cur = tsp - .05
        while cur < tsp and cur > tis:
            while True:
                u = np.random.uniform()
                v = np.random.uniform()
                t = cur - (1/maxMig) * log(1-u)
                if v > (t - tis)/(tsp-tis):
                    break
            migdict[cur] = [slope*cur - intercept, 'm' + timedict[tsp][1][1:]]
            cur -= t - cur
        migdict[mIso] = [0, 'm' + timedict[tsp][1][1:]]
    d = {**growthdict, **timedict, **migdict}
    od = collections.OrderedDict(sorted(d.items()))
    for k, v in od.items():
        if "t" in v[1]:
            dem_list.append(msp.MassMigration(time=k,
                                              source=v[1][1],
                                              destination=v[1][2],
                                              proportion=1.0))
            dem_list.append(msp.PopulationParametersChange(time=k,
                                                           initial_size=v[0],
                                                           population_id=v[1][2]))
        elif "g" in v[1]:
            dem_list.append(msp.PopulationParametersChange(time=k,
                                                           growth_rate=v[0],
                                                           population_id=v[1][1]))
        elif "m" in v[1]:
            dem_list.append(msp.MigrationRateChange(time=k,
                                                    rate=v[0],
                                                    matrix_index=(int(v[1][1]), int(v[1][2]))))
            dem_list.append(msp.MigrationRateChange(time=k,
                                                    rate=v[0],
                                                    matrix_index=(int(v[1][2]), int(v[1][1]))))
    return(dem_list)


def setupModel(m_index, params, growthdict):
    """Return model matching index

    model 1: divergence times and effective sizes
    model 2: divergence times and effective sizes with partial migration
    model 3: div times, eff sizes, partial migration, 1 admix: Van
    model 4: div times, eff sizes, partial migration, 2 admix: Van, Long
    model 5: div times, eff sizes, partial migration, 3 admix: Van, Long, Anc
    model 6: div times, eff sizes, partial migration, 4 admix: Van, Long, Anc, Like
    """

    if m_index == 1:
        dem_list = model1(params, growthdict)
    elif m_index == 2:
        dem_list = model2(params, growthdict)
    elif m_index == 3:
        dem_list = model2(params, growthdict)

    return(dem_list)
