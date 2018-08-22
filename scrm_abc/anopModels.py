#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  7 18:21:24 2018

@author: scott
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

    def genDem(self, ix, npops, params, demodict, parlist, Ne):
        """
        """
        # remove mMax and mIso from parlist and params
        if any(["mMax" in p for p in parlist]):
            mMax_ix = [i for i, p in enumerate(parlist) if "mMax" in p]
            mMaxdict = {parlist[i]: params[i] for i in mMax_ix}
            mIso_ix = [i for i, p in enumerate(parlist) if "mIso" in p]
            mIsodict = {parlist[i]: params[i] for i in mIso_ix}
            rmv = min(mMax_ix + mIso_ix)
            mMax = True
            parlist = parlist[0:rmv]
            params = params[0:rmv]
        else:
            # no gradual isolationo
            mMax = False
        timedict = defaultdict(list)
        for i, event in enumerate(parlist):
            if 'a' in event:
                # a here means append to previous time
                if params[i] <= 1:
                    if type(params[i-1]) is list:
                        a = params[i] * params[i-1][0]
                        timedict[a].append([event[1:]])
                    else:
                        a = params[i] * params[i-1]
                        timedict[a].append([event[1:]])
                else:
                    if type(params[i-1]) is list:
                        a = params[i] + params[i-1][0]
                        timedict[a].append([event[1:]])
                    else:
                        a = params[i] + params[i-1]
                        timedict[a].append([event[1:]])
            else:
                if type(params[i]) is list:
                    timedict[params[i][0]].append([event])
                    for p in params[i][1:]:
                        timedict[params[i][0]][0].append(p)
                else:
                    timedict[params[i]].append([event])
        if mMax:
            migdict = self.gradIso(timedict, Ne, mMaxdict, mIsodict)
            dd = defaultdict(list)
            # sort dict by times
            for d in (timedict, demodict, migdict):
                for key, value in d.items():
                    dd[key].extend(value)
            od = OrderedDict(sorted(dd.items()))
            # build models
            dem_list = self.modelMig_scrm(od, Ne, npops)
        else:
            dd = defaultdict(list)
            for d in (timedict, demodict):
                for key, value in d.items():
                    dd[key].extend(value)
            od = OrderedDict(sorted(dd.items()))
            # build models
            dem_list = self.model_scrm(od, Ne, npops)
        return(dem_list)

    def gradIso(self, timedict, Ne, mMaxdict, mIsodict, t_ints=10):
        """
        """
        migdict = defaultdict(list)
        # calculate these in coal time, then transform into gens
        if any(["X" in i for i in mMaxdict.keys()]):
            # single value for all
            mMax = list(mMaxdict.values())[0]
            mIso = list(mIsodict.values())[0]
            NemMax = mMax * 4*Ne
            tsc = mIso/(4*Ne)
            for td in list(timedict.keys()):
                if 'ej' in timedict[td][0][0] or 'es' in timedict[td][0][0]:
                    # ts = mIso_list  # list of different speciation/isolation times
                    tdc = td/(4*Ne)
                    if tdc-tsc < 0:
                        tlin = np.linspace(0, tdc, t_ints)
                        mlist = [(NemMax/(tdc+tsc))*(t+tsc) for t in tlin]
                    else:
                        tlin = np.linspace(tdc-tsc, tdc, t_ints)
                        # tlist = [np.round(t) for t in tlin]
                        mlist = [(NemMax/(tdc-(tdc-tsc)))*(t-(tdc-tsc)) for t in tlin]
                    for t, m in zip(tlin, mlist):
                        migdict[np.round(t*4*Ne)].append(['mg' + timedict[td][0][0][2:], m/(4*Ne)])
                else:
                    pass
        else:
            for mi in mMaxdict.keys():
                sp = mi[4:]
                mMax = mMaxdict[mi]
                mIso = mIsodict["mIso{}".format(sp)]
                ej = "ej{}".format(sp)
                es = "es{}".format(sp)
                for td in list(timedict.keys()):
                    if ej == timedict[td][0][0]:
                        NemMax = mMax * 4*Ne
                        tsc = mIso/(4*Ne)
                        tdc = td/(4*Ne)
                        if tdc-tsc < 0:
                            tlin = np.linspace(0, tdc, t_ints)
                            mlist = [(NemMax/(tdc+tsc))*(t+tsc) for t in tlin]
                        else:
                            tlin = np.linspace(tdc-tsc, tdc, t_ints)
                            # tlist = [np.round(t) for t in tlin]
                            mlist = [(NemMax/(tdc-(tdc-tsc)))*(t-(tdc-tsc)) for t in tlin]
                        for t, m in zip(tlin, mlist):
                            migdict[np.round(t*4*Ne)].append(['mg' + timedict[td][0][0][2:], m/(4*Ne)])
                    elif es == timedict[td][0][0]:
                        sp_pair = ["{}".format(es[2:-1]),"{}{}".format(es[2], es[-1])]
                        for i, mmax in enumerate(mMax):
                            NemMax = mmax * 4*Ne
                            tsc = mIso[i]/(4*Ne)
                            tdc = td/(4*Ne)
                            if tdc-tsc < 0:
                                tlin = np.linspace(0, tdc, t_ints)
                                mlist = [(NemMax/(tdc+tsc))*(t+tsc) for t in tlin]
                            else:
                                tlin = np.linspace(tdc-tsc, tdc, t_ints)
                                # tlist = [np.round(t) for t in tlin]
                                mlist = [(NemMax/(tdc-(tdc-tsc)))*(t-(tdc-tsc)) for t in tlin]
                            for t, m in zip(tlin, mlist):
                                migdict[np.round(t*4*Ne)].append(['mg' + sp_pair[i], m/(4*Ne)])
                    else:
                        pass
        return(migdict)

    def modelMig_scrm(self, od, Ne, npops):
        """
        """
        d = False
        dem_list = []
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
                    if len(v[0]) > 4:
                        # migration rates after hybridization
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], v[1]*4*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], v[1]*4*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][4], v[1]*4*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][4], v[0][2], v[1]*4*Ne))
                    else:
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            # symmetric since 1 rate
                            dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], v[1]*4*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], v[1]*4*Ne))
                elif "em" in v[0]:
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            if len(v) > 2:
                                # asymm since 2 rates
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], v[1]*4*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], v[2]*4*Ne))
                            elif len(v) == 2:
                                # symmetric since 1 rate
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], v[1]*4*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], v[1]*4*Ne))
                            else:
                                # stop migration
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], 0))
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], 0))
                elif "es" in v[0]:
                    if any(i in sourcelist for i in v[0]):
                        raise Exception("something wrong in es")
                    else:
                        npops += 1
                        dem_list.append("-es {} {} {}".format(k/(4*Ne), v[0][2], v[1]))
                        dem_list.append("-ej {} {} {}".format(k/(4*Ne), v[0][2], v[0][3]))
                        dem_list.append("-ej {} {} {}".format(k/(4*Ne), npops, v[0][4]))
                        dem_list.append("-eg {} {} {}".format(k/(4*Ne), v[0][2], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][4], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][4], v[0][2], 0))
                        sourcelist.append(v[0][2])
                elif "tm" in v[0]:
                    ev = k/(4*Ne)
                    evS = v[0][2:]
                    pm = v[1]
                    d = True
        if d:
            dem_list.append("-ev {} {} {} {}".format(ev, evS[0], evS[1], pm))
        return(dem_list)

    def model_scrm(self, od, Ne, npops):
        """
        """
        d = False
        dem_list = []
        sourcelist = []
        for k, event in od.items():
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
                elif "em" in v[0]:
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            if len(v) > 2:
                                # asymm since 2 rates
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], v[1]*4*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], v[2]*4*Ne))
                            elif len(v) == 2:
                                # symmetric since 1 rate
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], v[1]*4*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], v[1]*4*Ne))
                            else:
                                # stop migration
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], 0))
                                dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], 0))
                elif "es" in v[0]:
                    if any(i in sourcelist for i in v[0]):
                        raise Exception("something wrong in es")
                    else:
                        npops += 1
                        dem_list.append("-es {} {} {}".format(k/(4*Ne), v[0][2], v[1]))
                        dem_list.append("-ej {} {} {}".format(k/(4*Ne), v[0][2], v[0][3]))
                        dem_list.append("-ej {} {} {}".format(k/(4*Ne), npops, v[0][4]))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][3], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][3], v[0][2], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][2], v[0][4], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(4*Ne), v[0][4], v[0][2], 0))
                        sourcelist.append(v[0][2])
                elif "tm" in v[0]:
                    ev = k/(4*Ne)
                    evS = v[0][2:]
                    pm = v[1]
                    d = True
        if d:
            dem_list.append("-ev {} {} {} {}".format(ev, evS[0], evS[1], pm))
        return(dem_list)
