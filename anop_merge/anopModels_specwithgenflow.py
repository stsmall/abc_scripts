#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:31:12 2020
@author: Scott T. Small

class to build models for filet_sims.py

TODO: better documentation
"""
import numpy as np
from collections import OrderedDict
from collections import defaultdict


class Model(object):
    def __init__(self):
        """
        """
        return(None)

    def genDem(self, npops, params, demodict, parlist, Ne, ploidy):
        """
        """
        # speciation with gradual isolation
        if any(["mMax" in p for p in parlist]):
            # remove mMax and mIso from parlist and params
            mMax_ix = [i for i, p in enumerate(parlist) if "mMax" in p]
            mMaxdict = {parlist[i]: params[i] for i in mMax_ix}
            mIso_ix = [i for i, p in enumerate(parlist) if "mIso" in p]
            mIsodict = {parlist[i]: params[i] for i in mIso_ix}
            rmv = min(mMax_ix + mIso_ix)
            mMax = True
            parlist = parlist[0:rmv]
            params = params[0:rmv]
        # no gradual isolation
        else:
            mMax = False
        # parse events from demo file
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
                    l = []
                    l.append(event)
                    # timedict[params[i][0]].append([event])
                    for p in params[i][1:]:
                        # timedict[params[i][0]][0].append(p)
                        l.append(p)
                    timedict[params[i][0]].append(l)
                else:
                    timedict[params[i]].append([event])
        if mMax:  # run gradual isolation
            migdict = self.gradIso(timedict, Ne, mMaxdict, mIsodict, ploidy)
            dd = defaultdict(list)
            # sort dict by times
            for d in (timedict, demodict, migdict):
                for key, value in d.items():
                    dd[key].extend(value)
            od = OrderedDict(sorted(dd.items()))
            # build models
            dem_list = self.modelMig(od, Ne, npops, ploidy)
        else:  # strict speciation
            dd = defaultdict(list)
            for d in (timedict, demodict):
                for key, value in d.items():
                    dd[key].extend(value)
            od = OrderedDict(sorted(dd.items()))
            # build models
            dem_list = self.modelNoMig(od, Ne, npops, ploidy)
        return(dem_list)

    def gradIso(self, timedict, Ne, mMaxdict, mIsodict, ploidy, t_ints=10):
        """
        """
        migdict = defaultdict(list)
        # calculate these in coal time, then transform into gens
        if any(["X" in i for i in mMaxdict.keys()]):
            # single value for all
            mMax = list(mMaxdict.values())[0]
            mIso = list(mIsodict.values())[0]
            NemMax = mMax * ploidy * Ne
            tsc = mIso/(ploidy*Ne)
            for td in list(timedict.keys()):
                if 'ej' in timedict[td][0][0] or 'es' in timedict[td][0][0]:
                    # ts = mIso_list  # list of different speciation/isolation times
                    tdc = td/(ploidy*Ne)
                    if tdc-tsc < 0:
                        tlin = np.linspace(0, tdc, t_ints)
                        mlist = [(NemMax/(tdc+tsc))*(t+tsc) for t in tlin]
                    else:
                        tlin = np.linspace(tdc-tsc, tdc, t_ints)
                        # tlist = [np.round(t) for t in tlin]
                        mlist = [(NemMax/(tdc-(tdc-tsc)))*(t-(tdc-tsc)) for t in tlin]
                    for t, m in zip(tlin, mlist):
                        migdict[np.round(t*ploidy*Ne)].append(['mg' + timedict[td][0][0][2:], m/(ploidy*Ne)])
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
                        NemMax = mMax * ploidy*Ne
                        tsc = mIso/(ploidy*Ne)
                        tdc = td/(ploidy*Ne)
                        if tdc-tsc < 0:
                            tlin = np.linspace(0, tdc, t_ints)
                            mlist = [(NemMax/(tdc+tsc))*(t+tsc) for t in tlin]
                        else:
                            tlin = np.linspace(tdc-tsc, tdc, t_ints)
                            # tlist = [np.round(t) for t in tlin]
                            mlist = [(NemMax/(tdc-(tdc-tsc)))*(t-(tdc-tsc)) for t in tlin]
                        for t, m in zip(tlin, mlist):
                            migdict[np.round(t*ploidy*Ne)].append(['mg' + timedict[td][0][0][2:], m/(ploidy*Ne)])
                    elif es == timedict[td][0][0]:
                        sp_pair = ["{}".format(es[2:-1]),"{}{}".format(es[2], es[-1])]
                        for i, mmax in enumerate(mMax):
                            NemMax = mmax * ploidy * Ne
                            tsc = mIso[i]/(ploidy*Ne)
                            tdc = td/(ploidy*Ne)
                            if tdc-tsc < 0:
                                tlin = np.linspace(0, tdc, t_ints)
                                mlist = [(NemMax/(tdc+tsc))*(t+tsc) for t in tlin]
                            else:
                                tlin = np.linspace(tdc-tsc, tdc, t_ints)
                                # tlist = [np.round(t) for t in tlin]
                                mlist = [(NemMax/(tdc-(tdc-tsc)))*(t-(tdc-tsc)) for t in tlin]
                            for t, m in zip(tlin, mlist):
                                migdict[np.round(t*ploidy*Ne)].append(['mg' + sp_pair[i], m/(ploidy*Ne)])
                    else:
                        pass
        return(migdict)

    def modelMig(self, od, Ne, npops, ploidy):
        """
        """
        dem_list = []
        sourcelist = []
        for k, event in od.items():
            # event = [y for x in event for y in x]
            for v in event:
                if "ej" in v[0]:
                    dem_list.append("-ej {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3]))
                    dem_list.append("-eg {} {} {}".format(k/(ploidy*Ne), v[0][2], 0))
                    dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], 0))
                    dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], 0))
                    sourcelist.append(v[0][2])
                elif "Ne" in v[0]:
                    if v[0][2] not in sourcelist:
                        if v[1] == "None":
                            dem_list.append("-eg {} {} {}".format(k/(ploidy*Ne), v[0][2], ploidy*Ne*float(v[2])))
                        else:
                            dem_list.append("-en {} {} {}".format(k/(ploidy*Ne), v[0][2], int(v[1])/Ne))
                            dem_list.append("-eg {} {} {}".format(k/(ploidy*Ne), v[0][2], ploidy*Ne*float(v[2])))
                elif "mg" in v[0]:
                    if len(v[0]) > 4:
                        # migration rates after hybridization
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], v[1]*ploidy*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], v[1]*ploidy*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][4], v[1]*ploidy*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][4], v[0][2], v[1]*ploidy*Ne))
                    else:
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            # symmetric since 1 rate
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], v[1]*ploidy*Ne))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], v[1]*ploidy*Ne))
                elif "em" in v[0]:
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            if len(v) > 2:
                                # asymm since 2 rates
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], v[1]*ploidy*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], v[2]*ploidy*Ne))
                            elif len(v) == 2:
                                # symmetric since 1 rate
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], v[1]*ploidy*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], v[1]*ploidy*Ne))
                            else:
                                # stop migration
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], 0))
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], 0))
                elif "es" in v[0]:
                    if any(i in sourcelist for i in v[0]):
                        raise Exception("something wrong in es")
                    else:
                        npops += 1
                        dem_list.append("-es {} {} {}".format(k/(ploidy*Ne), v[0][2], v[1]))
                        dem_list.append("-ej {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3]))
                        dem_list.append("-ej {} {} {}".format(k/(ploidy*Ne), npops, v[0][4]))
                        dem_list.append("-eg {} {} {}".format(k/(ploidy*Ne), v[0][2], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][4], 0))
                        dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][4], v[0][2], 0))
                        sourcelist.append(v[0][2])
                elif "tm" in v[0]:
                    ev = k/(ploidy*Ne)
                    evS = v[0][2:]
                    pm = v[1]
                    dem_list.append("-ev {} {} {} {}".format(ev, evS[0], evS[1], pm))
        return(dem_list)

    def modelNoMig(self, od, Ne, npops, ploidy):
        """
        """
        dem_list = []
        sourcelist = []
        for k, event in od.items():
            for v in event:
                if "ej" in v[0]:
                    dem_list.append("-ej {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3]))
                    dem_list.append("-eg {} {} {}".format(k/(ploidy*Ne), v[0][2], 0))
                    dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], 0))
                    dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], 0))
                    sourcelist.append(v[0][2])
                elif "Ne" in v[0]:
                    if v[0][2] not in sourcelist:
                        if v[1] == "None":
                            dem_list.append("-eg {} {} {}".format(k/(ploidy*Ne), v[0][2], ploidy*Ne*float(v[2])))
                        else:
                            dem_list.append("-en {} {} {}".format(k/(ploidy*Ne), v[0][2], int(v[1])/Ne))
                            dem_list.append("-eg {} {} {}".format(k/(ploidy*Ne), v[0][2], ploidy*Ne*float(v[2])))
                elif "em" in v[0]:
                        if any(i in sourcelist for i in v[0]):
                            pass
                        else:
                            if len(v) > 2:
                                # asymm since 2 rates
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], v[1]*ploidy*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], v[2]*ploidy*Ne))
                            elif len(v) == 2:
                                # symmetric since 1 rate
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], v[1]*ploidy*Ne))
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], v[1]*ploidy*Ne))
                            else:
                                # stop migration
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], 0))
                                dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], 0))
                elif "es" in v[0]:
                    if any(i in sourcelist for i in v[0]):
                        raise Exception("something wrong in es")
                    else:
                        npops += 1
                        try:
                            dem_list.append("-es {} {} {}".format(k/(ploidy*Ne), v[0][2], v[1]))
                            dem_list.append("-ej {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3]))
                            dem_list.append("-ej {} {} {}".format(k/(ploidy*Ne), npops, v[0][4]))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][3], 0))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][3], v[0][2], 0))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][2], v[0][4], 0))
                            dem_list.append("-em {} {} {} {}".format(k/(ploidy*Ne), v[0][4], v[0][2], 0))
                            sourcelist.append(v[0][2])
                        except IndexError:
                            import ipdb;ipdb.set_trace()
                elif "tm" in v[0]:
                    ev = k/(ploidy*Ne)
                    evS = v[0][2:]
                    pm = v[1]
                    dem_list.append("-ev {} {} {} {}".format(ev, evS[0], evS[1], pm))
        return(dem_list)
