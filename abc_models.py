#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:31:12 2020
@author: Scott T. Small

class to build models for filet_sims.py

TODO: better documentation
"""
from collections import OrderedDict
from collections import defaultdict


class Model(object):
    """Creates command line for coalescent sims."""

    def __init__(self):
        """Create command line for coalescent sims."""
        return(None)

    def genDem(self, npops, params_list, demo_dict, event_list, scaled_Ne, anc_size, ploidy):
        """Parse and order events."""
        # parse events from demo file
        time_dict = defaultdict(list)
        for ix, event in enumerate(event_list):
            if type(params_list[ix]) is list:
                tmp_list = []
                tmp_list.append(event)
                for p in params_list[ix][1:]:
                    tmp_list.append(p)
                time_dict[params_list[ix][0]].append(tmp_list)
            else:
                time_dict[params_list[ix]].append([event])
        # write out model
        dd = defaultdict(list)
        for d in (time_dict, demo_dict):
            for key, value in d.items():
                dd[key].extend(value)
        od = OrderedDict(sorted(dd.items()))
        # build models
        dem_list = self.modelNoMig(od, scaled_Ne, npops, anc_size, ploidy)
        # TODO: add discoal parser
        # if discoal:
        #     dem_list = self.modelNoMig_discoal(od, scaled_Ne, npops, anc_size, ploidy)
        return(dem_list)

    def modelNoMig(self, od, scaled_Ne, npops, anc_size, ploidy):
        """Create model with no migration for msmove, em flag in ms."""
        dem_list = []
        sourcelist = []
        for time, event in od.items():
            for param in event:
                key = param[0]
                new_time = time / (ploidy*scaled_Ne)
                if "ej" in key:
                    # key = ej12
                    # pop1 -> pop2
                    pop1 = key[-2]
                    pop2 = key[-1]
                    dem_list.append(f"-ej {new_time} {pop1} {pop2}")
                    sourcelist.append(pop1)
                elif "Ne" in key:
                    if "_" in key:
                        # key = Ne_1_100000_0
                        _, pop, size, grow = key.split("_")
                    else:
                        pop = key[-1]
                        size = param[1]
                        grow = param[2]
                    if pop not in sourcelist:
                        new_Ne = int(size) / anc_size
                        dem_list.append(f"-en {new_time} {pop} {new_Ne}")
                        if float(grow) > 0:
                            dem_list.append(f"-eg {new_time} {pop} {grow}")
                elif "es" in key:
                    # es345
                    if any(i in sourcelist for i in key[2:]):
                        raise Exception("something wrong in es")
                    else:
                        pop1 = key[-3]
                        pop2 = key[-2]
                        pop3 = key[-1]
                        admx = param[1]
                        # admx is probability that each lineage stays in pop-i,1-p are admx
                        npops += 1
                        try:
                            dem_list.append(f"-es {new_time} {pop1} {admx}")
                            dem_list.append(f"-ej {new_time} {pop1} {pop2}")
                            dem_list.append(f"-ej {new_time} {npops} {pop3}")
                            sourcelist.append(pop1)
                        except IndexError:
                            breakpoint()
                elif "tm" in key:
                    # key = tm12
                    # -ev time pop_i pop_j prob_x; move i into j w/ prob x
                    # forward is j -> i
                    try:
                        pop1 = key[-2]
                        pop2 = key[-1]
                        prob = param[1]
                    except IndexError:
                        breakpoint
                    dem_list.append(f"-ev {new_time} {pop1} {pop2} {prob}")
        return(dem_list)

    def modelNoMig_discoal(self, od, scaled_Ne, npops, anc_size, ploidy):
        """Create model with no migration for discoal, em flag in ms."""
        dem_list = []
        sourcelist = []
        for time, event in od.items():
            for param in event:
                key = param[0]
                new_time = time / (ploidy*scaled_Ne)
                if "ej" in key:
                    # key = ej12
                    # pop1 -> pop2
                    pop1 = key[-2]
                    pop2 = key[-1]
                    dem_list.append(f"-ej {new_time} {pop1} {pop2}")
                    sourcelist.append(pop1)
                elif "Ne" in key:
                    # key = Ne_1_100000_0
                    _, pop, size, grow = key.split("_")
                    if pop not in sourcelist:
                        new_Ne = int(size) / anc_size
                        dem_list.append(f"-en {new_time} {pop} {new_Ne}")
                elif "es" in key:
                    # es345
                    if any(i in sourcelist for i in key[2:]):
                        raise Exception("something wrong in es")
                    else:
                        pop1 = key[-3]
                        pop2 = key[-2]
                        pop3 = key[-1]
                        admx = param[1]
                        # admx is probability that each lineage stays in pop-i,1-p are admx
                        npops += 1
                        try:
                            dem_list.append(f"-es {new_time} {pop1} {admx}")
                            dem_list.append(f"-ej {new_time} {pop1} {pop2}")
                            dem_list.append(f"-ej {new_time} {npops} {pop3}")
                            sourcelist.append(pop1)
                        except IndexError:
                            breakpoint()
                elif "tm" in key:
                    # key = tm12
                    # -ev time pop_i pop_j prob_x; move i into j w/ prob x
                    # forward is j -> i
                    try:
                        pop1 = key[-2]
                        pop2 = key[-1]
                        prob = param[1]
                    except IndexError:
                        breakpoint
                    dem_list.append(f"-ev {new_time} {pop1} {pop2} {prob}")
        return(dem_list)
