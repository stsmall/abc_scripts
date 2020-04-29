#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:31:12 2020
@author: Scott T. Small
class to build models for abc_sims.py
"""
from collections import OrderedDict


class Model(object):
    """Creates command line for coalescent sims."""

    def __init__(self):
        """Create command line for coalescent sims.

        Returns
        -------
        None.

        """
        return(None)

    def genDem(self, ms_dict, model_dict, demo_dict, event_dict, ms_path):
        """Parse and order events.

        Parameters
        ----------
        ms_dict : TYPE
            DESCRIPTION.
        model_dict : TYPE
            DESCRIPTION.
        demo_dict : TYPE
            DESCRIPTION.
        event_dict : TYPE
            DESCRIPTION.
        ms_path : TYPE
            DESCRIPTION.

        Returns
        -------
        dem_list : TYPE
            DESCRIPTION.

        """
        for event in event_dict.keys():
            time, *events = event_dict[event][0]
            demo_dict[time].append([event, events])
        ord_events = OrderedDict(sorted(demo_dict.items()))

        # build models
        if "discoal" in ms_path:
            dem_list = self.modelNoMig_discoal(ord_events, model_dict, ms_dict)
        elif "msprime" in ms_path:
            pass
            # dem_list = self.modelNoMig_msprime(ord_events, model_dict, ms_dict)
        else:
            dem_list = self.modelNoMig(ord_events, model_dict, ms_dict)
        return dem_list

    def modelNoMig(self, ord_events, model_dict, ms_dict):
        """Create model with no migration for ms/msmove.

        Parameters
        ----------
        ord_events : TYPE
            DESCRIPTION.
        model_dict : TYPE
            DESCRIPTION.
        ms_dict : TYPE
            DESCRIPTION.

        Returns
        -------
        dem_list : TYPE
            DESCRIPTION.

        """
        ploidy = model_dict["ploidy"] * 2
        npops = ms_dict["npops"]
        scaled_Ne = ms_dict["scaled_Ne"]
        init_size = list(model_dict["initialSize"])
        dem_list = []
        sourcelist = []
        for time, eve in ord_events.items():
            new_time = time / (ploidy*scaled_Ne)
            for keys in eve:
                if "Ne" in keys:
                    # key = Ne_1_100000_0
                    _, pop, size, grow = keys.split("_")
                    init_size[int(pop)-1] = size
                    if pop not in sourcelist:
                        new_Ne = int(size) / scaled_Ne
                        if float(grow) > 0:
                            dem_list.append(f"-en {new_time} {pop} {new_Ne} {grow}")
                        else:
                            dem_list.append(f"-en {new_time} {pop} {new_Ne}")
                else:
                    key = keys[0]
                    param = keys[1]
                    if "ej" in key:
                        # ej_12
                        pop1, pop2 = key.split("_")[1]
                        # pop1 -> pop2
                        if pop1 not in sourcelist:
                            dem_list.append(f"-ej {new_time} {pop1} {pop2}")
                            sourcelist.append(pop1)
                    elif "es" in key:
                        # es_34; es in ms/msmove
                        pop1, pop2 = key.split("_")[1]
                        if pop1 not in sourcelist:
                            prop = param[0]  # 1-prop are admixed from pop2
                            npops += 1
                            dem_list.append(f"-es {new_time} {pop1} {prop}")
                            dem_list.append(f"-ej {new_time} {npops} {pop2}")
                            sourcelist.append(npops)
                    elif "tm" in key:
                        # tm_12; ev in msmove
                        pop1, pop2 = key.split("_")[1]
                        if pop1 not in sourcelist:
                            prob = param[0]
                            dem_list.append(f"-ev {new_time} {pop1} {pop2} {prob}")
                    elif "em" in key:
                        pop1, pop2 = key.split("_")[1]
                        if not any(i in sourcelist for i in [pop1, pop2]):
                            mig = param[0]*ploidy*init_size[pop1]
                            dem_list.append(f"-em {new_time} {pop1} {pop2} {mig}")
        return dem_list

    def modelNoMig_discoal(self,  ord_events, model_dict, ms_dict):
        """Create model with no migration for discoal.

        Parameters
        ----------
        ord_events : TYPE
            DESCRIPTION.
        model_dict : TYPE
            DESCRIPTION.
        ms_dict : TYPE
            DESCRIPTION.

        Returns
        -------
        dem_list : TYPE
            DESCRIPTION.

        """
        ploidy = model_dict["ploidy"] * 2
        npops = ms_dict["npops"]
        scaled_Ne = ms_dict["scaled_Ne"]
        init_size = list(model_dict["initialSize"])
        dem_list = []
        sourcelist = []
        for time, eve in ord_events.items():
            new_time = time / (ploidy*scaled_Ne)
            for keys in eve:
                if "Ne" in keys:
                    # key = Ne_1_100000_0
                    _, pop, size, grow = keys.split("_")
                    pop = f"{int(pop)-1}"
                    init_size[int(pop)] = size
                    if pop not in sourcelist:
                        new_Ne = int(size) / scaled_Ne
                        if float(grow) > 0:
                            dem_list.append(f"-en {new_time} {pop} {new_Ne} {grow}")
                        else:
                            dem_list.append(f"-en {new_time} {pop} {new_Ne}")
                else:
                    key = keys[0]
                    param = keys[1]
                    if "ej" in key:
                        # ed_12
                        pop1, pop2 = key.split("_")[1]
                        pop1 = f"{int(pop1)-1}"
                        pop2 = f"{int(pop2)-1}"
                        if pop1 not in sourcelist:
                            dem_list.append(f"-ed {new_time} {pop1} {pop2}")
                            sourcelist.append(pop1)
                    elif "es" in key:
                        # es_343; daughter, parent1, parent2. can be same
                        pop1, pop2, pop3 = key.split("_")[1]
                        pop1 = f"{int(pop1)-1}"
                        pop2 = f"{int(pop2)-1}"
                        pop3 = f"{int(pop3)-1}"
                        if not any(i in sourcelist for i in [pop1, pop2, pop3]):
                            prop = 1 - param[0]
                            dem_list.append(f"-ea {new_time} {pop1} {pop2} {pop3} {prop}")
                    elif "em" in key:
                        pop1, pop2 = key.split("_")[1]
                        pop1 = f"{int(pop1)-1}"
                        pop2 = f"{int(pop2)-1}"
                        if not any(i in sourcelist for i in [pop1, pop2]):
                            mig = param[0]*ploidy*init_size[pop1]
                            dem_list.append(f"-m {new_time} {pop1} {pop2} {mig}")

        return dem_list
