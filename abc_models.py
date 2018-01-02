#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 15:44:48 2017

@author: scott
"""
import msprime as msp
import numpy as np
from math import log


def model1(configdict):
    """Model with no migration
    """
    generations = 1
    pops = configdict["initial_size"]
    np.random.seed()
    # anc population sizes in Ne
    N_HM = np.random.randint(100, 10000)
    N_HMK = np.random.randint(100, 10000)
    N_HMKP = np.random.randint(100, 10000)
    N_ANC = np.random.randint(100, 500000)
    # join times in gen
    T_HM = np.random.randint(100, 1000) / generations
    T_HMK = np.random.randint(100, 10000) / generations
    T_HMKP = np.random.randint(1000, 40000) / generations
    T_ANC = configdict["tmax"]
    while T_HM > T_HMK:
        T_HM = np.random.randint(100, 1000) / generations
        T_HMK = np.random.randint(100, 10000) / generations
    while T_HMK > T_HMKP:
        T_HMKP = np.random.randint(1000, 40000) / generations
    # intial growth rates
    r_H = -log(N_HM/pops[0]) / T_HM
    r_M = -log(N_HM/pops[1]) / T_HM
    r_K = -log(N_HMK/pops[2]) / T_HMK
    r_P = -log(N_HMKP/pops[3]) / T_HMKP
    # ancestral growth rates
    r_HM = -log(N_HMK/N_HM) / T_HMK
    r_HMK = -log(N_HMKP/N_HMK) / T_HMKP
    r_HMKP = -log(N_ANC/N_HMKP) / T_ANC
    priors = [0, N_HM, N_HMK, N_HMKP, N_ANC, T_HM, T_HMK, T_HMKP, r_H, r_M,
              r_K, r_P, r_HM, r_HMK, r_HMKP]
    dem = [
           msp.PopulationParametersChange(time=0, initial_size=pops[0], growth_rate=r_H, population_id=0),
           msp.PopulationParametersChange(time=0, initial_size=pops[1], growth_rate=r_M, population_id=1),
           msp.PopulationParametersChange(time=0, initial_size=pops[2], growth_rate=r_K, population_id=2),
           msp.PopulationParametersChange(time=0, initial_size=pops[3], growth_rate=r_P, population_id=3),
           msp.MassMigration(time=T_HM, source=0, destination=1, proportion=1.0),
           msp.PopulationParametersChange(time=T_HM, initial_size=N_HM, growth_rate=r_HM, population_id=1),
           msp.MassMigration(time=T_HMK, source=1, destination=2, proportion=1.0),
           msp.PopulationParametersChange(time=T_HMK, initial_size=N_HMK, growth_rate=r_HMK, population_id=2),
           msp.MassMigration(time=T_HMKP, source=2, destination=3, proportion=1.0),
           msp.PopulationParametersChange(time=T_HMKP, initial_size=N_HMKP, growth_rate=r_HMKP, population_id=3),
           msp.PopulationParametersChange(time=T_ANC, initial_size=N_ANC)]
    return(dem, priors, N_ANC)


def model2(configdict):
    """Model with migration from Mali to Haiti
    """
    # Migration rates during the various epochs.
    m_M_H = 0
    m_K_H = 0
    m_M_K = 0
    m_K_M = 0
    m_M_H = np.random.beta()  # (0, 1)
    migration_matrix = [
                        [0, m_M_H, m_K_H, 0],
                        [0, 0, m_K_M, 0],
                        [0, m_M_K, 0, 0],
                        [0, 0, 0, 0]
                        ]


def model3(configdict):
    """Model with migration from Mali to Haiti and Kenya to Haiti
    """
    # Migration rates during the various epochs.
    m_M_H = 0
    m_K_H = 0
    m_M_K = 0
    m_K_M = 0
    m_M_H = np.random.beta()  # (0, 1)
    m_K_H = np.random.beta()  # (0, 2)
    migration_matrix = [
                        [0, m_M_H, m_K_H, 0],
                        [0, 0, m_K_M, 0],
                        [0, m_M_K, 0, 0],
                        [0, 0, 0, 0]
                        ]
