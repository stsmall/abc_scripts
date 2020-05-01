#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:31:12 2020.
@author: Scott T. Small.

Main script for generating simulations for ABC and ML training.
depends: sim_models.py, sim_params.py

Example
-------

abc_sims.py -cfg examples/example.cfg -i 100000 --ms msmove --out test

generates a file with a random name that has 5,000 lines. each line is a call
to msmove under the CONFIG specifications. The --filet option tells the main
where to find the other modules. Basically an import issue.

Notes
-----
This relies on a config file and a model file. It can also use a distribution
of pi and rho (the population recombination rate). See the github for examples.


"""
import sys
import argparse
import os
import numpy as np
import re
import configparser
from collections import defaultdict
from tqdm import trange
from sim_params import drawParams
from sim_params import get_dist
from sim_models import Model


def selection_parse(model_dict, ms_dict):
    """Parse selection dict for discoal.

    Parameters
    ----------
    model_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    scaled_Ne = ms_dict["scaled_Ne"]
    rho_loc = ms_dict["rho_loc"],
    # sel params
    sel_dict = model_dict["sel_dict"]
    pop0_Ne = sel_dict["pop0_Ne"]
    if pop0_Ne == 0:
        pop0_Ne = scaled_Ne
    sweep_Ne = sel_dict["sweep_Ne"]
    alpha = sel_dict["alpha"]
    freq = sel_dict["freq"]
    sweep_stop = sel_dict["sweep_stop"]
    sweep_site = sel_dict["sweep_site"]
    part_freq = sel_dict["part_freq"]
    adapt = sel_dict["adapt"]
    hide = sel_dict["hide"]

    sel_list = []

    if sweep_Ne > 0:
        sel_list.append(f"-N {sweep_Ne}")
    # sweep time
    if type(sweep_stop) == list:
        ws_l = sweep_stop[0]/(4 * scaled_Ne)
        ws_h = sweep_stop[1]/(4 * scaled_Ne)
        sel_list.append(f"-ws 0 -Pu {ws_l} {ws_h}")
    else:
        tau = sweep_stop/(4 * scaled_Ne)
        sel_list.append(f"-ws {tau}")
    # sel coeff
    if type(alpha) == list:
        a_low = alpha[0] * 2 * pop0_Ne
        a_high = alpha[1] * 2 * pop0_Ne
        sel_list.append(f"-Pa {a_low} {a_high}")
    else:
        a = alpha * 2 * pop0_Ne
        sel_list.append(f"-a {a}")
    # offscreen sweep
    left_rho = sel_dict["left_rho"]
    if left_rho[1] > 0:
        time, scale = left_rho
        sel_list.append(f"-ls {time/(4*scaled_Ne)} {scale*rho_loc}")
    # recurrent to the left
    rrh_left = sel_dict["rrh_left"]
    if rrh_left > 0:
        sel_list.append(f"-L {rrh_left}")
    # recurrent at locus
    rrh_loc = sel_dict["rrh_loc"]
    if rrh_loc > 0:
        sel_list.append(f"-R {rrh_loc}")
    # starting freq
    if type(freq) == list:
        f_l, f_h = freq
        assert f_l >= 0
        assert f_h <= 1
        assert f_l < f_h
        sel_list.append(f"-Pf {f_l} {f_h}")
    elif freq > 0:
        assert freq <= 1
        sel_list.append(f"-f {freq}")
    else:
        # hard sweep
        pass
    # sweep site
    if type(sweep_site) == list:
        s_l, s_h = sweep_site
        assert s_l >= 0
        assert s_h <= 1
        assert s_l < s_h
        sel_list.append(f"-Px {s_l} {s_h}")
    else:
        assert 0 <= sweep_site <= 1
        sel_list.append(f"-x {sweep_site}")
    # partial sweep freq
    if type(part_freq) == list:
        p_l, p_h = part_freq
        assert p_l >= 0
        assert p_h <= 1
        assert p_l < p_h
        sel_list.append(f"-Pc {p_l} {p_h}")
    elif part_freq > 0:
        assert part_freq <= 1
        sel_list.append(f"-c {part_freq}")
    else:
        # no partial sweep
        pass
    # reccurrent adaptive mut
    if type(adapt) == list:
        a_l, a_h = adapt
        sel_list.append(f"-PuA {a_l} {a_h}")
    elif adapt > 0:
        sel_list.append(f"-uA {adapt}")
    else:
        # no recurrent adaptive mutation
        pass
    if hide:
        sel_list.append("-h")
    return sel_list


def coalSimSyntax(ms_path, model_dict):
    """Create parameters for specific model.

    Parameters
    ----------
    model_dict : TYPE
        DESCRIPTION.
    ms_path : str
        path to simulator
    Returns
    -------
    ms_dict : Dict
        BOO

    """
    ms_dict = {}

    # populations
    samples_sizes = model_dict["sampleSize"]
    npops = len(samples_sizes)
    ploidy = model_dict["ploidy"]
    theta_arr = model_dict["theta"]
    # effective pop size
    mut_rate = model_dict["mutation_rate"]
    if type(mut_rate) == list:
        low, high = mut_rate
        mut_rate = np.random.uniform(low, high)
    locus_len = model_dict["contig_length"]
    theta_nuc = np.random.choice(theta_arr)
    scaled_Ne = int(np.round((theta_nuc/(2*ploidy*mut_rate))))
    theta_loc = theta_nuc * locus_len

    # recombination rate
    rho_arr = model_dict["rho"]
    rho_mu = model_dict["rho_mu"]
    rec_rate = model_dict["recombination_rate"]
    if rec_rate:
        if type(rec_rate) == list:
            low, high = mut_rate
            rec_rate = np.random.uniform(low, high)
        else:
            rho = 4*scaled_Ne*rec_rate
    elif rho_mu:
        rho = theta_loc * rho_mu
    elif any(rho_arr):
        rho = np.random.choice(rho_arr) * locus_len
    else:
        rho = 0

    # gene conversion
    gen_conversion = model_dict["gene_conversion"][0]
    if gen_conversion > 0:
        tract = model_dict["gene_conversion"][1]
        if "discoal" in ms_path:
            gen_cov = f"-gr {gen_conversion*rho} {tract}"
        else:
            gen_cov = f"-c {gen_conversion*rho} {tract}"
    else:
        gen_cov = ""

    # subops
    init_sizes = model_dict["initialSize"]
    grow_rate = model_dict["growthRate"]
    mig_list = model_dict["migMat"].tolist()
    mig_mat = [val for sublist in mig_list for val in sublist]
    if "msprime" in ms_path:
        pass
    elif "discoal" in ms_path:
        subpops = f"-p {npops} {' '.join(map(str, sample_sizes))}"
        ne_sub_pops = [f"-en {i} {pop_ne/scaled_Ne}" for i, pop_ne in enumerate(init_sizes)]
        ne_subpop = " ".join(ne_sub_pops)
        grow_subpop = []
        if sum(mig_mat) == 0:
            mig_matrix = ''
        else:
            mig = []
            mig_matrix = zip(*mig_mat)
            for p, pop_m in enumerate(mig_matrix):
                for i, m in pop_m:
                    if p != i and m > 0:
                        mig.append(f"-m {p} {i} {m}")
    else:
        subpops = f"-I {npops} {' '.join(map(str, sample_sizes))}"
        ne_sub_pops = [f"-n {i+1} {pop_ne/scaled_Ne}" for i, pop_ne in enumerate(init_sizes)]
        ne_subpop = " ".join(ne_sub_pops)
        if sum(grow_rate) == 0:
            grow_subpop = ""
        else:
            grow_sub_pops = [f"-g {i+1} {g}" for i, g in enumerate(grow_rate)]
            grow_subpop = " ".join(grow_sub_pops)
        if sum(mig_mat) == 0:
            mig_matrix = ""
        else:
            mig_matrix = f"-ma {' '.join(map(str, mig_mat))}"

    ms_dict = {"npops": npops,
               "subpop": subpops,
               "theta_loc": theta_loc,
               "scaled_Ne": scaled_Ne,
               "rho_loc": rho,
               "gen_cov": gen_cov,
               "ne_subpop": ne_subpop,
               "grow_subpop": grow_subpop,
               "mig_matrix": mig_matrix}

    return ms_dict


def simulate(model_dict, demo_dict, par_dict, eventkey_dict, cond_list, ms_path):
    """Create a single instance of a call to ms/msmove.

    Parameters
    ----------
    model_dict: Dict
        contains info from config file
    demo_dict: Dict
        dict of pop sizes and times
    par_dict: Dict
        generator of distributions
    event_dict: List
        list of demographic events, column 2 in model
    ms_path: str
        location of ms exe

    Returns
    -------
    mscmd: str
        full call to ms/msmove
    params: list
        list of random parameters for that simulation

    """
    ms_dict = coalSimSyntax(ms_path, model_dict)
    # build selection command line
    if model_dict["sel_dict"]:
        sel_list = selection_parse(model_dict, ms_dict)
    else:
        sel_list = ''

    # draw and organize for parameters
    params_dict = get_dist(par_dict)
    for condit in cond_list:
        lt, gt = condit
        assert params_dict[lt][0] < params_dict[gt][0]

    time_dict = defaultdict(lambda: defaultdict(list))
    for time, event in demo_dict.items():
        time_dict[time]["Ne"].extend(event)

    for event in eventkey_dict.keys():
        for tbi in eventkey_dict[event]:
            if "Ne" in event:
                time, size, grow = params_dict[tbi]
                time_dict[time]["Ne"].append(f"{event}_{size}_{grow}")
            else:
                time, *params = params_dict[tbi]
                if "ej" in event:
                    params = [0]
                time_dict[time][event] = params
    # build demographic command line
    model = Model()
    dem_events = model.genDem(ms_dict, model_dict, time_dict, ms_path)
    # gather command line args
    ms_params = {
                'ms': ms_path,
                'nhaps': sum(model_dict["sampleSize"]),
                'loci': model_dict["loci"],
                'theta': ms_dict["theta_loc"],
                'rho': ms_dict['rho_loc'],
                'gen_cov': ms_dict['gen_cov'],
                'basepairs': model_dict["contig_length"] + 1,
                'subpops': ms_dict["subpop"],
                'ne_subpop': ms_dict["ne_subpop"],
                'growth_subpop': ms_dict["grow_subpop"],
                'migmat': ms_dict["mig_matrix"],
                'demo': " ".join(dem_events),
                'sel': " ".join(sel_list)
                }
    if "discoal" in ms_path:
        ms_base = ("{ms} {nhaps} {loci} {basepairs} -t {theta} -r {rho} "
                   "{gen_cov} {subpops} {ne_subpop} {demo} {sel}")
    elif "msprime" in ms_path:
        pass
    else:
        ms_base = ("{ms} {nhaps} {loci} -t {theta} -r {rho} {basepairs} "
                   "{gen_cov} {subpops} {ne_subpop} {growth_subpop} {migmat} {demo}")

    # ms/msmove/discoal/msprime command line
    mscmd = ms_base.format(**ms_params)
    ms_cmd = " ".join(mscmd.split())
    return(ms_cmd, params_dict)


def write_priors(eventkey_dict, params_dict):
    """Write priors to file."""
    if eventkey_dict:
        header_list = []
        for event in eventkey_dict.keys():
            for i, eve in enumerate(eventkey_dict[event]):
                if "tm" in event:
                    pop1 = event[-2]
                    pop2 = event[-1]
                    header_list.append(f"time_ev_{pop1}-{pop2}\tev_prop")
                elif "es" in event:
                    pop1 = event[-3]
                    pop2 = event[-2]
                    pop3 = event[-1]
                    header_list.append(f"time_es_{pop1}-{pop2}-{pop3}\tes_prop")
                elif "ej" in event:
                    pop1 = event[-2]
                    pop2 = event[-1]
                    header_list.append(f"time_ej_{pop1}-{pop2}")
                elif "Ne" in event:
                    pop1 = event[-1]
                    header_list.append(f"Ne_time_{pop1}_{i}\tNe_size_{pop1}_{i}\tNe_grow_{pop1}_{i}")
        return("\t".join(header_list))
    else:
        params_list = list(params_dict.keys())
        params_list.sort(key=lambda x: int(x[3:]))
        par_list = []
        for event in params_list:
            for params in params_dict[event]:
                if type(params) is list:
                    for param in params:
                        par_list.append(param)
                else:
                    par_list.append(params)
            priors_list = map(str, par_list)
        return("\t".join(priors_list))


def parse_args(args_in):
    """Parse args."""
    parser = argparse.ArgumentParser(prog=sys.argv[0],
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-cfg', "--configFile", required=True,
                        help="path to config file")
    parser.add_argument('-i', "--iterations", type=int, default=1,
                        help="number of iterations, number of lines")
    parser.add_argument("--ms", type=str, required=True,
                        help=" full path to ms/msmove/discoal exe")
    parser.add_argument("--out", type=str, required=True,
                        help="outfilename to write simulations")
    parser.add_argument("--ploidy", type=float, default=2,
                        help="options: hap=1, sex=1.5, auto=2")
    parser.add_argument("--rhomu", type=float,
                        help="ratio of rho/mu")
    return(parser.parse_args(args_in))


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    # =========================================================================
    #  Gather args
    # =========================================================================
    configFile = args.configFile
    config_dir = os.path.abspath(configFile)
    config_path = os.path.split(config_dir)[0]
    outfile = args.out
    sim_number = args.iterations
    ms_path = args.ms
    ploidy = args.ploidy
    rho_mu = args.rhomu

    # =========================================================================
    #  Config parser
    # =========================================================================
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(configFile)
    #
    sim = "simulation"
    contig_len = config.getint(sim, "contiglen")
    num_loci = config.getint(sim, "loci")
    recomb_rate = config.get(sim, "recombination_rate")
    if recomb_rate:
        if "," in recomb_rate:
            recomb_rate = list(map(float, recomb_rate.split(",")))
        else:
            recomb_rate = float(recomb_rate)
    mutation_rate = config.get(sim, "mutation_rate")
    if "," in mutation_rate:
        mutation_rate = list(map(float, mutation_rate.split(",")))
    else:
        mutation_rate = float(mutation_rate)
        assert mutation_rate > 0
    ancestral_size = config.getint(sim, "ancestral_population_size")
    #
    init = "initialize"
    sample_sizes = list(map(int, config.get(init, "sample_sizes").split(",")))
    initial_sizes = list(map(int, config.get(init, "initial_sizes").split(",")))
    growth_rate = list(map(float, config.get(init, "growth_rates").split(",")))
    gene_conversion = list(map(float, config.get(init, "gene_conversion").split(",")))
    #
    mig_file = config.get(init, "migration_matrix")
    mig_file_path = os.path.join(config_path, mig_file)
    migration_matrix = np.genfromtxt(mig_file_path, delimiter=",")
    if np.sum(migration_matrix) > 0:
        assert len(sample_sizes) == migration_matrix.shape[0]
    #
    sel = "selection"
    add_sel = config.getboolean(sel, "sel")
    hide = config.getboolean(sel, "hide")
    if add_sel:
        alpha = config.get(sel, "sel_coeff")
        freq = config.get(sel, "soft_freq")
        sweep_stop = config.get(sel, "sweep_time")
        sweep_site = config.get(sel, "sweep_site")
        part_freq = config.get(sel, "partial_sweep")
        adapt = config.get(sel, "adapt_mutrate")
        left_rho = config.get(sel, "leftRho")
        rrh_left = config.get(sel, "Lrecurrent")
        rrh_loc = config.get(sel, "Rrecurrent")
        pop0_Ne = config.get(sel, "pop1_effective_size")
        sweep_Ne = config.get(sel, "sweep_effective_size")

        sel_dict = {"alpha": alpha,
                    "freq": freq,
                    "sweep_stop": sweep_stop,
                    "sweep_site": sweep_site,
                    "part_freq": part_freq,
                    "adapt": adapt,
                    "left_rho": left_rho,
                    "rrh_left": rrh_left,
                    "rrh_loc": rrh_loc,
                    "pop0_Ne": pop0_Ne,
                    "sweep_Ne": sweep_Ne
                    }

        for key in sel_dict.keys():
            if "," in sel_dict[key]:
                sel_dict[key] = list(map(float, sel_dict[key].split(",")))
            else:
                if "." in sel_dict[key]:
                    sel_dict[key] = float(sel_dict[key])
                else:
                    sel_dict[key] = int(sel_dict[key])
        sel_dict["hide"] = hide
    else:
        sel_dict = {}

    #
    par = "parameters"
    theta_file = config.get(par, "theta_distribution")
    theta_file_path = os.path.join(config_path, theta_file)
    theta_array = np.loadtxt(theta_file_path)
    rho_file = config.get(par, "rho_distribution")
    rho_file_path = os.path.join(config_path, rho_file)
    rho_array = np.loadtxt(rho_file_path)
    model_file = config.get(par, "params")
    #
    model_dict = {"contig_length": contig_len,
                  "recombination_rate": recomb_rate,
                  "mutation_rate": mutation_rate,
                  "sampleSize": sample_sizes,
                  "initialSize": initial_sizes,
                  "growthRate": growth_rate,
                  "gene_conversion": gene_conversion,
                  "migMat": migration_matrix,
                  "loci": num_loci,
                  "anc_size": ancestral_size,
                  "ploidy": ploidy,
                  "rho_mu": rho_mu,
                  "theta": theta_array,
                  "rho": rho_array,
                  "sel_dict": sel_dict
                  }
    # =========================================================================
    #  Set Path and File variable
    # =========================================================================
    model_file_path = os.path.join(config_path, model_file)
    sim_path = os.path.join(config_path, f"{outfile}.{sim_number}.sims.out")

    # =========================================================================
    #  Main executions
    # =========================================================================
    #
    conditions, par_dict, eventkey_dict, demo_dict = drawParams(model_file_path)
    #
    with open(f"{sim_path}.priors", 'w') as priors_outfile:
        priors_list = write_priors(eventkey_dict, [])
        priors_outfile.write(f"{priors_list}\n")
        with open(sim_path, 'w') as sims_outfile:
            for i in trange(sim_number):
                mscmd, params_dict = simulate(model_dict,
                                              demo_dict.copy(),
                                              par_dict,
                                              eventkey_dict,
                                              conditions,
                                              ms_path)
                priors_list = write_priors([], params_dict)
                priors_outfile.write(f"{priors_list}\n")
                sims_outfile.write(f"{mscmd} >> {outfile}\n")
