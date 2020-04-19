#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 13:31:12 2020.
@author: Scott T. Small.

Main script for generating simulations for training the FILET model of
Schrider et al. 2018.

Example
-------

python filet_sims.py -cfg CONFIG -i 5000 --filet filet_sims/ --ms msmove
    --out noMig.sims.out --ploidy 1.5

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
import configparser
from abc_params import drawParams
from abc_params import get_dist
from abc_models import Model


def simulate(model_dict,
             demo_dict,
             par_gen,
             event_list,
             theta_array,
             rho_array,
             ms_path,
             ploidy,
             rho_mu,
             anc_size,
             tbd):
    """Create a single instance of a call to ms/msmove.

    Parameters defined by priors and distributions in model and config files.

    Parameters
    ----------
    model_dict: Dict
        contains info from config file
    demo_dict: Dict
        ?
    par_gen: generator
        generator of distributions
    par_list: list
        list of demographic events, column 2 in model
    theta_array: array
        array of theta values
    rho_array: array
        array of rho values
    ms_path: str
        location of ms exe
    ploidy: int
        ploidy of the chromosome
    rho_mu: float
        ratio of rho to mu
    effective_size: int
        effective population size

    Returns
    -------
    mscmd: str
        full call to ms/msmove
    params: ??
        ??

    """
    npops = len(model_dict["sampleSize"])
    subpops = "-I {} {}".format(npops, ' '.join(map(str, (model_dict["sampleSize"]))))

    # effective pop size
    theta_n = np.random.choice(theta_array)
    theta_l = theta_n * model_dict["contig_length"]
    scaled_Ne = int(np.round((theta_n/(2*ploidy*model_dict["mutation_rate"]))))

    # recombination rate
    if rho_mu:
        rho = theta_l * rho_mu
    else:
        rho = np.random.choice(rho_array) * model_dict["contig_length"]

    # pop Ne
    ne_sub_pops = ["-n {} {}".format(i+1, (pop_ne/anc_size)) for i, pop_ne in enumerate(model_dict["initialSize"])]
    grow_sub_pops = ["-g {} {}".format(i+1, g) for i, g in enumerate(model_dict["growthRate"])]

    # migration
    miglist = model_dict["migMat"].tolist()
    migmat = [val for sublist in miglist for val in sublist]

    # random draws for parameters
    params_list = get_dist(tbd, par_gen)
    params_list.reverse()
    # build model
    model = Model()
    dem_events = model.genDem(npops,
                              params_list,
                              demo_dict,
                              event_list,
                              scaled_Ne,
                              anc_size,
                              ploidy*2)
    # all to dict
    ms_params = {
                'ms': ms_path,
                'nhaps': sum(model_dict["sampleSize"]),
                'loci': model_dict["loci"],
                'theta': theta_l,
                'rho': rho,
                'basepairs': model_dict["contig_length"] + 1,
                'subpops': subpops,
                'ne_subpop': " ".join(ne_sub_pops),
                'growth_subpop': " ".join(grow_sub_pops),
                'migmat': " ".join(map(str, migmat)),
                'demo': " ".join(dem_events)
                 }
    # ms/msmove command line
    ms_base = ("{ms} {nhaps} {loci} -t {theta} -r {rho} {basepairs} "
               "{subpops} {ne_subpop} {growth_subpop} -ma {migmat} {demo}")
    mscmd = ms_base.format(**ms_params)
    return(mscmd, params_list)


def write_priors(par_list, params_list, header_only):
    """Write priors to file."""
    if header_only:
        header_list = []
        for event in par_list:
            if "tm" in event:
                pop1 = event[-2]
                pop2 = event[-1]
                header_list.append(f"time_ev_{pop1}-{pop2}\tev_prop\t")
            elif "es" in event:
                pop1 = event[-3]
                pop2 = event[-2]
                pop3 = event[-1]
                header_list.append(f"time_es_{pop1}-{pop2}-{pop3}\tes_prop\t")
            elif "ej" in event:
                pop1 = event[-2]
                pop2 = event[-1]
                header_list.append(f"time_ej_{pop1}-{pop2}\t")
            elif "Ne" in event:
                pop1 = event[-1]
                header_list.append(f"Ne_time_{pop1}\tNe_size_{pop1}\tNe_grow_{pop1}\t")
        return(header_list)
    else:
        event_list = []
        for event in params_list:
            if type(event) is list:
                for ev in event:
                    event_list.append(ev)
            else:
                event_list.append(event)
        return(map(str, event_list))


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
    parser.add_argument("--priors", type=str,
                        help="input a list of priors")
    parser.add_argument("--discoal", action="store_true",
                        help="use discoal w/ 0 indexing")
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
    priors_file = args.priors
    # =========================================================================
    #  Config parser
    # =========================================================================
    config = configparser.ConfigParser()
    config.read(configFile)
    #
    sim = "simulation"
    contig_len = config.getint(sim, "contiglen")
    num_loci = config.getint(sim, "loci")
    recomb_rate = config.getfloat(sim, "recombination_rate")
    mutation_rate = config.getfloat(sim, "mutation_rate")
    ancestral_size = config.getint(sim, "effective_population_size")
    #
    init = "initialize"
    sample_sizes = list(map(int, config.get(init, "sample_sizes").split(",")))
    initial_sizes = list(map(int, config.get(init, "initial_sizes").split(",")))
    growth_rate = list(map(float, config.get(init, "growth_rates").split(",")))
    #
    mig_file = config.get(init, "migration_matrix")
    mig_file_path = os.path.join(config_path, mig_file)
    migration_matrix = np.genfromtxt(mig_file_path, delimiter=",")
    assert len(sample_sizes) == migration_matrix.shape[0]
    #
    par = "parameters"
    theta_file = config.get(par, "theta_distribution")
    theta_file_path = os.path.join(config_path, theta_file)
    theta_array = np.loadtxt(theta_file_path)
    rho_file = config.get(par, "rho_distribution")
    rho_file_path = os.path.join(config_path, rho_file)
    rho_array = np.loadtxt(rho_file_path)
    model_file = config.get(par, "params")
    # =========================================================================
    #  Main executions
    # =========================================================================
    model_dict = {"contig_length": contig_len,
                  "recombination_rate": recomb_rate,
                  "mutation_rate": mutation_rate,
                  "sampleSize": sample_sizes,
                  "initialSize": initial_sizes,
                  "growthRate": growth_rate,
                  "migMat": migration_matrix,
                  "loci": num_loci
                  }

    model_file_path = os.path.join(config_path, model_file)
    tbd, par_gen, par_list, demo_dict = drawParams(model_file_path)
    tbd.reverse()
    sim_path = os.path.join(config_path, f"{outfile}.{sim_number}.sims.out")
    priors_outfile = open(f"{sim_path}.priors", 'w')
    priors_list = write_priors(par_list, [], header_only=True)
    priors_outfile.write("{}\n".format("\t".join(priors_list)))
    with open(sim_path, 'w') as f:
        for i in range(sim_number):
            mscmd, params = simulate(model_dict,
                                     demo_dict,
                                     par_gen,
                                     par_list,
                                     theta_array,
                                     rho_array,
                                     ms_path,
                                     ploidy,
                                     rho_mu,
                                     ancestral_size,
                                     tbd)
            priors_list = write_priors(par_list, params, header_only=False)
            priors_outfile.write("{}\n".format("\t".join(priors_list)))
            f.write("{} >> {}\n".format(mscmd, outfile))
