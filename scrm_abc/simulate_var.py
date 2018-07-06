#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 14:12:32 2018

@author: scott
"""

#def simulate_var(model, demodict, ix, parfx, parlist, thetaarray, rhoarray, its):
#    """Runs scrm by building it[eration]s number of files. Each line is then
#    executed in the shell and output is printed to the same file with
#    parameter choices in the file name. The loci parameter controls the size of
#    the fragment.
#    Example: its = 10,000, this will produce 10,000 input files each with
#    lines == loci. Each line is a separate execution of scrm with a different
#    value of Ne and Rho. This will produce 1 output files per input file.
#    """
#
#    params = [p() for p in parfx]
#    # avoid exactly equal times by small perturbation
#    while len(params) != len(set(params)):
#        params = [i + np.random.randint(-10, 10) for i in params]
#    # simulate
#    for i in range(model["loci"]):
#        Ne = int(np.round((np.random.choice(thetaarray)/(4*model["mutation_rate"]))))
#        theta = 4*Ne*model["mutation_rate"] * model["contig_length"]
#        rho = np.random.choice(rhoarray) * model["contig_length"]
#        subpops = "-I {} {}".format(len(model["sampleSize"]), ' '.join(map(str, (model["sampleSize"]))))
#        nesubpops = ["-n {} {}".format(i+1, ne/Ne) for i, ne in enumerate(model["initialSize"])]
#        gsubpops = ["-g {} {}".format(i+1, g*4*Ne) for i, g in enumerate(model["growthRate"])]
#        miglist = model["migMat"].tolist()
#        migmat = [val for sublist in miglist for val in sublist]
#        m = Model()
#        dem_events = m.genDem(ix, params, demodict, parlist, Ne, model["mMax"], model["mIso"])
#        # all to dict
#        ms_params = {
#                    'nhaps': sum(model["sampleSize"]),
#                    'loci': 1,
#                    'theta': theta,
#                    'rho': rho,
#                    'basepairs': model["contig_length"] + 1,
#                    'subpops': subpops,
#                    'ne_subpop': " ".join(nesubpops),
#                    'growth_subpop': " ".join(gsubpops),
#                    'migmat': " ".join(map(str, migmat)),
#                    'demo': " ".join(dem_events)
#                     }
#        # scrm command line
#        scrm_base = ("scrm {nhaps} {loci} -t {theta} -r {rho} {basepairs} "
#                     "{subpops} {ne_subpop} {growth_subpop} -ma {migmat} {demo} "
#                     "-SC abs -p 12")
#        mscmd = scrm_base.format(**ms_params)
#    return(mscmd, params)