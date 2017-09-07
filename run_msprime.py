#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 13:58:06 2017
run_msprime.py
Takes model object, reps, etc ...
@author: scott
"""

import argparse
import numpy as np
import msprime as ms


def variant_matrix_example():
    """Outputs msprime into format for scikit-allel
    """
    print("\nCreating full variant matrix")
    tree_sequence = msprime.simulate(
        sample_size=20, Ne=1e4, length=5e3, recombination_rate=2e-8,
        mutation_rate=2e-8, random_seed=10)
    shape = tree_sequence.get_num_mutations(), tree_sequence.get_sample_size()
    A = np.empty(shape, dtype="u1")
    for variant in tree_sequence.variants():
        A[variant.index] = variant.genotypes
    print(A)


tree_sequence = msp.simulate(sample_size=5, Ne=1000, length=1e4,
                                 recombination_rate=2e-8, mutation_rate=2e-8)
tree_sequence = msp.simulate(sample_size=5, recombination_rate=5, mutation_rate=5)

list(tree_sequence.haplotype())[i:j]  # haps similar to ms
tree_sequence.get_num_mutations()  # seg sites
tree_sequence.mutations()  # [0] is the location of the mutation
# have definition of pops as indexes 0:3, this takes an iterator
pop1 = tree_sequence.get_samples(1)  # pop1 is int, 1
tree_sequence.get_pairwise_diversity(samples=pop1)
tree_sequence.variants()  # position and string of column haps
with open("output.vcf", "w") as vcf_file:
    tree_sequence.write_vcf(vcf_file, 2)