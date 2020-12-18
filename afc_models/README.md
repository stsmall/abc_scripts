# Models used in our PNAS paper on An. funestus complex (AFC)
[Radiation with reticulation marks the origin of a major malaria vector.](LINK)

Our paper on the An funestus complex deals with trying to place reticulations on a species tree when many of the reticulations are recursive (kind of). We approached this problem by using population genomic data and model selection. Here we present the three steps of our simulations.

[First](abcrf), we sought to find the species tree with reticulations that best explained the evolutionary history of the AFC. We focuses on three candidate models based on results of admxituregraph. These models are named according to the tree topology that they represent: topo1, topo3, topo7 (see [Fig2](LINK) for more details). Each of these models is run the same way using the code in the repository.  

`abc_sims.py -cfg AFC.abc.cfg -i 100 --ms msmove -m topo.model7.txt --out m7.msSims`

This will create a file with '100' lines (iterations) that can be run from the command line where 'model.topo7.txt' is used to set up the simulation demography, split times, and admxiture rates. It will use the coalescent simulator [msmove](https://github.com/geneva/msmove) and the input config file 'AFC.abc.cfg'. The output file will be named 'm7.smSims'. Statistics for the output file can then be calculated using 'abc_stats.py'.  

`abc_stats.py sim --infile ms7.msSims --outfile m7.stats.txt --pops 4 --pairs 0-1 0-2 0-3 --stats sfs jsfs --mask FOO.mask --gff FOO.gff --mode split-run`

If you need to calculate statistics from data to serve as the target for the approximate bayesian computation, you can use.  

`abc_stats.py obs YOUR.VCF CHR CHROM_LENGTH OUTPUT sample_pops.txt --pairs Moz Lik Van Lon Par --stats sfs jsfs --mask_file AFC.mask.fasta`


[Second](abc_tree7), we determined that the best-fit topology with reticulations was topology/tree vii (7). To generate simulations for approximate bayesian computation (ABC) to estimate parameters such as introgression times and divergence times simultaneously you can once again you the above procedure for generating simulations and then calculating summary statistics.  

`abc_sims.py -cfg AFC.abc.cfg -i 100 --ms msmove -m topo.model7.txt --out m7.msSims`

`abc_stats.py sim --infile ms7.msSims --outfile m7.stats.txt --pops 4 --pairs 0-1 0-2 0-3 --stats sfs jsfs --mask FOO.mask --gff FOO.gff --mode split-run`


[Third](filet), we used the ML classifier [FILET](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007341) to estimate which regions of the genome were introgressed and in which direction. In the filet directory you can find models for generating training data for each pairwise introgression event. Each file contains the three testing directions.

**mig12.msOut** migration from species 1 into species 2 backwards in time, representing introgression from species 2 into species 1 forward in time.

**mig21.msOut** migration from species 2 into species 1 backwards in time, representing introgression from species 1 into species 2 forward in time.

**noMig.msOut** no migration between the species pairs.

As above simulations are run using the `python abc_sims.py`. However, FILET has only been tested on the feature vector of summary statistics listed in the official FILET publication. Thus, I would recommend following the turorials and using the provided tools listed at www.github.com/kr-colab/FILET once you have the simulations.

*helper scripts* I have also written some helper scripts to take you from the simulations to the final training set [introgression](https://github.com/stsmall/An_funestus/tree/master/introgression)


## Other content in this directory

### AFC.pi.dist.txt
Distribution of nucleotide diversity in 100kb windows along the genome using [pi_xy](https://github.com/ksamuk/pixy)
### AFC.pi.noncoding.dist.txt
Distribution of nucleotide diversity in 100kb windows along the genome using pi_xy, retaining only windows that are 5kb from coding loci
### AFC.rho.dist.txt
Distribution of the population recombination rate in 100kb windows along the genome using [LDJump](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12994)
### AFC.rho.noncoding.dist.txt
Distribution of the population recombination rate in 100kb windows along the genome using LDJump retaining only windows that are 5kb from coding loci
### inversion_coordinates.Afunpv21.bed
A bed file listing the coordinates of the main polymorphic inversion in An. funestus: 3Ra, 3Rb, 2Ra, and 3La.
### migration_mat.txt
An example migration matrix file for using with the abc_sims.py
### AFC.abc.cfg
The generic config file used with abc_sims.py that lists the sample sizes and effective population sizes used for the PNAS publication.
