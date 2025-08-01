#!/bin/bash

#Get an interactive job if you have relatively large genomes (>100Mbp)
salloc --account=rrg-rieseber-ac --time=3:00:00 --ntasks=1 --cpus-per-task=4 --mem=20G

#in the interactive job, run the following: 
module load StdEnv/2023 bbmap/39.06 #we'll use the software called BBMap https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/
your_genome=Sclerotinia_minor_hap

stats.sh -Xmx20g $your_genome.fa format=2

#results will be printed on the screen, or you can save it as a text file. 
stats.sh -Xmx20g $your_genome.fa format=2 2> bbmap_output_${your_genome}.txt
