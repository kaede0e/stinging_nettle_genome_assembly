#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=30G
#SBATCH --output=
#SBATCH --error=

######### Genome repeat detection and soft-masking with RED ###########

#First, set up the working directory.
module load StdEnv/2020 gcc python/3.8 augustus hmmer blast+ metaeuk prodigal r
#virtualenv red_env
source ~/bin/red_env/bin/activate
#pip install biopython
#pip install natsort

export PATH=$PATH:/home/~/bin/RED_RepeatDetector/redUnix64
/home/~/bin/redmask/redmask.py \
-i Nettle_female_H2_Round_5_chrname_reordered_genome.fa \
-o Nettle_female_H2_Round_5_chrname_reordered_genome_repeatmsaked.fa
