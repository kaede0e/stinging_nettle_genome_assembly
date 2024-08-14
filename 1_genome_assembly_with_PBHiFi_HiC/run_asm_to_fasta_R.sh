#!/bin/bash
#SBATCH --account=def-gowens
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/asm_to_fasta_nettle_female.06Aug2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/asm_to_fasta_nettle_female.06Aug2024.err

module load r/4.3

##### Convert .assembly to .fasta with R script by Moshutava and Eric ######

#mkdir FASTA
Rscript asm_to_fasta_me.R Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap1.p_ctg.0.review.4.assembly Round_3 ../references/*.fa FASTA
