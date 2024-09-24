#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=192000M
#SBATCH --account=def-gowens
#SBATCH --cpus-per-task=48
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_male_kmer_jellyfish.02May2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_male_kmer_jellyfish.02May2024.err

### kmer plot with jellyfish ###

#####################################
##### Variables / data ##############
#####################################

fastq_gz=/home/kaedeh/scratch/Nettle/Pacbio_hifi/Nettle_male_Pacbio_hifi.fastq.gz
zcat $fastq_gz > Nettle_male_Pacbio_hifi.fastq
export PATH=PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/jellyfish-2.3.1/jellyfish

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

module --force purge
module load StdEnv/2020 jellyfish/2.3.0

jellyfish count -C -m 21 -s 18000000 -t 46 *.fastq -o Nettle_male_hifi_reads.jf

jellyfish histo -t 40 Nettle_male_hifi_reads.jf > Nettle_male_hifi_reads.histo
