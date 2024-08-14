#!/bin/bash
#SBATCH --time=1-0
#SBATCH --mem=20G
#SBATCH --account=def-gowens

module load python/3.10
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/python_env/bin/activate

fastq-filter -e 0.01 -o Nettle_female_Pacbio_hifi_Q20_filtered.fastq Nettle_female_Pacbio_hifi.fastq
