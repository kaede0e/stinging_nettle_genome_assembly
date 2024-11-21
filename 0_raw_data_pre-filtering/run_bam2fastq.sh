#!/bin/bash
#SBATCH --account=$accountname
#SBATCH --time=1-0
#SBATCH --mem=30Gb
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/bam2fastq_nettle_male.12Mar2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/bam2fastq_nettle_male.12Mar2024.err

SMRT_ROOT=/home/~/bin/pacbio/smrtlink
export PATH=$PATH:/home/~/bin/pacbio/smrtlink/smrtcmds/bin/

bam2fastq -o Pacbio_hifi/Nettle_male_Pacbio_hifi Pacbio_hifi/Urtica_male_m84185_240301_023219_s3.hifi_reads.bam
