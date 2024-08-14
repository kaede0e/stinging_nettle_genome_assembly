#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=192000M
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=48
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hifiasm_asm_8_with_HiC_and_ONT.08Jul2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hifiasm_asm_8_with_HiC_and_ONT.08Jul2024.err

### HiFiasm + Hi-C genome assembly pipeline ###

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

module load StdEnv/2023
export PATH=PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/Hifiasm/hifiasm-0.19.8

#####################################
##### Variables / data ##############
#####################################

HiC_1=/home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/HiC/HiC_Nettle_FemQC1_FKDL240020059-1A_223GFKLT4_L1_1.fq.gz
HiC_2=/home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/HiC/HiC_Nettle_FemQC1_FKDL240020059-1A_223GFKLT4_L1_2.fq.gz
Hifi_fastq=/home/kaedeh/scratch/Nettle/Pacbio_hifi/Nettle_female_Pacbio_hifi_Q30_filtered.fastq
ONT_30kbp=/home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/ONT_nettle_baecalled_reads_passed/nettle_basecalled_reads_combined_30kb_or_longer.fastq

hifiasm -o Nettle_female_hifi_hifiasm_0.19.8_ONT30kb_HiC.asm -t 48 --h1 $HiC_1 --h2 $HiC_2  $Hifi_fastq

# ---------------------------------------------------------------------
echo "Done assembly with Hifiasm. Use manual programs to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
