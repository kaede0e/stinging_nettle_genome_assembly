#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --mem=125G
#SBATCH --nodes=1
#SBATCH --account=def-gowens
#SBATCH --ntasks-per-node=32
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_juicer_hifiasm_4_hap1.11May2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_juicer_hifiasm_4_hap1.11May2024.err

## Juicer and 3D-DNA pipelne for genome scaffolding ##

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

#####################################
### Execution of programs ###########
#####################################

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1

#run juicer
bash scripts/juicer.sh -D $PWD \
-g Nettle_f_hifiasm_4_hap1 -s DpnII \
-z references/Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap1.p_ctg.fa \
-y restriction_sites/Nettle_f_hifiasm_4_hap1_DpnII.txt \
-p restriction_sites/Nettle_f_hifiasm_4_hap1_DpnII.chrom.sizes \
-t 28 -S early
