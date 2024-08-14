#!/bin/bash
#SBATCH --account=def-mtodesco
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i NettleFemale -f Round_2_hap1.reviewed.chr_assembled.chr.fasta -h H1 -c 15 -m 128000 -g FALSE &> RepeatOBserver_Hap1_log.txt
