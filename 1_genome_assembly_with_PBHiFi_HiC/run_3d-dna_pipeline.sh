#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH --account=def-gowens
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_3ddna_hifiasm_4_hap1.13May2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_3ddna_hifiasm_4_hap1.13May2024.err

### 3D-DNA pipeline after getting Hi-C contact merge_nodups.txt file ###

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

#####################################
### Execution of programs ###########
#####################################

module load StdEnv/2020 python/3.11.2 java/17.0.2 lastz/1.04.03
#virtualenv 3ddna_env
#pip install scipy numpy matplotlib #libraries required for 3d-dna
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/3ddna_env/bin/activate
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/3d-dna #3D de novo assembly: version 180114

#####################################
##### Variables / data ##############
#####################################
contig_fasta=/home/kaedeh/scratch/Nettle/HiC_hap1/references/*.fa
merged_nodups=/home/kaedeh/scratch/Nettle/HiC_hap1/aligned/merged_nodups.txt

#run script; -r 0 runs only the first scaffolding, no polishing.
run-asm-pipeline.sh -m haploid -r 0 $contig_fasta $merged_nodups

deactivate

# ---------------------------------------------------------------------
echo "Done 3D-DNA pipeline scaffolding.  Use Juicebox to visualize and manually fix scaffolds."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
