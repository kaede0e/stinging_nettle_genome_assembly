#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --account=
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=GENESPACE_comparative_genomics_cont
#SBATCH --output=
#SBATCH --error=


# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

##### Make sure required modules are installed for Orthofiner, MCScan, to run GENESPACE #####
module load StdEnv/2020 gcc/9.3.0 python/3.8 java/13.0.2 diamond/2.0.15 r/4.3.1
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/smcpp_python_env/bin/activate
export PATH=$PATH:/home/~/bin/OrthoFinder_source
export PATH=$PATH:/home/~/bin/MCScanX
export PATH=$PATH:/home/~/bin/MCScanX/downstream_analyses

Rscript genespace_with_orthofinder_pre-parsed.R
