#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000M

# BUSCO with contig file (fasta)
# BUSCO requires you to work in a Python Wheel environment.
# So first, you create a Python environment in your workspace.
# you're suggested to work in a new screen so that all softwares are properly loaded.
module --force purge
module load StdEnv/2020 gcc python augustus hmmer blast+ metaeuk prodigal r
source /home/~/scripts/busco_env/bin/activate

# Run BUSCO job:
busco --offline \
--in /home/kaedeh/scratch/Nettle/assembly/hifiasm_asm_4/Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap1.p_ctg.fa \
--out Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap1._BUSCO_output \
--lineage_dataset /home/~/busco_downloads/lineages/eudicots_odb10 \
--mode genome \
--cpu ${SLURM_CPUS_PER_TASK-1}


# Once job is done, this will generate a directory *BUSCO_output
# in the directory you're calling this command in.

# ---------------------------------------------------------------------
echo "Finished BUSCO analysis: `date`"
# ---------------------------------------------------------------------
echo ""

# ---------------------------------------------------------------------
echo "Check scores in *program_assembly_BUSCO_output/short_summary"
# ---------------------------------------------------------------------
echo ""
