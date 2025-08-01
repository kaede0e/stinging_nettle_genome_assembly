# BUSCO with contig file (fasta)
# BUSCO requires you to work in a Python Wheel environment.
# So first, you create a Python environment in your workspace.
# you're suggested to work in a new screen so that all softwares are properly loaded.
module load StdEnv/2020 gcc/9.3.0 python/3.10 augustus/3.5.0 hmmer/3.3.2 blast+/2.13.0 metaeuk/6 prodigal/2.6.3 r/4.3.1 bbmap/38.86
virtualenv busco_env_v2
source busco_env/bin/activate
pip install --no-index biopython==1.81 pandas==2.1.0 busco==5.5.0 #as of July 2025 (thanks Vedin for figuring out!)

# If working properly, you will be directed to a python environment like this:
## (busco_env) [kaedeh@cdr767 strawberry_genome]$

# Download the appropriate lineage dataset by: 
busco --download eudicots_odb10


# You're good to proceed.

#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000M

# Run BUSCO job: (fungal genome example)
busco --offline \
--in with_newhic/Sclerotinia_minor_hifiasm_1.asm.hic.hap1.p_ctg.fa \
--out Sclerotinia_minor_hifiasm1_newhic_H1.BUSCO_output \
--lineage_dataset /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/busco_env/reference_data/busco_downloads/lineages/helotiales_odb10 \
--mode genome --cpu 4

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
