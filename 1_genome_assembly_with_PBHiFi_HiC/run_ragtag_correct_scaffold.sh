#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --ntasks=1
#SBATCH --account=rrg-rieseber-ac
#SBATCH --cpus-per-task=4
#SBATCH --job-name=Nettle_assembly_ragtag_with_DToL_refgenome
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female.ragtag.correct.scaffold.hap1.20Sep2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female.ragtag.correct.scaffold.hap1.20Sep2024.err

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

##### Ragtag - reference genome based scaffolding #####

module load StdEnv/2020 minimap2/2.24
source /home/kaedeh/projects/def-gowens/kaedeh/Nettle/scripts/busco_env/bin/activate

refgenome=GCA_964188135.1_drUrtDioi1.1_genomic.chrname_modified.chr.fa
qrygenome_H1=~/projects/def-gowens/kaedeh/Nettle/output/assembly/Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap1.p_ctg.fa
qrygenome_H2=~/projects/def-gowens/kaedeh/Nettle/output/assembly/Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap2.p_ctg.fa

ragtag.py correct \
$refgenome \
$qrygenome_H1 \
-u -o ragtag/Nettle_female_hifi_hifiasm_hap1_DToL
echo "Ragtag correct finished"

ragtag.py scaffold \
$refgenome \
ragtag/Nettle_female_hifi_hifiasm_hap1_DToL/ragtag.correct.fasta \
-u -o ragtag/Nettle_female_hifi_hifiasm_hap1_DToL
echo "Ragtag scaffold finished"
