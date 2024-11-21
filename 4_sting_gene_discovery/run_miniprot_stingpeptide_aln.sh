#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6000M
#SBATCH --cpus-per-task=16
#SBATCH --job-name=run_miniprot_stingpeptide_aln
#SBATCH --output=
#SBATCH --error=

#----------------------------------------------
###### using miniprot and tblastn to find stinging peotides in nettle genome #########
export PATH=$PATH:/home/~/bin/miniprot
#----------------------------------------------

miniprot -t16 -d Nettle_female_H2_Round_5_chrname_reordered_genome.mpi Nettle_female_H2_Round_5_chrname_reordered_genome.fa
#miniprot -Iut16 --gff Nettle_female_H2_Round_5_chrname_reordered_genome.mpi Urticaceae_uniprot_reviewed_proteins.fasta > Nettle_female_H2_Urticaceae_uniprot_reviewed_proteins_aln.gff

miniprot -Iut16 --aln Nettle_female_H1_Round_5_chrname_reordered_genome.mpi Urtica_toxins_peptides.fasta > Nettle_female_H1_Urtica_toxins_peptides.aln
miniprot -Iut16 --aln Nettle_female_H2_Round_5_chrname_reordered_genome.mpi Urtica_toxins_peptides.fasta > Nettle_female_H2_Urtica_toxins_peptides.aln
