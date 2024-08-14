#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=rrg-rieseber-ac
#SBATCH --ntasks=9
#SBATCH --mem-per-cpu=50G
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_to_hap2_anchorwave_with_hap1genes.01Aug2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_to_hap2_anchorwave_with_hap1genes.01Aug2024.err


#----------- Anchorwave whole genome alignment from chromosome-level assembly data -----------#

# Anchorwave - whole genome alignment from chromosome-level assembly data
# This software allows you to align two genomes using coding sequence (CDS),
# and the output can be tabulated and exported, visualized on R.

module load StdEnv/2023 minimap2/2.28
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/anchorwave/
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/anchorwave/scripts

#0 Set variables
Genome=Nettle_female
Hap1=/home/kaedeh/scratch/Nettle/annotation/final_annotation_files/Round_2_hap1.reviewed.chr_assembled.fasta
Hap1_gene_anno_gff3=/home/kaedeh/scratch/Nettle/annotation/final_annotation_files/Nettle_female_hap1_genes_anno.gff3
Hap2=/home/kaedeh/scratch/Nettle/HiC_hap2/3d_dna_pipeline/FASTA/Round_2_hap2.reviewed.chr_assembled.fasta

#1 extracting CDS on both haplotypes
#anchorwave gff2seq \
#-i $Hap1_gene_anno_gff3 -r $Hap1 -o anchorwave_gff2seq_${Genome}_H1.cds.fasta
#anchorwave gff2seq \
#-i $Hap1_gene_anno_gff3 -r $Hap2 -o anchorwave_gff2seq_${Genome}_H2.cds.fasta

#2 aligning CDS on genome
#minimap2 -x splice -t 4 -k 12 -a -p 0.4 -N 20 \
#$Hap1 anchorwave_gff2seq_${Genome}_H1.cds.fasta > anchorwave_minimap2_${Genome}_H1.cds.sam
#minimap2 -x splice -t 4 -k 12 -a -p 0.4 -N 20 \
#$Hap2 anchorwave_gff2seq_${Genome}_H1.cds.fasta > anchorwave_minimap2_${Genome}_H2.cds.sam

#3 visualize the alignment to assess synteny on R by transforming SAM file to Tabulated file
#perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/anchorwave/scripts/alignmentToDotplot.pl \
#$Hap1_gene_anno_gff3 anchorwave_minimap2_${Genome}_H1.cds.sam > anchorwave_minimap2_${Genome}_H1.cds.tab
#perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/anchorwave/scripts/alignmentToDotplot.pl \
#$Hap1_gene_anno_gff3 anchorwave_minimap2_${Genome}_H2.cds.sam > anchorwave_minimap2_${Genome}_H2.cds.tab

#4 genome alignment with relocation variation, chromosome fusion or whole genome duplication (proali)
anchorwave proali \
-i $Hap1_gene_anno_gff3 \
-as anchorwave_gff2seq_${Genome}_H1.cds.fasta \
-r $Hap1 \
-a anchorwave_minimap2_${Genome}_H2.cds.sam \
-ar anchorwave_minimap2_${Genome}_H1.cds.sam \
-s $Hap2 \
-n ${Genome}_Hap1_vs_Hap2.anchors \
-R 1 -Q 1 -t 6 \
-o anchorwave_proali_${Genome}_Hap1_vs_Hap2.maf \
-f anchorwave_proali_${Genome}_Hap1_vs_Hap2.f.maf > anchorwave_proali_log.txt
