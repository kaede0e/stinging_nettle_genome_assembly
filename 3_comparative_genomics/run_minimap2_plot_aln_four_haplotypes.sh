#!/bin/bash

#----------------------- This is used in the male nettle genome manuscript to plot the recombination landscape using alignment to female assembly. --------------# 

#Two genomes you're trying to align:
#1) FH1 (female haplotype 1)
#2) FH2 (female haplotype 2)

#Step 1: do pairwise alignment using minimap2
# Align female H1 to refgenome = $genome1 (FH1; only use 13 chromosomes, no scaffolds)
genome1=Nettle_female_H1_Round_8_genome.chr.fa
genome2=Nettle_female_H1_Round_8_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H1_chr_vs_female_Round8_H1.paf

# Align female H2
genome1=Nettle_female_H1_Round_8_genome.chr.fa
genome2=Nettle_female_H2_Round_8_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H1_chr_vs_female_Round8_H2.paf

# Align male H1
genome1=Nettle_female_H1_Round_8_genome.chr.fa
genome2=Nettle_male_H1_Round_3_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H1_chr_vs_male_Round3_H1.paf

# Align male H2
genome1=Nettle_female_H1_Round_8_genome.chr.fa
genome2=Nettle_male_H2_Round_3_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H1_chr_vs_male_Round3_H2.paf

# Repeat with female H2 to refgenome = $genome1 (FH2; only use 13 chromosomes, no scaffolds)
# Align female H1
genome1=Nettle_female_H2_Round_8_genome.chr.fa
genome2=Nettle_female_H1_Round_8_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H2_chr_vs_female_Round8_H1.paf

# Align female H2
genome1=Nettle_female_H2_Round_8_genome.chr.fa
genome2=Nettle_female_H2_Round_8_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H2_chr_vs_female_Round8_H2.paf

# Align male H1
genome1=Nettle_female_H2_Round_8_genome.chr.fa
genome2=Nettle_male_H1_Round_3_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H2_chr_vs_male_Round3_H1.paf

# Align male H2
genome1=Nettle_female_H2_Round_8_genome.chr.fa
genome2=Nettle_male_H2_Round_3_genome_renamed.fa

minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H2_chr_vs_male_Round3_H2.paf


#Step 2: convert paf to table to plot and select alignments with MQ > 30, and add a percent alignment column calculated by $10/$11*100 (PAF file matching nucleotides / alignment length)
for i in $(ls Nettle_female_H*_chr_vs_*.paf)
do
        name=$(echo $i | cut -d "." -f 1);
        awk '$12 >= 30 {print $8, $9, $3, $4, ($10*100)/$11, $1, $6, $5 }' $i > "$name""MQ30_with_percentID_tab.txt";
done


#Step 3: Format tables for plots
#combine the alignment files into one, so that you can make a single plot comparing four haplotypes as long as the reference genome used in alignment is the same for all haplotypes.
cat Nettle_female_H1_chr_vs_*MQ30_with_percentID_tab.txt > Nettle_female_H1_chr_vs_four_haplotypes_with_percentID_tab.txt
cat Nettle_female_H2_chr_vs_*MQ30_with_percentID_tab.txt > Nettle_female_H2_chr_vs_four_haplotypes_with_percentID_tab.txt

#for recombination plot, filter the alignment by >95% to show the difference.
for i in $(ls *_vs_four_haplotypes_with_percentID_tab.txt);
do
        name=$(echo $i | cut -d "." -f 1);
        cat $i | awk '{if ($5 >= 95){print}}' > "$name""95.txt";
done


#Step 4: make plots with Eric's Rscript for alignment coloured by percent ID
for j in $( ls *95.txt)
do
        Rscript PLOT_cont2_allCHR_PER.R $j
done
