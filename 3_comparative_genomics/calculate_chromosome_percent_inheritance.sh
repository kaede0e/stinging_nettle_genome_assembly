#!/bin/bash

# This script is used to calculate the proportion of chromosome that got inherited to the offspring by using nucleotide alignment.

# Step 1: run minimap2 to align and produce .paf file (follow "run_minimap2_plot_aln_four_haplotypes.sh)
#eg) # Align female H1 to refgenome = $genome1 (FH1; only use 13 chromosomes, no scaffolds)
genome1=Nettle_female_H1_Round_8_genome.chr.fa
genome2=Nettle_female_H1_Round_8_genome_renamed.fa

#minimap2 -cx asm5 $genome1 $genome2 -t 40 > Nettle_female_H1_chr_vs_female_Round8_H1.paf

# Step 2: calculate the sum of alingnment length, only if alignment quality is >30 and % identity is >99.7% (these cutoffs are adjustable)
cat Nettle_female_H2_chr_vs_male_Round3_H1.paf \
| awk '{if ($1 == "M_H1_chr_01" && $6 == "Urtica_dioica_female_chr_01"){print}}' \ #makes sure you're only comparing homologous pair
| awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5 }' \ #calculates % match based on (number of matching nucleitides/alignment lengthx100%)
| sed s/' '/'  '/g | awk '{if ($6 >= 99.7){print}}' \ #if % nucleotide match for teh alignment is >99.7%
| cut -f 5 | paste -sd+ | bc #sum of the alignment length

#let me do this in a loop for all chromosomes... 
for chr in {1..9};
do
  cat Nettle_female_H1_chr_vs_male_Round3_H1.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H1_chr_0"chr"" && $6 == "Urtica_dioica_female_chr_0"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH1vsMH1", "chr"chr"", $1)}' > nucl_match/FH1_vs_MH1_${chr}_matching_nucl.txt;
  
  cat Nettle_female_H2_chr_vs_male_Round3_H1.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H1_chr_0"chr"" && $6 == "Urtica_dioica_female_chr_0"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH2vsMH1", "chr"chr"", $1)}' > nucl_match/FH2_vs_MH1_${chr}_matching_nucl.txt;
  
  cat Nettle_female_H1_chr_vs_male_Round3_H2.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H2_chr_0"chr"" && $6 == "Urtica_dioica_female_chr_0"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH1vsMH2", "chr"chr"", $1)}' > nucl_match/FH1_vs_MH2_${chr}_matching_nucl.txt;
  
  cat Nettle_female_H2_chr_vs_male_Round3_H2.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H2_chr_0"chr"" && $6 == "Urtica_dioica_female_chr_0"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH2vsMH2", "chr"chr"", $1)}' > nucl_match/FH2_vs_MH2_${chr}_matching_nucl.txt;
done
for chr in {10..13};
do
  cat Nettle_female_H1_chr_vs_male_Round3_H1.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H1_chr_"chr"" && $6 == "Urtica_dioica_female_chr_"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH1vsMH1", "chr"chr"", $1)}' > nucl_match/FH1_vs_MH1_${chr}_matching_nucl.txt;
  
  cat Nettle_female_H2_chr_vs_male_Round3_H1.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H1_chr_"chr"" && $6 == "Urtica_dioica_female_chr_"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH2vsMH1", "chr"chr"", $1)}' > nucl_match/FH2_vs_MH1_${chr}_matching_nucl.txt;
  
  cat Nettle_female_H1_chr_vs_male_Round3_H2.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H2_chr_"chr"" && $6 == "Urtica_dioica_female_chr_"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH1vsMH2", "chr"chr"", $1)}' > nucl_match/FH1_vs_MH2_${chr}_matching_nucl.txt;
  
  cat Nettle_female_H2_chr_vs_male_Round3_H2.paf \
  | awk -v chr="$chr" '{if ($1 == "M_H2_chr_"chr"" && $6 == "Urtica_dioica_female_chr_"chr""){print}}' \
  | awk '$12 >= 30 {print $8, $9, $3, $4, $4-$3, ($10*100)/$11, $1, $6, $5}' | sed s/' '/'       '/g \
  | awk '{if ($6 >= 99.7){print}}' | cut -f 5 | paste -sd+ | bc \
  | awk -v chr="$chr" '{print ("FH2vsMH2", "chr"chr"", $1)}' > nucl_match/FH2_vs_MH2_${chr}_matching_nucl.txt;
done

cat nucl_match/*_matching_nucl.txt | sort > all_chr_nucl_match.txt
rm -r nucl_match/




