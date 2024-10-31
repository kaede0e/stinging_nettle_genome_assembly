#!/bin/bash

####### SyRI for visualization of SVs from minimap2 reasults #######

#This was downloaded in salal using conda. 

#conda create -n syri_env -c bioconda syri #installed to /home/khirabayashi/miniforge3/envs
#conda deactivate # to exit from current conda env. 
conda activate syri_env

cwd="."     # Change to working directory
cd $cwd
PATH_TO_SYRI="/home/khirabayashi/miniforge3/envs/syri_env/bin/syri" #Change the path to point to syri

## NOTE: IDEALLY, SYRI EXPECTS THAT THE HOMOLOGOUS CHROMOSOMES IN THE TWO GENOMES
## WOULD HAVE EXACTLY SAME CHROMOSOME ID. THEREFORE, IT IS HIGHLY RECOMMENDED THAT
## THE USER PRE-PROCESSES THE FASTA FILES TO ENSURE THAT HOMOLOGOUS CHROMOSOMES
## HAVE EXACTLY THE SAME ID IN BOTH FASTA FILES CORRESPONDING TO THE TWO GENOMES.
## IN CASE, THAT IS NOT THE CASE, SYRI WOULD TRY TO FIND HOMOLOGOUS GENOMES USING
## WHOLE GENOME ALIGNMENTS, BUT THAT METHOD IS HEURISTICAL AND CAN RESULT IN
## SUBOPTIMAL RESULTS.
## ALSO, IT IS REQUIRED THAT THE TWO GENOMES (FASTA FILES) HAVE THE SAME NUMBER
## OF CHROMOSOMES.

ln -sf GCA_964188135.1_drUrtDioi1.1_genomic.fna refgenome
ln -sf Nettle_female_H1_Round_5_chrname_reordered_genome.chr.fa qrygenome_H1
ln -sf Nettle_female_H2_Round_5_chrname_reordered_genome.chr.fa qrygenome_H2

## Perform whole genome alignment
# Using minimap2 for generating alignment. Any other whole genome alignment tool can also be used.
minimap2 -ax asm5 --eqx refgenome qrygenome_H1 > Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_aln.sam
minimap2 -ax asm5 --eqx refgenome qrygenome_H2 > Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_aln.sam

## Note: It is recommended that the user tests different alignment settings to find what
## alignment resolution suits their biological problem. Some alignment tools find
## longer alignments (with lots of gaps) while other find smaller more fragmented
## alignments. The smaller alignments generally have higher alignment identity scores
## and are more helpful in identifying smaller genomic structural rearrangments.
## But they could also lead to significant increase in redundant alignments which
## leads to increase in runtime of the alignment tool and SyRI.

## Run SyRI round 1 with SAM or BAM file as input
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_aln.sam -r refgenome -q qrygenome_H1 --prefix DToLpri_Round5_H1 -k -F S 2> syri_error1_H1_out.txt
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_aln.sam -r refgenome -q qrygenome_H2 --prefix DToLpri_Round5_H2 -k -F S 2> syri_error1_H2_out.txt

# OR
#samtools view -b out.sam > out.bam
#python3 $PATH_TO_SYRI -c out.bam -r refgenome -q qrygenome -k -F B

## Note: SyRI requires your chromosomes to be aligned in the syntenic direction. If one is reversed, you need to reverse-complement it to proceed.
cat syri_error1_H1_out.txt | grep "ERROR" | cut -f 11 -d " " | sed s/"\."/""/g > chr_to_rev_H1.txt ;
mkdir qrygenome_H1_chr
cat Nettle_female_H2_Round_5_syri_input_genome_chrnames.txt | grep -v -f chr_to_rev_H1.txt > chr_to_keep_H1.txt
  for chr in `cat chr_to_keep_H1.txt` ;
  do
    echo $chr > qrygenome_H1_chr/$chr.txt ;
    samtools faidx qrygenome_H1 -r qrygenome_H1_chr/$chr.txt > qrygenome_H1_chr/$chr.fwd.fa ;
  done #forward strand fasta for each chromosome
  for chr in `cat chr_to_rev_H1.txt` ;
  do
    echo $chr > qrygenome_H1_chr/rev_$chr.txt ;
    samtools faidx qrygenome_H1 -i -r qrygenome_H1_chr/rev_$chr.txt > qrygenome_H1_chr/$chr.rev.fa ;
  done
  cat qrygenome_H1_chr/*.fa > qrygenome_H1_chr_rev.fa
cat qrygenome_H1_chr_rev.fa | sed s/'\/rc'/''/g > qrygenome_H1_chr_rev_renamed.fa

cat syri_error1_H2_out.txt | grep "ERROR" | cut -f 11 -d " " | sed s/"\."/""/g > chr_to_rev_H2.txt ;
mkdir qrygenome_H2_chr
cat Nettle_female_H2_Round_5_syri_input_genome_chrnames.txt | grep -v -f chr_to_rev_H2.txt > chr_to_keep_H2.txt
  for chr in `cat chr_to_keep_H2.txt` ;
  do
    echo $chr > qrygenome_H2_chr/$chr.txt ;
    samtools faidx qrygenome_H2 -r qrygenome_H2_chr/$chr.txt > qrygenome_H2_chr/$chr.fwd.fa ;
  done #forward strand fasta for each chromosome
  for chr in `cat chr_to_rev_H2.txt` ;
  do
    echo $chr > qrygenome_H2_chr/rev_$chr.txt ;
    samtools faidx qrygenome_H2 -i -r qrygenome_H2_chr/rev_$chr.txt > qrygenome_H2_chr/$chr.rev.fa ;
  done
  cat qrygenome_H2_chr/*.fa > qrygenome_H2_chr_rev.fa
cat qrygenome_H2_chr_rev.fa | sed s/'\/rc'/''/g > qrygenome_H2_chr_rev_renamed.fa

###This qrygenome fasta has weirdly retained chromosome names from empty liens also and messed up the fasta.
###Manually inspect the fasta so that headers are all correctly in fasta format.

#Re-do whole genome alignment with the reverse complemented fasta
minimap2 -ax asm5 --eqx refgenome qrygenome_H1_chr_rev_renamed.fa -t 10 > Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_rev2_aln.sam
minimap2 -ax asm5 --eqx refgenome qrygenome_H2_chr_rev_renamed.fa -t 10 > Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_rev2_aln.sam

## Run SyRI round 2 with the new SAM file as input NOTE: SYRI doesn't take sam files with "." in in front of .sam
#python3 $PATH_TO_SYRI -c Nettle_female_Round_5_H1_H2_minimap2_rev2_aln.sam -r refgenome -q qrygenome_chr_rev_renamed.fa -k -F S 2> syri_error2_out.txt
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_rev2_aln.sam -r refgenome -q qrygenome_H1_chr_rev_renamed.fa --prefix DToLpri_Round5_H1_rev -k -F S 2> syri_error2_H1_out.txt
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_rev2_aln.sam -r refgenome -q qrygenome_H2_chr_rev_renamed.fa --prefix DToLpri_Round5_H2_rev -k -F S 2> syri_error2_H2_out.txt

## SyRI would report genomic structural differences in syri.out and syri.vcf.

## Plot the syri.out using plotsr - also installed in the SyRI conda env.
#plotsr --sr syri.out --genomes genomes.txt -o syri_plotsr_output_plot.png -S 0.5 -W 7 -H 10 -f 8


## Using SyRI to identify genomic rearrangements from whole-genome alignments
## generated using MUMmer. A .tsv (out.filtered.coords) file is used as the input.
#nucmer --maxmatch -c 100 -b 500 -l 50 refgenome qrygenome       # Whole genome alignment. Any other alignment can also be used.
#delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta     # Remove small and lower quality alignments
#show-coords -THrd out.filtered.delta > out.filtered.coords      # Convert alignment information to a .TSV format as required by SyRI
#python3 $PATH_TO_SYRI -c out.filtered.coords -d out.filtered.delta -r refgenome -q qrygenome
