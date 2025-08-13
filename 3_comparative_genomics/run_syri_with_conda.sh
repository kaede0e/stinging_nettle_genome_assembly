#!/bin/bash

####### SyRI for visualization of SVs from minimap2 reasults #######
#This was downloaded in salal using conda.

#conda create -n syri_env -c bioconda syri #installed to /home/khirabayashi/miniforge3/envs
#conda deactivate # to exit from current conda env.
#conda activate syri_env

refgenome=Nettle_male_H1_Round_3_genome_reoriented_chrrenamed_for_syri.chr.fa
qrygenome=Nettle_male_H2_Round_3_genome_reoriented_chrrenamed_for_syri.chr.fa
prefix=maleH1_vs_maleH2

cwd="."     # Change to working directory
cd $cwd
PATH_TO_SYRI="/home/khirabayashi/miniforge3/envs/syri_env/bin/syri" #Change the path to point to syri

ln -sf Nettle_male_H1_Round_3_genome_reoriented_chrrenamed_for_syri.chr.fa maleH1 #for labelling what you're alining
ln -sf Nettle_male_H2_Round_3_genome_reoriented_chrrenamed_for_syri.chr.fa maleH2
#ln -sf Salmonberry_redmorph_H1.v1.chrrenamed_for_syri.chr.fa RedH1
#ln -sf Salmonberry_redmorph_H2.v1.chrrenamed_for_syri.chr.fa RedH2

minimap2 -ax asm5 --eqx $refgenome $qrygenome -t 60 > ${prefix}_minimap2_aln.sam
#minimap2 -ax asm5 --eqx $refgenome $qrygenome -t 60 > ${refgenome}_${qrygenome}_minimap2_aln.sam #if you don't like using prefix variable; but gives you very long names on the output files

## Run SyRI round 1 with SAM or BAM file as input
python3 $PATH_TO_SYRI -c ${prefix}_minimap2_aln.sam -r $refgenome -q $qrygenome --prefix ${prefix} -k -F S 2> syri_error1_${prefix}_out.txt
#python3 $PATH_TO_SYRI -c ${refgenome}_${qrygenome}_minimap2_aln.sam -r $refgenome -q $qrygenome --prefix ${refgenome}_vs_${qrygenome} -k -F S 2> syri_error1_${refgenome}_${qrygenome}_out.txt
