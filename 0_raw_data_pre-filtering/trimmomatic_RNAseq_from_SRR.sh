#!/bin/bash

#cDNA available from Nettle Xu et al. 2021 transcriptoe study (https://www.mdpi.com/1422-0067/22/22/12319)
#SRR8511647- Leaf
java -Xmx100g -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE Urtica_dioica_leaves_SRR8511647_1.fastq Urtica_dioica_leaves_SRR8511647_2.fastq \
Urtica_dioica_leaves_SRR8511647_1_paired_trimmomatic.fastq.gz Urtica_dioica_leaves_SRR8511647_1_unpaired_trimmomatic.fastq.gz \
Urtica_dioica_leaves_SRR8511647_2_paired_trimmomatic.fastq.gz Urtica_dioica_leaves_SRR8511647_2_unpaired_trimmomatic.fastq.gz \
-trimlog Urtica_dioica_leaves_Illumina_trimmomatic_output.log -summary Urtica_dioica_leaves_Illumina_trimmomatic_output_stats.txt \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

#SRR8511647- Roots
java -Xmx100g -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE Urtica_dioica_roots_SRR8511650_1.fastq Urtica_dioica_roots_SRR8511650_2.fastq \
Urtica_dioica_roots_SRR8511650_1_paired_trimmomatic.fastq.gz Urtica_dioica_roots_SRR8511650_1_unpaired_trimmomatic.fastq.gz \
Urtica_dioica_roots_SRR8511650_2_paired_trimmomatic.fastq.gz Urtica_dioica_roots_SRR8511650_2_unpaired_trimmomatic.fastq.gz \
-trimlog Urtica_dioica_roots_Illumina_trimmomatic_output.log -summary Urtica_dioica_roots_Illumina_trimmomatic_output_stats.txt \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

#SRR8511647- Fibres
java -Xmx100g -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE Urtica_dioica_fibres_SRR8511695_1.fastq Urtica_dioica_fibres_SRR8511695_2.fastq \
Urtica_dioica_fibres_SRR8511695_1_paired_trimmomatic.fastq.gz Urtica_dioica_fibres_SRR8511695_1_unpaired_trimmomatic.fastq.gz \
Urtica_dioica_fibres_SRR8511695_2_paired_trimmomatic.fastq.gz Urtica_dioica_fibres_SRR8511695_2_unpaired_trimmomatic.fastq.gz \
-trimlog Urtica_dioica_fibres_Illumina_trimmomatic_output.log -summary Urtica_dioica_fibres_Illumina_trimmomatic_output_stats.txt \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
