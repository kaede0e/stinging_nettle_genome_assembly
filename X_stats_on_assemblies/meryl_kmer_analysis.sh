#!/bin/bash

#### Analyzing polished assembly reuslts with Merqury ####

export PATH=$PATH:/home/khirabayashi/bin/merqury/eval
export PATH=$PATH:/home/khirabayashi/bin/merqury


#bash best_k.sh 574795222 kmer = 19.5314 so maybe go with 19 
meryl k=19 count output Nettle_f_Illumina_reads.meryl Nettle_f_Illumina_reads_fwd.1.fq
#calculate QV score for the assembly based on Illumina reads
qv.sh Nettle_f_Illumina_reads.meryl Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap1.p_ctg.fa
qv.sh Nettle_f_Illumina_reads.meryl Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap2.p_ctg.fa

ln -s /home/khirabayashi/bin/merqury/merqury.sh
meryl k=19 count /path/to/illumina/reads.fq output Nettle_female_Illumina_fwd.meryl
mkdir Nettle_female_Round_2_asm_trimmomatic_meryl_output
merqury.sh Nettle_female_Illumina_filtered_fwd.meryl \
Round_2_hap1.reviewed.chr_assembled.fasta \
Round_2_hap2.reviewed.chr_assembled.fasta \
Nettle_female_Round_2_asm_trimmomatic_meryl_output
