#!/bin/bash
#SBATCH --account=
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --output=
#SBATCH --error=

module load r/4.3

##### Convert .assembly to .fasta with R script by Moshutava and Eric ######

#mkdir FASTA
Rscript asm_to_fasta_me.R Nettle_female_hifi_hifiasm_0.19.8_Q20reads.asm.hic.hap1.p_ctg.0.review.4.assembly Round_3 ../references/*.fa FASTA


#-------WARNING: The asm_fasta_me.R alone doesn't produce an accurate fasta when the each separate contig was assigned to a separate "superscaffold". 
#ie. when there are small contigs in Juicebox that are grouped into one artificial scaffolds, the R script will  merge them together in a single line in fasta starting with >, nmessing up my contig/scafofld counts...
#solution: 
#every time a contig is merged in a scaffold, 500 N's are added as an error range by their R script. 
#so what I am trying to do is to replace the N's with a new line in fasta. 
#--------

module load bbmap

cat random_section_with_Ns.fa | \
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
sed s/'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'/'\n>newline\n'/g | \
sed s/'  '/'\n'/g | fold -w 80 > random_section_replaced_Ns.fa 

rename.sh in=random_section_replaced_Ns.fa out=random_section_replaced_Ns_renamed.fa prefix=Urtica_dioica_female_debris

