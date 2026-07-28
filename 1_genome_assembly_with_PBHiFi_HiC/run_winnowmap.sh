#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=24
#SBATCH --account=
#SBATCH --output=
#SBATCH --error=

### quality check the Juicer+3ddna pipeline output by aligning two haplotypes ###

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

#####################################
### Execution of programs ###########
#####################################

module load samtools

export PATH=$PATH:/home/~/bin/Winnowmap/bin
HIFI_reads=/home/kaedeh/scratch/Nettle/Pacbio_hifi/Nettle_female_Pacbio_hifi_Q30_filtered.fastq
ONT_reads=/home/kaedeh/scratch/Nettle/HiC_hap1_hap2/references/Nettle_female_canu_correctedReads_10kb.fq
genome_H1=/home/kaedeh/scratch/Nettle/HiC_hap1_hap2/references/round_5/Nettle_female_H1_Round_5_syri_input_genome.fa
genome_H2=/home/kaedeh/scratch/Nettle/HiC_hap1_hap2/references/round_5/Nettle_female_H2_Round_5_syri_input_genome.fa

#1. build a merylDB from your ONT/HIFI reads that you want to map
meryl k=19 count output meryl_k19_nettle_female_H1.merylDB $genome_H1
meryl k=19 count output meryl_k19_nettle_female_H2.merylDB $genome_H2

#2. prepare repeats for winnowmap
meryl print greater-than distinct=0.9998 meryl_k19_nettle_female_H1.merylDB > repetitive_k19_H1.txt
meryl print greater-than distinct=0.9998 meryl_k19_nettle_female_H2.merylDB > repetitive_k19_H2.txt

#3. run winnowmap to map HiFi and ONT reads to the assembly
winnowmap -k 19 -W repetitive_k19_H1.txt -ax map-ont $genome_H1 $ONT_reads > winnowmap/Nettle_female_H1_v1_winnowmap_ONT.sam
samtools sort Nettle_female_H1_v1_winnowmap_ONT.sam --threads 22 > Nettle_female_H1_v1_winnowmap_ONT.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H1_v1_winnowmap_ONT.sorted.bam.bai Nettle_female_H1_v1_winnowmap_ONT.sorted.bam
winnowmap -k 19 -W repetitive_k19_H2.txt -ax map-ont $genome_H2 $ONT_reads > winnowmap/Nettle_female_H2_v1_winnowmap_ONT.sam
samtools sort Nettle_female_H2_v1_winnowmap_ONT.sam --threads 22 > Nettle_female_H2_v1_winnowmap_ONT.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H2_v1_winnowmap_ONT.sorted.bam.bai Nettle_female_H2_v1_winnowmap_ONT.sorted.bam
winnowmap -k 19 -W repetitive_k19_H1.txt -ax map-pb $genome_H1 $HIFI_reads > winnowmap/Nettle_female_H1_v1_winnowmap_HIFI.sam
samtools sort Nettle_female_H1_v1_winnowmap_HIFI.sam --threads 22 > Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam.bai Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam
winnowmap -k 19 -W repetitive_k19_H1.txt -ax map-pb $genome_H2 $HIFI_reads > winnowmap/Nettle_female_H2_v1_winnowmap_HIFI.sam
samtools sort Nettle_female_H2_v1_winnowmap_HIFI.sam --threads 22 > Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam


#4. you can download the reference genome you want to visualize at this stage, making sure that the indexing is done properly.

#5. extract bam file by chromosome
for chrnum in {1..9}; 
do
  samtools view -b Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_0${chrnum} > Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam.bai Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam; 
done

for chrnum in {10..13}; 
do
  samtools view -b Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_${chrnum} > Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam.bai Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam; 
done

for chrnum in {1..9}; 
do
  samtools view -b Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_0${chrnum} > Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam; 
done

for chrnum in {10..13}; 
do
  samtools view -b Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_${chrnum} > Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam; 
done

#6. prepare a bed file of inversion positions
nano Nettle_female_H1_syriINV_chr02.bed
cat syri.out | awk '{if ($11 == "INV"){print}}' | grep "Urtica_dioica_female_chr_02" | awk '{if ($3-$2 >= 10000){print}}' | cut -f -3
Urtica_dioica_female_chr_02     13675776        13687734
Urtica_dioica_female_chr_02     19833948        19844872
Urtica_dioica_female_chr_02     20209329        20227555
Urtica_dioica_female_chr_02     20228137        21550767
Urtica_dioica_female_chr_02     21721136        22012180
Urtica_dioica_female_chr_02     40088164        43483608
#also if you want, you can create a table for INV >10,000bp like this: 
cat syri.out | awk '{if ($11 == "INV"){print}}' | grep "Urtica_dioica_female_chr_02" | awk '{if ($3-$2 >= 10000){print}}'
Urtica_dioica_female_chr_02     13675776        13687734        -       -       Urtica_dioica_female_chr_02     13685943        13697901  INV547   -       INV     -
Urtica_dioica_female_chr_02     19833948        19844872        -       -       Urtica_dioica_female_chr_02     19126056        19134042  INV549   -       INV     -
Urtica_dioica_female_chr_02     20209329        20227555        -       -       Urtica_dioica_female_chr_02     19163924        19185438  INV550   -       INV     -
Urtica_dioica_female_chr_02     20228137        21550767        -       -       Urtica_dioica_female_chr_02     19476892        20954093  INV551   -       INV     -
Urtica_dioica_female_chr_02     21721136        22012180        -       -       Urtica_dioica_female_chr_02     21254530        21466722  INV552   -       INV     -
Urtica_dioica_female_chr_02     40088164        43483608        -       -       Urtica_dioica_female_chr_02     40601525        44371308  INV555   -       INV     -

#7. load the genome.fasta, then INV.bam and INV.bed on IGV. 


