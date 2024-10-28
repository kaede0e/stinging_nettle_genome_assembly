#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=24
#SBATCH --account=def-mtodesco
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_ONT_HIFI_winnowmap_aln_samtobam.13Sep2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_ONT_HIFI_winnowmap_aln_samtobam.13Sep2024.err

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

export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/Winnowmap/bin
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
