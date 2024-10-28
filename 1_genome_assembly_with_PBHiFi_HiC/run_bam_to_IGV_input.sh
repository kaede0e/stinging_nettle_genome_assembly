#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=rrg-rieseber-ac
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000M

module load samtools

#make sure bam files are sorted and indexed
#samtools view
for chrnum in {1..9};
do
  samtools view -b Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_0${chrnum} -q 60 > Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam;
done

for chrnum in {11..13};
do
  samtools view -b Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_${chrnum} -q 60 > Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam;
done
