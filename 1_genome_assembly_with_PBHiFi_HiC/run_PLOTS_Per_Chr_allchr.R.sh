#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000M

module load StdEnv/2023 r/4.4.0

for j in $( ls *tab_fp.txt)
do
   for i in $(cut -f 7 $j | sort | uniq)
   do
   Rscript PLOTS_Per_Chr_allchr.R $j $i
   done
done
