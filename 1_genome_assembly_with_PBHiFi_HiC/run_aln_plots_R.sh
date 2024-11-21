#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH --account=
#SBATCH --output=
#SBATCH --error=

module load r/4.4

### conver the aln.paf from minimap2 into a visual plot in ggplot2 from R script ###

#Step3 R plots

for i in $(ls *_tab.txt)
do
#Number_ch=`printf %02d $i`

Rscript aln_plots.R $i

done
