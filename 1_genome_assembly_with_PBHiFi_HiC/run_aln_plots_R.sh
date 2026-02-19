#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH --account=
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_hap2_minimap2_aln.23May2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_hap2_minimap2_aln.23May2024.err

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

module load StdEnv/2023 minimap2/2.28

#Steo1 run minimap2 alignment
minimap2 -cx asm5 Round_2_hap1.reviewed.chr_assembled.fasta Round_2_hap2.reviewed.chr_assembled.fasta > Nettle_female_hap1_hap2_aln.paf

#Step2 convert paf to table to plot and select alignments with MQ > 30
for i in $(ls *paf)
do
        name=$(echo $i | cut -d "." -f 1)
        awk '$12 >= 30 {print $8, $9, $3, $4, "30", $1, $6, $5 }' $i > "$name""MQ30_tab.txt"
done


#------------------------------------------# 
#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH --account=
#SBATCH --output=
#SBATCH --error=

module load r/4.4

### conver the aln.paf from minimap2 into a visual plot in ggplot2 from R script ###
# Make sure you have the aln_plots.R script from Eric copied to your directory. 

#Step3 R plots

for i in $(ls *_tab.txt)
do
#Number_ch=`printf %02d $i`

Rscript aln_plots.R $i

done
