cat << EOF > SPP_repeats_H1.sh
#!/bin/bash
#SBATCH --account=def-mtodesco
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=128000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i NettleFemaleR5 -f Nettle_female_H1_Round_5_chrname_reordered_genome.chr.fa -h H1 -c 15 -m 128000 -g FALSE &> RepeatOBserver_Hap1_R5_log.txt

EOF

sbatch SPP_repeats_H1.sh

#- using repeat_seq_finder.sh to find sequences of the repeats found by RepeatOBserver 
module load StdEnv/2020 seqkit/2.3.1 emboss/6.6.0
#telomeric repeat: 
./repeat_seq_finder.sh "/home/kaedeh/scratch/Nettle/repeat_annotation/RepeatObserver/output/input_chromosomes/NettleFemalev2_H1-AT/chromosome_files" \
"NettleFemalev2_H1-AT_Chr6part01.fasta" 100 100000 175 60 200 > telomeric_repeat_175bp_chr06.txt
