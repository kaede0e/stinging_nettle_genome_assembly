#!/bin/bash
#SBATCH --account=
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=17
#SBATCH --mem=30Gb
#SBATCH --output=
#SBATCH --error=

#0. before starting, check ONT reads yield and quality and filter as necessary
#also the input ONT reads has to be in .fasta format so convert that.
module load StdEnv/2023 bbmap/39.06
#reformat.sh in=Nettle_female_canu_correctedReads.fasta \
#minlength=30000 out=Nettle_female_canu_correctedReads_30kb.fq
#lhist=Nettle_female_ONT_length_histogram aqhist=Nettle_female_ONT_avgquality_histogram

module load racon/1.5.0

export PATH=$PATH:/home/~/bin/TGS-GapCloser
#TGS_READS_FILE=
#tgsgapcloser --scaff SCAFF_FILE --reads TGS_READS_FILE --output OUT_PREFIX
tgsgapcloser  \
        --scaff  Round_2_hap1.reviewed.chr_assembled.fasta \
        --reads  Nettle_female_canu_correctedReads_30kb.fq \
#        --output canu_tgsgapcloser_hap1_30kb_longer_ne \
        --ne \
        --thread 16 \
        >pipe.log 2>pipe.err
