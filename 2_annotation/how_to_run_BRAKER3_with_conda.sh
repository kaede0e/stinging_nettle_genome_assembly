#!/bin/bash
####### Run BRAKER3 pipeline in salal ########

# Miniforge is wonderful! https://github.com/conda-forge/miniforge
# 1. install with curl and .sh script
# 2. you need to
source /home/khirabayashi/miniforge3/bin/activate
(base) #this indicates that you are in the mamba/conda module basically

# my conda environnments are found in /home/khirabayashi/miniforge3/envs

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/khirabayashi/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/khirabayashi/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/home/khirabayashi/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/home/khirabayashi/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/home/khirabayashi/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "/home/khirabayashi/miniforge3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<

# I have added the above to activate miniforge3 from .bashrc.

# In addition, BRAKER3 requires many of its dependencies to be explicitly stated in .bashrc. Add below also:
export PATH=$PATH:/home/khirabayashi/bin/BRAKER/scripts
export PATH=$PATH:/home/khirabayashi/bin/BRAKER/GeneMark-ETP/bin
export PATH=$PATH:/home/khirabayashi/bin/BRAKER/GeneMark-ETP/tools
export GENEMARK_PATH=/home/khirabayashi/bin/BRAKER/GeneMark-ETP/bin
export AUGUSTUS_CONFIG_PATH=/home/khirabayashi/miniforge3/envs/augustus/config/
export AUGUSTUS_SCRIPTS_PATH=/home/khirabayashi/miniforge3/envs/augustus/bin/
export AUGUSTUS_BIN_PATH=/home/khirabayashi/miniforge3/envs/augustus/bin/
export PATH=$PATH:/home/khirabayashi/miniforge3/envs/augustus/config/
export PATH=$PATH:/home/khirabayashi/miniforge3/envs/augustus/bin/
export BAMTOOLS_PATH=/home/khirabayashi/miniforge3/envs/augustus/include/bamtools/api/
export PROTHINT_PATH=/home/khirabayashi/bin/BRAKER/GeneMark-ETP/bin/gmes/ProtHint/bin/
export PATH=$PATH:/home/khirabayashi/bin/BRAKER/GeneMark-ETP/bin/gmes/TSEBRA/bin/
export TSEBRA_PATH=/home/khirabayashi/bin/BRAKER/GeneMark-ETP/bin/gmes/TSEBRA/bin/

# exit the screen, open a new screen and you're in the mamba container
# or you can do source ~/.bashrc to activate the changes in the current screen.

# install all dependencies in appropriate conda/anaconda/bioconda/mamba settings - perl module 2.26 was sufficient for most programs EXCEPT Augustus (2.5.0) which required perl >2.32.
#So for this, I made a separate conda for Augustus and configured it there, and then add that configurations to the PATH.

# now you use BRAKER3 like this:

conda activate BRAKER3_env
(BRAKER3_env) # indicates you're in the BRAKER3_env container.

#Test dataset example:
cd ~/bin/BRAKER/example/tests/test3.sh
cat test3.sh:

#/-----------------------------------------------------
wd=test3

if [ -d $wd ]; then
    rm -r $wd
fi

# The expected runtime of this test is ~20 minutes.

# --gm_max_intergenic 10000 option is used here only to make the test run faster.
# It is not recommended to use this option in real BRAKER runs. The speed increase
# achieved by adjusting this option is negligible on full-sized genomes. Also,
# --skipOptimize is used here only to make the test run faster.
# Use an appropriate BUSCO lineage in real life use case. eukaryota_odb10 is used here
# only to make the test run faster.

# For instructions on how to prepare the proteins.fa input file from OrthoDB,
# see https://github.com/gatech-genemark/ProtHint#protein-database-preparation

( time braker.pl --genome=../genome.fa --prot_seq=../proteins.fa --bam=../RNAseq.bam --workingdir=$wd --threads 8 --skipOptimize ) &> test3.log
#/-----------------------------------------------------

#BRAKER3 run:
# pre-requisite files
#1. softmask your genome (I used RED - Repeat Detector in Cedar)
#2. align RNAseq to your genome and export filel in BAM format
#3. ensure desired protein database from OrthoDB is installed
# then run braker (This is the test3.sh pipeline that uses RNAseq+protein search and TSEBRA to filter gene hits):

#braker.pl --PROTHINT_PATH=/home/ProtHint/bin --species "stinging_nettle" --genome=Nettle_female_H1_v1_genome.softmasked.fa --prot_seq=Viridiplantae.fa --bam=Nettle_female_H1_v1_genome_red_softmasked_RNAseq_on_genome_hisat2_aln.bam --softmasking --workingdir=braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB --threads 20 &

braker.pl --genome=BRAKER3_input/Nettle_female_H1_Round_5_chrname_reordered_genome_repeatmsaked.fa.softmasked.fa \
--species "Urtica_dioica" \
--prot_seq=BRAKER3_input/Viridiplantae.fa \
--bam=BRAKER3_input/Nettle_female_H1_RNAseq_on_genome_hisat2_aln.bam \
--softmasking \
--workingdir=braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB --threads 20 &> Nettle_female_H1_braker3_v1.log

braker.pl --genome=BRAKER3_input/Nettle_female_H2_Round_5_chrname_reordered_genome_repeatmsaked.fa.softmasked.fa \
--species "Urtica_dioica" \
--useexisting \
--prot_seq=BRAKER3_input/Viridiplantae.fa \
--bam=BRAKER3_input/Nettle_female_H2_RNAseq_on_genome_hisat2_aln.bam \
--softmasking \
--workingdir=braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB --threads 20 &> Nettle_female_H2_braker3_v1.log
