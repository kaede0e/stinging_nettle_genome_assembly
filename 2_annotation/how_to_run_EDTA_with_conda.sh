#!/bin/bash
####### Conda & mamba in salal ########

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

# I have added the above to activate mamba from .bashrc. However, I probably should remove these lines...

# EDTA requirede mamba activate so I added mamba init which gives you the above.
# exit the screen, open a new screen and you're in the mamba container

#install EDTA: (original mamba install EDTA did't work for some reason so do this instead so you install all specific dependencideds.)
mamba create -n EDTA2.2 -c conda-forge -c bioconda -c r annosine2 biopython cd-hit coreutils genericrepeatfinder genometools-genometools glob2 h5py==3.9 keras==2.11 ltr_finder ltr_retriever mdust multiprocess muscle openjdk pandas perl perl-text-soundex pyarrow python r-base r-dplyr regex repeatmodeler r-ggplot2 r-here r-tidyr scikit-learn swifter tensorflow-cpu==2.11 tesorter samtools bedtools grep

# now you use EDTA like this:

mamba activate EDTA2.2
(EDTA2.2) # indicates you're in the EDTA2.2 container installed with mamba.

#Test dataset example:
perl EDTA.pl --genome test/genome.fa --cds test/genome.cds.fa --exclude test/genome.exclude.bed --sensitive 1 --anno 1 --threads 8

perl /home/khirabayashi/bin/EDTA/EDTA.pl \
--genome EDTA_input/EDTA_input_Nettle_female_H1_v1_genome.chr.fa \
--sensitive 1 --anno 1 --overwrite 0 --threads 32 &> EDTA_error2.txt

#Threads 8 took way too long to finish LINE search so I bumped it up to 32 threads (discuss with your lab mates for resource allocation) and finishsed the job in 2~3 days.
