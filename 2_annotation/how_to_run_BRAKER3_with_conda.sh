#!/bin/bash
####### Run BRAKER3 pipeline in salal ########

# Miniforge is wonderful! https://github.com/conda-forge/miniforge
# 1. install with curl and .sh script
# 2. you need to
source /home/khirabayashi/miniforge3/bin/activate
(base) #this indicates that you are in the mamba/conda module basically

# my conda environnments are found in /home/khirabayashi/miniforge3/envs

# You need the following directory structure (in my case at least this worked)
#### .bashrc contains ####
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
#Installed BRAKER3 from github: /home/khirabayashi/bin/BRAKER
export PATH=$PATH:/home/khirabayashi/bin/BRAKER
export PATH=$PATH:/home/khirabayashi/bin/BRAKER/scripts

#Then install manually all the perl dependencies with conda install xxx (have to do one by one if any requirement files fail)
#
conda install --yes --file requirements.txt
#
#BRAKER3 requirements for my case: run the following as 
#"conda install -c bioconda --yes --file requirements_braker_with_bioconda.txt" 
echo "perl-app-cpanminus
perl-file-spec
perl-hash-merge
perl-list-util
perl-module-load-conditional
perl-posix
perl-file-homedir
perl-parallel-forkmanager
perl-scalar-util-numeric
perl-yaml
perl-class-data-inheritable
perl-exception-class
perl-test-pod
perl-file-which # skip if you are not comparing to reference annotation
perl-mce
perl-threaded
perl-list-util
perl-math-utils
cdbtools
perl-data-dumper" > requirements_braker_with_bioconda.txt

#"conda install -c --yes --file bioconda requirements_braker_with_anaconda.txt" 
echo "perl
biopython" > requirements_braker_with_bioconda.txt

#Augustus requirements for my case: run the following as
#"mamba install --yes --file requirements_augustus_mamba.txt"
echo "zlib
boost-cpp
gsl
lp_solve
biopython
perl
perl-app-cpanminus
perl-module-build
perl-dbi
perl-scalar-list-utils
perl-file-which
sqlite
suitesparse
tar" > requirements_augustus_mamba.txt
#"conda install -c --yes --file bioconda requirements_augustus_bioconda.txt"
echo "ucsc-twobitinfo
ucsc-fatotwobit
perl-yaml
perl-parallel-forkmanager
htslib
diamond
cdbtools
bamtools" > requirements_augustus_bioconda.txt


#now install rest of the programs that could not be found in the above list individually. 
conda install -c eumetsat perl-yaml-xs

#GeneMarkET downloaded from github inside BRAKER
export PATH=$PATH:/home/khirabayashi/bin/BRAKER/GeneMark-ETP/bin
export PATH=$PATH:/home/khirabayashi/bin/BRAKER/GeneMark-ETP/tools

#### check for dependencies with ./check_install.pl 
conda install hcc::perl-statistics-linefit #this was needed for GeneMarkETP , but problem
Could not solve for environment specs
The following packages are incompatible
├─ perl-statistics-linefit is installable and it requires
│  └─ perl >=5.26.2,<5.26.3.0a0 , which can be installed;
└─ perl 5.32.1  is not installable because it conflicts with any installable versions previously reported.

#as long as this check passes, you're probably okay to run the test script for BRAKER! 


### now you use BRAKER3 like this: ###

conda activate BRAKER3_env
(BRAKER3_env) # indicates you're in the BRAKER3_env container.

#Test dataset example:
cd ~/bin/BRAKER/example/tests/test3.sh
cat test3.sh:

#-----------------------------------------------------
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
#-----------------------------------------------------

#BRAKER3 run for nettle genome:
# pre-requisite files
#1. softmask your genome (I used RED - Repeat Detector)
```
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=10000
#SBATCH --output=
#SBATCH --error=

######### Genome repeat detection and soft-masking with RED ###########

#First, set up the working directory.
module load StdEnv/2020 gcc python/3.8 augustus hmmer blast+ metaeuk prodigal r
#virtualenv red_env
source /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/red_env/bin/activate
#pip install biopython
#pip install natsort

export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/RED_RepeatDetector/redUnix64
/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/redmask/redmask.py \
-i Nettle_male_hap1_v1.fa \
-o Nettle_male_hap1_v1

/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/redmask/redmask.py \
-i Nettle_male_hap2_v1.fa \
-o Nettle_male_hap2_v1
```

#2. align RNAseq to your genome and export filel in BAM format
```
#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=125000M
#SBATCH --output=
#SBATCH --error=

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

## Preparing bam file for RNAseq alignment to the refgenome to be used as BRAKER input ##

#####################################
#### Execution of programmes ########
#####################################

export PATH=$PATH:/home/kaedeh/scratch/Cannabis/annotation
export PATH=$PATH:/home/kaedeh/scratch/Cannabis/annotation/reference_RNAseq #store my Illumina reads here
module load StdEnv/2020 gcc/9.3.0 blast+/2.13.0
module load hisat2 stringtie samtools

# ---------------------------------------------------------------------

#### 0. QC RNA reads & check all files are there ########
# Make sure to run fastqc/0.11.9 or fastp before using the data.
# Make sure to copy one final draft genome to working directory
# Make sure to have all RNAseq .fastq files (paired end, stranded) in ~/raw_fastq_files/*_1.fastq or *_2.fastq, matching prefix for paired reads.

# ---------------------------------------------------------------------

input_refgenome=Nettle_male_hap2_v1.softmasked.fa
assembly_name=Nettle_male_hap2_v1.softmasked

#### Transcript alignment to genome ####
hisat2-build $input_refgenome ${assembly_name}_hisat2_index #took about 20min
ls /home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/reference_data/RNAseq/*_1_paired*.fastq.gz | tr '\n' ',' | sed 's/.$//' > RNAseq_m1_list_of_files.txt
ls /home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/reference_data/RNAseq/*_2_paired*.fastq.gz | tr '\n' ',' | sed 's/.$//' > RNAseq_m2_list_of_files.txt
hisat2 \
-q -x ${assembly_name}_hisat2_index -1 `cat RNAseq_m1_list_of_files.txt` -2 `cat RNAseq_m2_list_of_files.txt` -S ${assembly_name}_RNAseq_on_genome_hisat2_aln.sam --threads 12
samtools sort --threads 12 -o ${assembly_name}_RNAseq_on_genome_hisat2_aln.bam ${assembly_name}_RNAseq_on_genome_hisat2_aln.sam


#echo "Finished Hisat2 alignment: `date`"
```

#3. ensure desired protein database from OrthoDB is installed
# then run braker (This is the test3.sh pipeline that uses RNAseq+protein search and TSEBRA to filter gene hits):

#braker.pl --PROTHINT_PATH=/home/ProtHint/bin --species "stinging_nettle" --genome=Nettle_female_H1_v1_genome.softmasked.fa --prot_seq=Viridiplantae.fa --bam=Nettle_female_H1_v1_genome_red_softmasked_RNAseq_on_genome_hisat2_aln.bam --softmasking --workingdir=braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB --threads 20 &

braker.pl --genome=BRAKER3_input/Nettle_male_hap1_v1.softmasked.fa \
--species "Urtica_dioica_male" \
--prot_seq=BRAKER3_input/Viridiplantae.fa \
--bam=BRAKER3_input/Nettle_male_hap2_v1_RNAseq_on_genome_hisat2_aln.bam \
--softmasking \
--workingdir=braker3.Nettle_male_H1_v1_RNA_ProtsViridiplantaeOrthoDB --threads 20 &> Nettle_male_H1_braker3_v1.log

braker.pl --genome=BRAKER3_input/Nettle_male_hap2_v1.softmasked.fa \
--species "Urtica_dioica_male" \
--useexisting \
--prot_seq=BRAKER3_input/Viridiplantae.fa \
--bam=BRAKER3_input/Nettle_male_hap2_v1_RNAseq_on_genome_hisat2_aln.bam \
--softmasking \
--workingdir=braker3.Nettle_male_H2_v1_RNA_ProtsViridiplantaeOrthoDB --threads 20 &> Nettle_male_H2_braker3_v1.log
