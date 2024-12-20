#!/bin/bash
#SBATCH --account=
#SBATCH --time=24:00:00
#SBATCH --mem=50G

module load StdEnv/2020 python/2.7.18

export PATH=$PATH:/home/~/bin/anchorwave/
export PATH=$PATH:/home/~/bin/anchorwave/scripts

python2 maf-convert.py tab anchorwave_proali_Nettle_female_Hap1_vs_Hap2.maf > anchorwave_proali_Nettle_female_Hap1_vs_Hap2.tab
