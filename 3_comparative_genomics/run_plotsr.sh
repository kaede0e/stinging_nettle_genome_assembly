#!/bin/bash
## I also installed plotsr in the same conda environment as syri /home/khirabayashi/miniforge3/envs. So, you have to do: 
#conda deactivate #if you're in another conda. 
#conda activate syri_env

plotsr --sr Nettle_male_H1_vs_H2_Round_3_genome_reoriented_chrrenamed_for_syri_minimap2_syri.out \
--genomes genomes.txt \
--cfg plotsr_custom.cfg \
-o syri_plotsr_output_plot.png -S 0.5 -W 7 -H 10 -f 8

### Customization of syri plotsr plot using the cfg file that looks like this:
#----------- plotsr_custom.cfg -----------#
```
## COLOURS and transparency for alignments (syntenic, inverted, translocated, and duplicated)
syncol:#CCCCCC
invcol:#AA4499
tracol:#9ACD32
dupcol:#00BBFF
synlwd:0                ## Line width for syntenic annotations
invlwd:0.1              ## Line width for inversions
tralwd:0.1              ## Line width for translocations
duplwd:0.1              ## Line width for duplications
alpha:0.8

## Margins and dimensions:
chrmar:0.1              ## Adjusts the gap between chromosomes and tracks. Higher values leads to more gap
exmar:0.1               ## Extra margin at the top and bottom of plot area
marginchr:0.1           ## Margin between adjacent chromosomes when using --itx

## Legend
legend:T                ## To plot legend use T, use F to not plot legend
genlegcol:-1            ## Number of columns for genome legend, set -1 for automatic setup
bbox:0,1.01,0.5,0.3             ## [Left edge, bottom edge, width, height]
bbox_v:0,1.1,0.5,0.3    ## For vertical chromosomes (using -v option)
bboxmar:0.5             ## Margin between genome and annotation legends

## Tracks
norm:T                  ## For each chromosome, independently normalise the y-axis of tracks. Use T for normalising independently, and F to normalise based on max value across all chromosomes

## Axis
maxl:-1                 ## Manually set maximum chromosome position. Use `-1` for automatic selection. Does not work with --itx
genname:F               ## Write genome names adjacent to the chromosome (T) or not (F)
```

#----------- genomes.txt -----------#
```
#file   #name   #tags
Nettle_male_H1_Round_3_genome_reoriented_chrrenamed_for_syri.chr.fa     H1      lw:1.5;lc:#44AA99
Nettle_male_H2_Round_3_genome_reoriented_chrrenamed_for_syri.chr.fa     H2      lw:1.5;lc:#CC6677
```
