### HiFiasm + Hi-C genome assembly pipeline ###
This page shows steps taken to with commandilne tools in chronological order;j pipeline for a chromosome-scale phased genome assembly of stinging nettle (_Urtica dioica_).
To replicate this pipeline in your genome using your cluster, see the individual script separated by types of jobs/programs. I tried to explicitly label which scripts were used. 

### ---------Useful code for pre-assembly with Hi-Fi data----------- ###
```
#!/bin/bash
#SBATCH --account=
#SBATCH --time=1-0
#SBATCH --mem=30Gb
#SBATCH --output=
#SBATCH --error=

SMRT_ROOT=/home/~bin/pacbio/smrtlink
export PATH=$PATH:/home/~bin/pacbio/smrtlink/smrtcmds/bin/

bam2fastq -o Pacbio_hifi/Nettle_female_Pacbio_hifi Pacbio_hifi/Urtica_QC_m84074_230923_140523_s2.hifi_reads.bam

#to check that this conversion was complete 
zgrep . Nettle_female_Pacbio_hifi.fastq.gz | awk 'NR%4==2{c++; l+=length($0)}
END{
	print "Number of reads: "c;
	print "Number of bases in reads: "l
}'

#Number of reads: 4115159
#Number of bases in reads: 66308008222

#- for general info on a bam file: 
samtools view -H Urtica_QC_m84074_230923_140523_s2.hifi_reads.bam
```
#output: 

#nettle_f
```
@HD     VN:1.6  SO:unknown      pb:5.0.0
@RG     ID:e50f14c1/70--70      PL:PACBIO       DS:READTYPE=CCS;Ipd:Frames=ip;PulseWidth:Frames=pw;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;BASECALLERVERSION=5.0;FRAMERATEHZ=100.000000;BarcodeFile=metadata/m84074_230923_140523_s2.barcodes.fasta;BarcodeHash=4ce24eed41885057053eb1865e933bdd;BarcodeCount=96;BarcodeMode=Symmetric;BarcodeQuality=Score        LB:FPAC230363934-1A     PU:m84074_230923_140523_s2      SM:Bio Sample 71        PM:REVIO        BC:GACGTGTATGTGTGAG CM:R/P1-C1/5.0-25M
@PG     ID:ccs  PN:ccs  VN:7.0.0 (commit v7.0.0)        DS:Generate circular consensus sequences (ccs) from subreads.       CL:/opt/pacbio/tag-ccs-current/bin/ccs --streamed --log-level INFO --stderr-json-log --kestrel-files-layout --movie-name m84074_230923_140523_s2 --log-file metadata/m84074_230923_140523_s2.ccs.log --min-rq 0.9 --non-hifi-prefix fail --knrt-ada --pbdc-model /opt/pacbio/tag-ccs-current/bin/../models/revio_v1.onnx --alarms metadata/m84074_230923_140523_s2.ccs.alarms.json
@PG     ID:lima VN:2.7.1 (commit v2.7.1-1-gf067520)     CL:/opt/pacbio/tag-lima-current/bin/lima --movie-name m84074_230923_140523_s2 --kestrel-files-layout --quality hifi --output-missing-pairs --shared-prefix --hifi-preset SYMMETRIC-ADAPTERS --store-unbarcoded --split-named --reuse-source-uuid --reuse-biosample-uuids --stderr-json-log --alarms metadata/m84074_230923_140523_s2.hifi_reads.lima.alarms.json --log-file metadata/m84074_230923_140523_s2.hifi_reads.lima.log pb_formats/m84074_230923_140523_s2.hifi_reads.consensusreadset.primrose.xml metadata/m84074_230923_140523_s2.barcodes.fasta hifi_reads/m84074_230923_140523_s2.hifi_reads.demux.bam
@PG     ID:primrose     VN:1.4.0 (commit v1.4.0)        CL:/opt/pacbio/tag-primrose-current/bin/primrose --movie-name m84074_230923_140523_s2 --kestrel-files-layout --quality hifi --reuse-source-uuid --stderr-json-log --log-file metadata/m84074_230923_140523_s2.hifi_reads.primrose.log --alarms metadata/m84074_230923_140523_s2.hifi_reads.primrose.alarms.json
```

#nettle_m
```
@HD     VN:1.6  SO:coordinate   pb:5.0.0
@RG     ID:27e97fc0     PL:PACBIO       DS:READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;BASECALLERVERSION=5.0;SMRTCELLKIT=102-202-200;SMRTCELLID=EA074655;RUNID=r84185_20240229_221539;ICSVERSION=13.0.0.212033;MOVIELENGTH=1800.0;FRAMERATEHZ=100.000000   LB:pla2144033   PU:m84185_240301_023219_s3      SM:F140575      PM:REVIO        CM:R/P1-C1/5.0-25M
@PG     ID:ccs  PN:ccs  VN:8.0.0 (commit v8.0.0)        DS:Generate circular consensus sequences (ccs) from subreads.   CL:/opt/pacbio/tag-ccs-current/bin/ccs --streamed --log-level INFO --log-file metadata/m84185_240301_023219_s3.ccs.log --stderr-json-log --instrument-files-layout --movie-name m84185_240301_023219_s3 --min-rq -1 --non-hifi-prefix fail --knrt-ada --pbdc-model /opt/pacbio/tag-ccs-current/bin/../models/revio_13_0.onnx --alarms metadata/m84185_240301_023219_s3.ccs.alarms.json
@PG     ID:pbtrim       VN:1.0.0 (commit v1.0.0)        CL:/opt/pacbio/tag-trim-current/trim --streamed --movie-name m84185_240301_023219_s3 --instrument-files-layout --reuse-source-uuid --quality hifi --log-level INFO --log-file metadata/m84185_240301_023219_s3.hifi_reads.trim.log --alarms metadata/m84185_240301_023219_s3.hifi_reads.trim.alarms.json
@PG     ID:jasmine      PN:jasmine      VN:2.2.0 (commit v2.2.0)        DS:Predict 5mC in PacBio HiFi reads.    CL:/opt/pacbio/tag-jasmine-current/bin/jasmine --streamed --movie-name m84185_240301_023219_s3 --instrument-files-layout --quality hifi --reuse-source-uuid --stderr-json-log --gpu-devices 0,0,0,0 --log-level INFO --log-file metadata/m84185_240301_023219_s3.hifi_reads.jasmine.log --alarms metadata/m84185_240301_023219_s3.hifi_reads.jasmine.alarms.json
@PG     ID:samtools     PN:samtools     PP:jasmine      VN:1.17 CL:samtools view -H m84185_240301_023219_s3.hifi_reads.bam
```

#SMRT tools commandline in Cedar 
#- gives you basic information on your HiFi reads (read length average, Qscore, yield in bp)
```
SMRT_ROOT=/home/~bin/pacbio/smrtlink
export PATH=$PATH:/home/~bin/pacbio/smrtlink/smrtcmds/bin/
dataset create --name Nettle_female_Pacbio_hifi --type ConsensusReadSet Nettle_female_Pacbio_hifi.xml Urtica_QC_m84074_230923_140523_s2.hifi_reads.bam
mkdir Nettle_female_Pacbio_hifi_runqc
runqc-reports -o Nettle_female_Pacbio_hifi_runqc Nettle_female_Pacbio_hifi.xml
```

#------------------------------------------------------------------------#
```
#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --mem=192000M
#SBATCH --account=
#SBATCH --cpus-per-task=48
#SBATCH --output=
#SBATCH --error=

### Draft genome assembly with HiFiasm + Hi-C data ###

#####################################
### Execution of programs ###########
#####################################

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

module load hifiasm/0.19.5

#####################################
##### Variables / data ##############
#####################################

HiC_1=/home/kaedeh/scratch/Nettle/HiC/HiC_Nettle_FemQC1_FKDL240020059-1A_223GFKLT4_L1_1.fq.gz
HiC_2=/home/kaedeh/scratch/Nettle/HiC/HiC_Nettle_FemQC1_FKDL240020059-1A_223GFKLT4_L1_2.fq.gz
Hifi_fastq=/home/kaedeh/scratch/Nettle/Pacbio_hifi/Nettle_female_Pacbio_hifi.fastq.gz

hifiasm -o Nettle_female.asm -t 48 --h1 $HiC_1 --h2 $HiC_2 $Hifi_fastq -s 0.4 --hom-cov 128

# ---------------------------------------------------------------------
echo "Done assembly with Hifiasm. Use 3D-DNA to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
```

# ---------------------------------------------------------------------
## Juicer and 3D-DNA pipelne for genome scaffolding ## 
#1. Run juicer to produce .hic and .assembly
#2. Open the .hic file in Juicebox (can use cloud version for bigger files)
#3. 3d-dna first step keeps all contigs, 2nd and 3rd step polishes - removes contigs to trash (on debris. 
#4. For placing cnntigs in the right place, look for lines of red/white. If there is a thick white line, move contigs. 

#do this to run it as a pipeline based on Eric's script (https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/Juicer)
#haplotype 1 
ln -s /home/~bin/juicer/CPU/ scripts
cd scripts/common
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
cd ../..
mkdir fastq
cd fastq
ln -s /home/Nettle/raw_data/HiC/HiC_Nettle_FemQC1_FKDL240020059-1A_223GFKLT4_L1_1.fq.gz HiC_Nettle_FemQC1_R1.fastq.gz
ln -s /home/Nettle/raw_data/HiC/HiC_Nettle_FemQC1_FKDL240020059-1A_223GFKLT4_L1_2.fq.gz HiC_Nettle_FemQC1_R2.fastq.gz
cd ../..
mkdir references
cd references/
ln -s /home/Nettle/output/assembly/Nettle_female.asm.hic.hap1.p_ctg.fa Nettle_female.asm.hic.hap1.p_ctg.fa
module load bwa
bwa index Nettle_female.asm.hic.hap1.p_ctg.fa #Run this in a job
cd ..
mkdir restriction_sites 
cd restriction_sites
scp /home/~bin/juicer/misc/generate_site_positions.py .
#---- here download the juicer/misc/generate_site_positions.py and edit accordingly
  filenames = {
    'hg19': '/seq/references/Homo_sapiens_assembly19.fasta',
    'mm9' : '/seq/references/Mus_musculus_assembly9.fasta',
    'mm10': '/seq/references/Mus_musculus_assembly10.fasta',
    'hg18': '/seq/references/Homo_sapiens_assembly18.fasta',
    'Nettle_female.asm.hic.hap1.p_ctg': '../references/Nettle_female.asm.hic.hap1.p_ctg.fa', #here you put your contig/assembly and its path
  }
module load python
python generate_site_positions.py DpnII Nettle_female.asm.hic.hap1_hap2.p_chr #Run this in a job

#generate a file Chromosome_sizes.sh with this:
for i in $(ls *DpnII.txt)
do
name=$(echo $i | cut -d "." -f 1 )
awk 'BEGIN{OFS="\t"}{print $1, $NF}'  $i > "$name"".chrom.sizes"
done

#---- now run this script in a core with many CPUs ----------#
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=510000M
#SBATCH --cpus-per-node=32
#SBATCH --account=
#SBATCH --output=
#SBATCH --error=

#####################################
### Execution of programs ###########
#####################################

module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
export PATH=$PATH:/home/~bin/juicer/CPU

#run juicer
bash scripts/juicer.sh -g Nettle_female.asm.hic.hap1.p_ctg -s DpnII \
-z references/Nettle_female.asm.hic.hap1.p_ctg.fa -y restriction_sites/Nettle_female.asm.hic.hap1.p_ctg_DpnII.txt \
-p restriction_sites/Nettle_female.chrom.sizes \
-t 28 -S early

# ---------------------------------------------------------------------
echo "Done Juicer Hi-C analysis.  Use 3D-DNA to scaffold contigs further."
# ---------------------------------------------------------------------

echo "Finished job at `date`"

#job efficiency: 
Nodes: 1
Cores per node: 32
CPU Utilized: 3-16:30:59
CPU Efficiency: 35.76% of 10-07:33:20 core-walltime
Job Wall-clock time: 07:44:10
Memory Utilized: 61.89 GB
Memory Efficiency: 49.51% of 125.00 GB
# -------------------------------------------------------------------------


#---- next we run the 3D-DNA pipeline to visualize and fix misassemblies ------------#
#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH --account=
#SBATCH --output=
#SBATCH --error=

### 3D-DNA pipeline after getting Hi-C contact merge_nodups.txt file ### 

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

#####################################
### Execution of programs ###########
#####################################

module load StdEnv/2020 python/3.11.2 java/17.0.2 lastz/1.04.03
#virtualenv 3ddna_env
#pip install scipy numpy matplotlib #libraries required for 3d-dna 
source /home/~bin/3ddna_env/bin/activate
export PATH=$PATH:/home/~bin/3d-dna #3D de novo assembly: version 180114


#####################################
##### Variables / data ##############
#####################################
contig_fasta=/home/Nettle/HiC_hap1/references/*.fa
merged_nodups=/home/Nettle/HiC_hap1/aligned/merged_nodups.txt

#run script; -r 0 runs only the first scaffolding, no polishing. 
run-asm-pipeline.sh -m haploid -r 0 $contig_fasta $merged_nodups

deactivate


# ---------------------------------------------------------------------
echo "Done 3D-DNA pipeline scaffolding.  Use Juicebox to visualize and manually fix scaffolds."
# ---------------------------------------------------------------------

echo "Finished job at `date`"


##### Convert .assembly to .fasta with R script by Moshutava and Eric ######

#!/bin/bash
#SBATCH --account=
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --output=
#SBATCH --error=

module load r

mkdir FASTA
Rscript asm_to_fasta_me.R Nettle_female.asm.hic.hap1.p_ctg.0.review2.assembly Round_1 Nettle_female.asm.hic.hap1.p_ctg.fa FASTA

#----- Here, your assembly is almost ready. Use the chromosome-level assemblies to perform a crude gene annotation. (see my RNAseq based TransDecoder pipeline) -----# 




#----- To check if SVs between haplotypes are real ---------# 

1) (after annotation) run anchorwave
2) run minimap2 alignment + SyRI, check for those regions (using the output coordinates) on Hi-C map
3) Visualize if long-reads map to breakpoints of SVs

#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=20G
#SBATCH --output=
#SBATCH --error=


#----------- Anchorwave whole genome alignment from chromosome-level assembly data -----------#

# Anchorwave - whole genome alignment from chromosome-level assembly data
# This software allows you to align two genomes using coding sequence (CDS),
# and the output can be tabulated and exported, visualized on R.

mopdule load StdEnv/2023 minimap2/2.28
export PATH=$PATH:/home/~bin/anchorwave/
export PATH=$PATH:/home/~bin/anchorwave/scripts

#0 Set variables
Genome=Nettle_female
Hap1=/home/kaedeh/scratch/Nettle/annotation/final_annotation_files/Round_2_hap1.reviewed.chr_assembled.fasta
Hap1_gene_anno_gff3=/home/kaedeh/scratch/Nettle/annotation/final_annotation_files/Nettle_female_hap1_genes_anno.gff3
Hap2=/home/kaedeh/scratch/Nettle/HiC_hap2/3d_dna_pipeline/FASTA/Round_2_hap2.reviewed.chr_assembled.fasta

#1 extracting CDS on one haplotype
anchorwave gff2seq \
-i $Hap1_gene_anno_gff3 -r $Hap1 -o anchorwave_gff2seq_${Genome}_H1.cds.fasta

#2 aligning CDS on both haplotypes
minimap2 -x splice -t 4 -k 12 -a -p 0.4 -N 20 \
$Hap1 anchorwave_gff2seq_${Genome}_H1.cds.fasta > anchorwave_minimap2_${Genome}_H1.cds.sam
minimap2 -x splice -t 4 -k 12 -a -p 0.4 -N 20 \
$Hap2 anchorwave_gff2seq_${Genome}_H1.cds.fasta > anchorwave_minimap2_${Genome}_H2.cds.sam

#3 (Not necessary) visualize the alignment to assess synteny on R by transforming SAM file to Tabulated file 
#perl /home/~bin/anchorwave/scripts/alignmentToDotplot.pl \
#$Hap1_gene_anno_gff3 anchorwave_minimap2_${Genome}_H1.cds.sam > anchorwave_minimap2_${Genome}_H1.cds.tab
#perl /home/~bin/anchorwave/scripts/alignmentToDotplot.pl \
#$Hap1_gene_anno_gff3 anchorwave_minimap2_${Genome}_H2.cds.sam > anchorwave_minimap2_${Genome}_H2.cds.tab

#4 genome alignment with relocation variation, chromosome fusion or whole genome duplication - proali
anchorwave proali \
-i $Hap1_gene_anno_gff3 \
-as anchorwave_gff2seq_${Genome}_H1.cds.fasta \
-r $Hap1 \
-a anchorwave_minimap2_${Genome}_H2.cds.sam \
-ar anchorwave_minimap2_${Genome}_H1.cds.sam \
-s $Hap2 \
-n ${Genome}_Hap1_vs_Hap2.anchors \
-R 1 -Q 1 -t 9 \
-o anchorwave_proali_${Genome}_Hap1_vs_Hap2.maf \
-f anchorwave_proali_${Genome}_Hap1_vs_Hap2.f.maf > anchorwave_proali_log.txt

#5 convert .maf to plot table
#!/bin/bash
#SBATCH --account=def-mtodesco
#SBATCH --time=24:00:00
#SBATCH --mem=50G

module load StdEnv/2020 python/2.7.18

export PATH=$PATH:/home/~bin/anchorwave/
export PATH=$PATH:/home/~bin/anchorwave/scripts

python2 maf-convert.py tab anchorwave_proali_Nettle_female_Hap1_vs_Hap2.maf > anchorwave_proali_Nettle_female_Hap1_vs_Hap2.tab

#6 convert table to format to plot
for i in $(ls *_Hap1_vs_Hap2.tab)
do
       name=$(echo $i | cut -d "." -f 1 )

       grep -v "#" $i | awk -F "\t" -v OFS="\t" '{print $3, $3 + $4, $8, $8 + $9, 30, $7, $2, $10}' > "$name""_tab_fp.txt"

done


#----- 2. To check SVs on Hi-C heatmap using precise coordinates (especially if the region doesn't fall between contigs and requires you to cut contigs for misassembly), run minimap2+SyRI ---------# 
cwd="."     # Change to working directory
cd $cwd
PATH_TO_SYRI="/home/khirabayashi/miniforge3/envs/syri_env/bin/syri" #Change the path to point to syri

## NOTE: IDEALLY, SYRI EXPECTS THAT THE HOMOLOGOUS CHROMOSOMES IN THE TWO GENOMES
## WOULD HAVE EXACTLY SAME CHROMOSOME ID. THEREFORE, IT IS HIGHLY RECOMMENDED THAT
## THE USER PRE-PROCESSES THE FASTA FILES TO ENSURE THAT HOMOLOGOUS CHROMOSOMES
## HAVE EXACTLY THE SAME ID IN BOTH FASTA FILES CORRESPONDING TO THE TWO GENOMES.
## IN CASE, THAT IS NOT THE CASE, SYRI WOULD TRY TO FIND HOMOLOGOUS GENOMES USING
## WHOLE GENOME ALIGNMENTS, BUT THAT METHOD IS HEURISTICAL AND CAN RESULT IN
## SUBOPTIMAL RESULTS.
## ALSO, IT IS REQUIRED THAT THE TWO GENOMES (FASTA FILES) HAVE THE SAME NUMBER
## OF CHROMOSOMES.

ln -sf GCA_964188135.1_drUrtDioi1.1_genomic.fna refgenome
ln -sf Nettle_female_H1_Round_5_chrname_reordered_genome.chr.fa qrygenome_H1
ln -sf Nettle_female_H2_Round_5_chrname_reordered_genome.chr.fa qrygenome_H2

## Perform whole genome alignment
# Using minimap2 for generating alignment. Any other whole genome alignment tool can also be used.
minimap2 -ax asm5 --eqx refgenome qrygenome_H1 > Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_aln.sam
minimap2 -ax asm5 --eqx refgenome qrygenome_H2 > Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_aln.sam

## Run SyRI round 1 with SAM or BAM file as input
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_aln.sam -r refgenome -q qrygenome_H1 --prefix DToLpri_Round5_H1 -k -F S 2> syri_error1_H1_out.txt
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_aln.sam -r refgenome -q qrygenome_H2 --prefix DToLpri_Round5_H2 -k -F S 2> syri_error1_H2_out.txt

## Note: SyRI requires your chromosomes to be aligned in the syntenic direction. If one is reversed, you need to reverse-complement it to proceed.
cat syri_error1_H1_out.txt | grep "ERROR" | cut -f 11 -d " " | sed s/"\."/""/g > chr_to_rev_H1.txt ;
mkdir qrygenome_H1_chr
cat Nettle_female_H2_Round_5_syri_input_genome_chrnames.txt | grep -v -f chr_to_rev_H1.txt > chr_to_keep_H1.txt
  for chr in `cat chr_to_keep_H1.txt` ;
  do
    echo $chr > qrygenome_H1_chr/$chr.txt ;
    samtools faidx qrygenome_H1 -r qrygenome_H1_chr/$chr.txt > qrygenome_H1_chr/$chr.fwd.fa ;
  done #forward strand fasta for each chromosome
  for chr in `cat chr_to_rev_H1.txt` ;
  do
    echo $chr > qrygenome_H1_chr/rev_$chr.txt ;
    samtools faidx qrygenome_H1 -i -r qrygenome_H1_chr/rev_$chr.txt > qrygenome_H1_chr/$chr.rev.fa ;
  done
  cat qrygenome_H1_chr/*.fa > qrygenome_H1_chr_rev.fa
cat qrygenome_H1_chr_rev.fa | sed s/'\/rc'/''/g > qrygenome_H1_chr_rev_renamed.fa

cat syri_error1_H2_out.txt | grep "ERROR" | cut -f 11 -d " " | sed s/"\."/""/g > chr_to_rev_H2.txt ;
mkdir qrygenome_H2_chr
cat Nettle_female_H2_Round_5_syri_input_genome_chrnames.txt | grep -v -f chr_to_rev_H2.txt > chr_to_keep_H2.txt
  for chr in `cat chr_to_keep_H2.txt` ;
  do
    echo $chr > qrygenome_H2_chr/$chr.txt ;
    samtools faidx qrygenome_H2 -r qrygenome_H2_chr/$chr.txt > qrygenome_H2_chr/$chr.fwd.fa ;
  done #forward strand fasta for each chromosome
  for chr in `cat chr_to_rev_H2.txt` ;
  do
    echo $chr > qrygenome_H2_chr/rev_$chr.txt ;
    samtools faidx qrygenome_H2 -i -r qrygenome_H2_chr/rev_$chr.txt > qrygenome_H2_chr/$chr.rev.fa ;
  done
  cat qrygenome_H2_chr/*.fa > qrygenome_H2_chr_rev.fa
cat qrygenome_H2_chr_rev.fa | sed s/'\/rc'/''/g > qrygenome_H2_chr_rev_renamed.fa

###This qrygenome fasta has weirdly retained chromosome names from empty lines also and messed up the fasta.
###Manually inspect the fasta so that headers are all correctly in fasta format.

# Re-do whole genome alignment with the reverse complemented fasta
minimap2 -ax asm5 --eqx refgenome qrygenome_H1_chr_rev_renamed.fa -t 10 > Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_rev2_aln.sam
minimap2 -ax asm5 --eqx refgenome qrygenome_H2_chr_rev_renamed.fa -t 10 > Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_rev2_aln.sam

## Run SyRI round 2 with the new SAM file as input NOTE: SYRI doesn't take sam files with "." in in front of .sam
#python3 $PATH_TO_SYRI -c Nettle_female_Round_5_H1_H2_minimap2_rev2_aln.sam -r refgenome -q qrygenome_chr_rev_renamed.fa -k -F S 2> syri_error2_out.txt
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H1_minimap2_rev2_aln.sam -r refgenome -q qrygenome_H1_chr_rev_renamed.fa --prefix DToLpri_Round5_H1_rev -k -F S 2> syri_error2_H1_out.txt
python3 $PATH_TO_SYRI -c Nettle_female_DToL_primary_Round_5_reordered_H2_minimap2_rev2_aln.sam -r refgenome -q qrygenome_H2_chr_rev_renamed.fa --prefix DToLpri_Round5_H2_rev -k -F S 2> syri_error2_H2_out.txt

## SyRI would report genomic structural differences in syri.out and syri.vcf.

## Plot the syri.out using plotsr - also installed in the SyRI conda env.
#plotsr --sr syri.out --genomes genomes.txt -o syri_plotsr_output_plot.png -S 0.5 -W 7 -H 10 -f 8

## then from syri.out file, you can select for regions you want to check (eg. "INV") and copy those to Excel spreadsheet for record keeping. Go to Juicebox and check which orientations makes more sense.


#----- 3. To check if SVs between haplotypes are real, visualize HIFI reads on IGV ---------# 
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=24
#SBATCH --account=def-mtodesco
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_ONT_HIFI_winnowmap_aln_samtobam.13Sep2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1_ONT_HIFI_winnowmap_aln_samtobam.13Sep2024.err

module load samtools

export PATH=$PATH:/home/~bin/Winnowmap/bin
HIFI_reads=/home/kaedeh/scratch/Nettle/Pacbio_hifi/Nettle_female_Pacbio_hifi_Q30_filtered.fastq
ONT_reads=/home/kaedeh/scratch/Nettle/HiC_hap1_hap2/references/Nettle_female_canu_correctedReads_10kb.fq
genome_H1=/home/kaedeh/scratch/Nettle/HiC_hap1_hap2/references/round_5/Nettle_female_H1_Round_5_syri_input_genome.fa
genome_H2=/home/kaedeh/scratch/Nettle/HiC_hap1_hap2/references/round_5/Nettle_female_H2_Round_5_syri_input_genome.fa

#1. build a merylDB from your ONT/HIFI reads that you want to map
#meryl k=19 count output meryl_k19_nettle_female_H1.merylDB $genome_H1
#meryl k=19 count output meryl_k19_nettle_female_H2.merylDB $genome_H2

#2. prepare repeats for winnowmap
#meryl print greater-than distinct=0.9998 meryl_k19_nettle_female_H1.merylDB > repetitive_k19_H1.txt
#meryl print greater-than distinct=0.9998 meryl_k19_nettle_female_H2.merylDB > repetitive_k19_H2.txt

#3. run winnowmap to map HiFi and ONT reads to the assembly
#winnowmap -k 19 -W repetitive_k19_H1.txt -ax map-ont $genome_H1 $ONT_reads > winnowmap/Nettle_female_H1_v1_winnowmap_ONT.sam
samtools sort Nettle_female_H1_v1_winnowmap_ONT.sam --threads 22 > Nettle_female_H1_v1_winnowmap_ONT.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H1_v1_winnowmap_ONT.sorted.bam.bai Nettle_female_H1_v1_winnowmap_ONT.sorted.bam
#winnowmap -k 19 -W repetitive_k19_H2.txt -ax map-ont $genome_H2 $ONT_reads > winnowmap/Nettle_female_H2_v1_winnowmap_ONT.sam
samtools sort Nettle_female_H2_v1_winnowmap_ONT.sam --threads 22 > Nettle_female_H2_v1_winnowmap_ONT.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H2_v1_winnowmap_ONT.sorted.bam.bai Nettle_female_H2_v1_winnowmap_ONT.sorted.bam
#winnowmap -k 19 -W repetitive_k19_H1.txt -ax map-pb $genome_H1 $HIFI_reads > winnowmap/Nettle_female_H1_v1_winnowmap_HIFI.sam
samtools sort Nettle_female_H1_v1_winnowmap_HIFI.sam --threads 22 > Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam.bai Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam
#winnowmap -k 19 -W repetitive_k19_H1.txt -ax map-pb $genome_H2 $HIFI_reads > winnowmap/Nettle_female_H2_v1_winnowmap_HIFI.sam
samtools sort Nettle_female_H2_v1_winnowmap_HIFI.sam --threads 22 > Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam
samtools index -M -b --threads 22 -o Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam

#4. you can download the reference genome you want to visualize at this stage, making sure that the indexing is done properly.

#5. extract bam file by chromosome
for chrnum in {1..9}; 
do
  samtools view -b Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_0${chrnum} > Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam.bai Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam; 
done

for chrnum in {10..13}; 
do
  samtools view -b Nettle_female_H1_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_${chrnum} > Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam.bai Nettle_female_H1_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam; 
done

for chrnum in {1..9}; 
do
  samtools view -b Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_0${chrnum} > Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_0${chrnum}.bam; 
done

for chrnum in {10..13}; 
do
  samtools view -b Nettle_female_H2_v1_winnowmap_HIFI.sorted.bam Urtica_dioica_female_chr_${chrnum} > Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam;
  samtools index -M -b -o Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam.bai Nettle_female_H2_v1_winnowmap_HIFI.sorted_chr_${chrnum}.bam; 
done

#6. prepare a bed file of inversion positions
nano Nettle_female_H1_syriINV_chr02.bed
cat syri.out | awk '{if ($11 == "INV"){print}}' | grep "Urtica_dioica_female_chr_02" | awk '{if ($3-$2 >= 10000){print}}' | cut -f -3
Urtica_dioica_female_chr_02     13675776        13687734
Urtica_dioica_female_chr_02     19833948        19844872
Urtica_dioica_female_chr_02     20209329        20227555
Urtica_dioica_female_chr_02     20228137        21550767
Urtica_dioica_female_chr_02     21721136        22012180
Urtica_dioica_female_chr_02     40088164        43483608
#also if you want, you can create a table for INV >10,000bp like this: 
cat syri.out | awk '{if ($11 == "INV"){print}}' | grep "Urtica_dioica_female_chr_02" | awk '{if ($3-$2 >= 10000){print}}'
Urtica_dioica_female_chr_02     13675776        13687734        -       -       Urtica_dioica_female_chr_02     13685943        13697901  INV547   -       INV     -
Urtica_dioica_female_chr_02     19833948        19844872        -       -       Urtica_dioica_female_chr_02     19126056        19134042  INV549   -       INV     -
Urtica_dioica_female_chr_02     20209329        20227555        -       -       Urtica_dioica_female_chr_02     19163924        19185438  INV550   -       INV     -
Urtica_dioica_female_chr_02     20228137        21550767        -       -       Urtica_dioica_female_chr_02     19476892        20954093  INV551   -       INV     -
Urtica_dioica_female_chr_02     21721136        22012180        -       -       Urtica_dioica_female_chr_02     21254530        21466722  INV552   -       INV     -
Urtica_dioica_female_chr_02     40088164        43483608        -       -       Urtica_dioica_female_chr_02     40601525        44371308  INV555   -       INV     -

#7. load the genome.fasta, then INV.bam and INV.bed on IGV. 


#----- Finally, combine ONT reads with TGS GapCloser to see if we can close any gaps in the assembly ---------# 

#!/bin/bash
#SBATCH --account=
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=250Gb
#SBATCH --output=
#SBATCH --error=

#0. before starting, check ONT reads yield and quality and filter as necessary
#also the input ONT reads has to be in .fasta format so convert that. 
#reformat.sh in=nettle_basecalled_reads_combined_53122.fastq \
#lhist=Nettle_female_ONT_length_histogram aqhist=Nettle_female_ONT_avgquality_histogram

module load racon/1.5.0

export PATH=$PATH:/home/~bin/TGS-GapCloser
#TGS_READS_FILE=
#tgsgapcloser --scaff SCAFF_FILE --reads TGS_READS_FILE --output OUT_PREFIX

#mkdir TGSgapclosing_hap1 #you have to start with a fresh new directory and run the script from here. 
tgsgapcloser --scaff ../Round_2_hap2.reviewed.chr_assembled.fasta --reads ../Nettle_female_canu_correctedReads_30kb.fq --output ../Round_2_hap2_canu_corrected_ONT_30kbp_TGSgapcloser --ne --thread 16

#well.... I did this and run the program but it wasn't able to patch any gaps with the ONT reads I gave. 

