## HiFiasm + Hi-C genome assembly pipeline using YaHS + 3ddna (Feb 2025) ##
This page shows steps taken to with commandilne tools in chronological order; pipeline for a chromosome-scale phased genome assembly of stinging nettle (_Urtica dioica_) male.
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
## Juicer and YaHS + 3ddna for genome scaffolding ##
1. Run juicer to produce merged_nodups.txt file
2. Use Eric's awk script to convert merged_nodups.txt to a .bed file formatted for YaHS.
3. Run YaHS to scaffold contigs, and convert output file to .hic and .assembly.
4. Open the .hic and .assembly file in Juicebox and do manual curation per haplotype.
5. Export .reviewed.assembly and generate corrected fasta file. 
6. Combine H1 + H2 corrected fasta file and use this to run a new round of juicer using the same Hi-C reads.
7. From merged_nodups.txt, run 3ddna pipeline, first stage only (which only does scaffolding - no correction, and create .hic and .assembly file).
8. Open the .hic and .assembly file in Juicebox and do final manual curation with minimal cutting of contigs.
9. Generate fasta file from .review.assembly using Eric's script (asm_to_fasta.R).
10. Run juicer again on the final assembly, YaHS, then check on Juicebox per haplotype as a final review. 

#do this to run it as a pipeline based on Eric's script (https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/Juicer)
#haplotype 1 
```
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
```
Run BWA indexing in a job inside /references : 
```
#!/bin/bash
#SBATCH --account=
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb

module load bwa

bwa index *.fna
```
once indexing is done, continue with /restriction_sites
```
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
```

In the end, you want your directory structure to look like this:
```
/home/~/bin/juicer
	-- CPU/juicer.sh
 		/common
   			-- juicer_tools.jar -> juicer_tools.1.1.9_jcuda.0.8.ja
/scratch/HiC
	-- scripts -> /home/~/bin/juicer/CPU
 	-- /references
  		-- contig_assembly.fa -> ~/home/~/asm.fa
    		-- contig_assembly.fa.bwt
      		-- contig_assembly.fa.pac
    		-- contig_assembly.fa.ann
    		-- contig_assembly.fa.amb
      		-- contig_assembly.fa.sa
	-- /fastq
 		-- HiC_read_R1.fastq.gz -> ~/home/data/~/rawreads1.fq.gz
   		-- HiC_read_R2.fastq.gz -> ~/home/data/~/rawreads2.fq.gz
     	-- /restriction_sites 
      		-- generate_positions.py (edited)
		-- contig_assembly.DpnII.txt
  		-- conting_assembly.chrom.sizes
```
#---- now run this Juicer script in a core with many CPUs ----------#
_run_jucier.sh_
```
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
echo "Done Juicer Hi-C analysis.  Use YaHS to scaffold contigs further."
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
```

#---- next we run the YaHS to scaffold and fix misassemblies ------------#

#After running juicer, we use the file merged_nodups.txt to create a bed file for yahs (thanks to Eric for the script!)
```
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=0-3
#SBATCH --mem=128000

awk '
{
     if ($9 > 0 && $12 > 0) #this filters out mapping quality 0 reads
    {
    # Calculate end position for read 1 using pos1 and cigar1
    cigar = $10
    pos_start1 = $3
    sum1 = 0
    while (match(cigar, /([0-9]+)([MDNX=])/, arr)) {
        sum1 += arr[1]
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    pos_end1 = pos_start1 + sum1  # End position for read 1

    # Calculate end position for read 2 using pos2 and cigar2
    cigar = $13
    pos_start2 = $7
    sum2 = 0
    while (match(cigar, /([0-9]+)([MDNX=])/, arr)) {
        sum2 += arr[1]
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    pos_end2 = pos_start2 + sum2  # End position for read 2

    # Print interleaved output for read 1 and read 2
    print $2, pos_start1, pos_end1, $15"/1", $9
    print $6, pos_start2, pos_end2, $16"/2", $12
    }
}
' merged_nodups.txt > merged_nodups_for_yahs.bed
```

Now we run YaHS like this:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=257000
#SBATCH --output=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_hap1__hic.21Mar2024.out
#SBATCH --error=/home/kaedeh/scratch/Nettle/log_file/Nettle_female_juicer_hic.21Mar2024.err

#---------- This is a scaffolding program that takes Juicer merge_no_dup in bed format --------#

module load  StdEnv/2023  gcc/12.3  samtools/1.20
export PATH=$PATH:/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/yahs

#1. give indexed genome
input_refgenome=/home/kaedeh/scratch/Salmonberry/HiC_hap1/references/Salmonberry_redmorph_hifi_hifiasm_0.19.8_homcov_Q20_HiC.asm.hic.hap1.p_ctg.fa
samtools faidx $input_refgenome

#2. run yahs with the bed file indicating alignment of Hi-C reads to the refgenome.
yahs $input_refgenome merged_nodups_for_yahs.bed
```

Now we going to create our hic and assembly file to observe it in Juicebox
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=64
#SBATCH --mem=248G

module load StdEnv/2020 python/3.11.2 java/17.0.2 lastz/1.04.03

/home/egonza02/scratch/SOFTWARE/YAHS/yahs/juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp Harg2202_HIC_oldhifiasm.hic.hap1.fasta.fai >out_JBAT.log 2>&1

(java -jar -Xmx240G /home/egonza02/scratch/ASSEMBLIES/SCAFFOLDING/CPU_MODE/JUICER_NEW_FASTAS/ICGS/Harg2202_HIC_oldhifiasm.hic.hap1/aligned/Create_output_for_Yahs/RUN_JUICERAGAIN/scripts/common/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
```
#####################################

#to run yahs in salal (lab server), I run YaHS like this: 

#it is required that the /aligned/merged_nodups_for_yahs.bed is there, and also the juicer_tools.1.9.9_jcuda.0.8.jar is there in the directory you're calling this script. 
#_run_yahs.sh_
```
#!/bin/bash
export PATH=$PATH:/home/khirabayashi/bin/yahs
ref_asm=SAT.hifiasm.hic.hap1_hap2.p_ctg_inbred.fa

#1. scaffolding by yahs
yahs references/${ref_asm} \
aligned/merged_nodups_foryahs_onlyperfectmatches_softclip.bed

#2. convert yahs output to juicer_toos compatible format for manual checking. ### make sure the java juicer_tools is this v1.9.9!!! or else it doesn't work
juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp references/${ref_asm}.fai >out_JBAT.log 2>&1
mv out_JBAT.assembly ${ref_asm}.out_JBAT.assembly

#3. convert output to .hic so you can visualize in juicebox.
(java -jar -Xmx32G juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
mv out_JBAT.hic ${ref_asm}.out_JBAT.hic

mkdir yahs
mv yahs.out* out_JBAT* yahs

#4. Open Juicebox to correct contigs/scaffolds as necessary. 

#5. once you have the .review file again, run below to reproduce .fasta file. 
juicer post -o out_JBAT Salmonberry_redmorph_hifi_hifiasm_0.19.8_homcov_Q20_HiC.asm.hic.hap2.p_ctg_out_JBAT.review.1.assembly \
out_JBAT.liftover.agp ../references/*.fa     #the contig/previous assembly file you used to generate juicer data
```

#######################################################################################
#at this point change scaffold names to chromosome names that match the order and orientation of published reference genome.

### The step below further confirms SVs between haplotypes by combining H1 and H2 assembly, map Hi-C reads simultaneously, and visualize on Juicebox to manually check the Hi-C interactions. 
- You first cat concatenate haplotype assemblies
```
cat Round1_reviewed_H1.fasta Round1_reviewed_H2.fasta > Round1_reviewed_H1_H2_v1.fasta
```
- Run Juicer with Round1_reviewed_H1_H2_v1.fasta as the input reference genome.
_run_jucier.sh_ above. 

- Afterwards, I have tested the scaffolding with YaHS but this works poorly because YaHS doesn't recognize chromosome boundaries correctly and starts dividing them up and makes the assembly worse. 
So, we like to use 3ddna pipeline following the previous pipeline used in female nettle assembly. 

```
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

#run script; -r 0 runs only the first scaffolding, no polishing. -q 0 allows for mapping quality zero reads to be visible. 
run-asm-pipeline.sh -m haploid -r 0 -q 0 $contig_fasta $merged_nodups

deactivate


# ---------------------------------------------------------------------
echo "Done 3D-DNA pipeline scaffolding.  Use Juicebox to visualize and manually fix scaffolds."
# ---------------------------------------------------------------------

echo "Finished job at `date`"
```
Here, you load and check .hic and .assembly file with Juicebox. Manually correct assemblies H1 and H2 simultaneously if necesasry, but keep it minimal and only correct if they are obvious. 

#-------- Convert .assembly to .fasta with R script by Moshutava and Eric ---------#
```
#!/bin/bash
#SBATCH --account=
#SBATCH --time=1:00:00
#SBATCH --mem=30Gb
#SBATCH --output=
#SBATCH --error=

module load r

mkdir FASTA
Rscript asm_to_fasta_me.R Round_1_Nettle_male_hifiasm_yahs_hap1_hap2_scaff.0.review.4.assembly Round_1_Nettle_male_hifiasm_yahs_hap1_hap2_scaff.fa FASTA

```


#----- To further check if SVs between haplotypes are real ---------# 

#----- using precise coordinates (especially if the region doesn't fall between contigs and requires you to cut contigs for misassembly), run minimap2+SyRI ---------# 
```
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
```
Then from syri.out file, you can select for regions you want to check (eg. "INV") and copy those to Excel spreadsheet for record keeping. Go to Juicebox and check which orientations makes more sense.

