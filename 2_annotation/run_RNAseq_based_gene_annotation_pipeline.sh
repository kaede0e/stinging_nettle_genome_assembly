#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --account=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12000M
#SBATCH --output=
#SBATCH --error=

# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
echo "SLURM_JOBID: " $SLURM_JOBID
# ---------------------------------------------------------------------
echo ""

## RNAseq evidence-based gene annotation pipeline for novel genome assembly ##

# ---------------------------------------------------------------------

# Description of this pipeline:
#1. Map RNAseq/cDNA seq to draft genome assembly (HISAT2)
#2. Assemble transcript, retaining hits with BLAST alignment (StringTie)
#3. Find gene structure including gene, mRNA, CDS, UTR, etc. within transcripts (TransDecoder)
#4. Annotate genes with orthologs (EggNOG-mapper)

# ---------------------------------------------------------------------

#####################################
#### Execution of programmes ########
#####################################

export PATH=$PATH:/home/kaedeh/scratch/Nettle/annotation
export PATH=$PATH:/home/~/raw_data/reference_data/RNAseq #store my Illumina reads here
export PATH=$PATH:/home/~/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/
export PATH=$PATH:/home/~/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/ #this is not working so you need full path to invoke scripts
module load StdEnv/2020 gcc/9.3.0 blast+/2.13.0
module load hisat2 stringtie samtools bedtools

# ---------------------------------------------------------------------

#### 0. QC RNA reads & check all files are there ########
# Make sure to run fastqc/0.11.9 or fastp before using the data.
# Make sure to copy one final draft genome to working directory 
# Make sure to have all RNAseq .fastq files (paired end, stranded) in ~/raw_fastq_files/*_1.fastq or *_2.fastq, matching prefix for paired reads.

# ---------------------------------------------------------------------

input_refgenome=Round_2_hap1.reviewed.chr_assembled.fasta
assembly_name=Nettle_female_hap1

#### 1. Transcript alignment to genome ####
hisat2-build $input_refgenome ${assembly_name}_hisat2_index #took about 20min
ls /home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/reference_data/RNAseq/*_1_paired*.fastq.gz | tr '\n' ',' | sed 's/.$//' > RNAseq_m1_list_of_files.txt
ls /home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/reference_data/RNAseq/*_2_paired*.fastq.gz | tr '\n' ',' | sed 's/.$//' > RNAseq_m2_list_of_files.txt
hisat2 \
-q -x ${assembly_name}_hisat2_index -1 `cat RNAseq_m1_list_of_files.txt` -2 `cat RNAseq_m2_list_of_files.txt` -S ${assembly_name}_RNAseq_on_genome_hisat2_aln.sam
samtools sort -o ${assembly_name}_RNAseq_on_genome_hisat2_aln.bam ${assembly_name}_RNAseq_on_genome_hisat2_aln.sam

#echo "Finished Hisat2 alignment: `date`"

# ---------------------------------------------------------------------

#### 2. Transcript assembly and extract fasta sequence into fasta format ####
stringtie -o ${assembly_name}_stringtie_gene_structure_output.gtf ${assembly_name}_RNAseq_on_genome_hisat2_aln.bam #assembles the aligned transcripts into a file with structural definitions
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/EVidenceModeler-v2.0.0/EvmUtils/misc/cufflinks_gtf_to_alignment_gff3.pl \
${assembly_name}_stringtie_gene_structure_output.gtf > ${assembly_name}_stringtie_transcript_asm_output.gff #this changes structure file into gff

#echo "Finished StringTie transcript assembly: `date`"

# ---------------------------------------------------------------------

#### 3. Find CDS within transcripts ####
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl \
${assembly_name}_stringtie_gene_structure_output.gtf $input_refgenome > ${assembly_name}_stringtie_transcripts.fasta #extracts fasta seq of assembled transcript, in StringTie transcript coordinate (i.e. STR1.1)
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl \
${assembly_name}_stringtie_gene_structure_output.gtf > ${assembly_name}_stringtie_transcripts.gff3 #still in transcript coordinate
TransDecoder.LongOrfs -t ${assembly_name}_stringtie_transcripts.fasta

#echo "Finished TransDecoder long ORF prediction: `date`"

# BlastP library is prepared from UniProt Arabidopsis and Vaccinium known proteins#

blastp -query ${assembly_name}_stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep \
    -db /home/kaedeh/scratch/Nettle/blastp_lib/uniprot_Urticaceae_Arabidopsis_known_proteins.fasta \
    -outfmt 6 -evalue 1e-5 -num_threads 8 -out ${assembly_name}_blastp_Vaccinium_Arabidopsis_homology_based_genes.txt

#echo "Finished BlastP gene search for orthologous genes in Vaccinium/Arabidopsis: `date`"

TransDecoder.Predict -t ${assembly_name}_stringtie_transcripts.fasta --retain_blastp_hits ${assembly_name}_blastp_Vaccinium_Arabidopsis_homology_based_genes.txt
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl \
     ${assembly_name}_stringtie_transcripts.fasta.transdecoder.gff3 \
     ${assembly_name}_stringtie_transcripts.gff3 \
     ${assembly_name}_stringtie_transcripts.fasta > ${assembly_name}_stringtie_transcripts.fasta.transdecoder.genome.gff3 #save this final .gff3 file and make sure it is accessible to you.

awk -F ";Name" '{print $1}' ${assembly_name}_stringtie_transcripts.fasta.transdecoder.genome.gff3 > ${assembly_name}_stringtie_transcripts.fasta.transdecoder.genome.clean.gff3 #cleans up the annotation file with long headline.

echo "Finished TransDecoder CDS prediction: `date`"

mkdir final_annotation_files
scp $input_refgenome final_annotation_files/ #genome assembly
scp ${assembly_name}_stringtie_transcripts.fasta.transdecoder.genome.clean.gff3 final_annotation_files/${assembly_name}_genes_anno.gff3 #gene positions and features annotated on the above genome

echo "Finished TranDecoder genome annotation with genes: `date`"

# ---------------------------------------------------------------------

#### 4. Gene functional annotation on EggNOG-mapper web ####
# get a protein.fasta from annotation (.gff3 mapped to genome)
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gff3_file_to_bed.pl ${assembly_name}_stringtie_transcripts.fasta.transdecoder.genome.clean.gff3 > ${assembly_name}_stringtie_transcripts.fasta.transdecoder.genome.clean.bed
perl /home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/TransDecoder/TransDecoder-TransDecoder-v5.5.0/util/gff3_file_to_proteins.pl \
--gff3 final_annotation_files/${assembly_name}_genes_anno.gff3 \
--fasta $input_refgenome > final_annotation_files/${assembly_name}_genes_anno.faa
# take the longest gene from isoforms/splicing variants
cat final_annotation_files/${assembly_name}_anno.faa | perl /home/kaedeh/projects/rrg-gowens/kaedeh/Lingonberry/scripts/fasta2longestisoform.pl > final_annotation_files/${assembly_name}_longest_genes_anno.faa

echo "Finished converting nt .fasta into proteins .fasta: `date`"
echo "Congrats, your final_annotation_files/*anno.faa is ready for submission to EggNOG-mapper!: `date`"

# ---------------------------------------------------------------------
