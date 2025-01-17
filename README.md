# Stinging nettle genome assembly
Codes/pipeline used in genome assembly report on _Urtica dioica_ ssp. _dioica_ (stinging nettle): https://www.mdpi.com/2223-7747/14/1/124

This assembly report consists of: 
- _De novo_ genome assembly with PacBio HiFi reads + Hi-C
- Gene annotation with BRAKER3 
- Repeat annotation with EDTA (TEs) and RepeatOBserver (telomeric sequences, centromeric repeats)
- Describe polycentric behaviour of nettle chromosomes
- Describe genes/proteins that are potentially responsible for stinging hairs
- Describe fibre synthesis genes presence
- Describe basic phylogenetics with related taxa in Urticaceae family. 

## _De novo_ assembly pipeline (exact steps I've taken)
1) Filter HiFi reads by Q20 [fastq.filter](https://github.com/LUMC/fastq-filter)
2) Run [Hififasm](https://hifiasm.readthedocs.io/en/latest/) with HiFi + Hi-C
3) Run [Juicer](https://github.com/aidenlab/juicer) to produce .hic and .assembly file - works as a matrix for Hi-C reads to then map onto. 
4) Use [3d-dna_pipeline](https://github.com/aidenlab/3d-dna) to map Hi-C reads to assembly
5) Open the output in Juicebox, manually fix contigs & assign chromosomes - see https://aidenlab.org/assembly/manual_180322.pdf for details
6) Repeat 3-5 on each haplotype separately - Round_1_H1.fa, Round_1_H2.fa
7) Combine H1 and H2 and repeat 3-5 - H1_H2_v1.fa
8) Repeat 3-5 on each haplotype separately again on the original contig-level assembly file - Round_2_H1.fa, Round_2_H2.fa
9) Combine H1 and H2 and repeat 3-5 - H1_H2_v2.fa
10) Perform a draft gene annotation with the Transdecoder pipeline used in lingonberry paper, ensuring that the annotated genes covered most of the chromosomes (BUSCO: 61% but genes were homogeneously distributed across the genome so I used it)
11) Error correct ONT reads using Canu, filter the output by Q10 and length 30kbp, then map in a further effort to close gaps in the assembly using TGS_gap_closer - did nothing. 
12) Use Anchorwave to align H1 and H2, verify any SVs that are present and whether they are real on Hi-C map
13) Use minimap2 to align H1 and H2 also, similarly do the same to check SVs --> looks like minimap2 identifies a lot more dubious organization of contigs in chromosomes
14) Use RepeatOBserver to check if some of the "difficult to resolve" regions are overlaping with centromeres or other repetitive regions - in that case, there is not much you can do to resolve the assembly... Also double check that homeologous chromosomes have reasonably similar centromeric positions. 
15) Perform final polishing with Juicebox on H1_H2_v2.fa file and ensure all SVs are real - Round_3_H1.fa, Round_3_H2.fa
16) Assign chromosome numbers that are homeologous between two haplotypes (13 chr each) - Nettle_female_H1_chr_renamed.fa, Nettle_female_H2_chr_renamed.fa
17) Based on SyRI output in step 13, manually check for INV/TRANS regions if they are real or not
18) First identify exact coordinates of such variations. Second change these in H1 / H2 assembly file in Juicebox. Third, if these are unclear, check in H1+H2 asesmbly file in Juicebox.
19) To check the uncertainties remain, I mapped HiFi reads and >30kb ONT reads to the assembly using Winnowmap, visualized the alignment in IGV - the point is to see if any of those long-reads span the breakpoints of such SVs to support the orientation/position of contigs.
20) Finalize the best guesses --> Round_5_H1.fa, Round_5_H2.fa
21) Align to the published _U. dioica_ genome from Darwin Tree of Life - if any variations are found, we make notes of it. I changed my chromosome numbering based on homologous relationship to this published genome (Udio_DToL). --> Nettle_female_Round_5_H1_chrname_renamed.fa, Nettle_female_Round_5_H2_chrname_renamed.fa
22) Due to some regions still showing messy Hi-C heatmap, I chopped off the obviously misplaced regions within a contig into debris, suspecting that some of these regions might be switch errors. --> Round_6_H1.fa, Round_6_H2.fa
23) Combined the H1+H2 --> H1_H2_v4.fa and run 3d-dna pipeline again, but with -q 0 setting to allow Hi-C reads to be mapped at repetitive regions. On the Hi-C heatmap with H1+H2 assembly, modify misassemblies as necessary. A lot of small contigs that were moved to debris in the previous steps could be placed back in chromosomes, in the correct haplotype (haplotype with more interaction than the other). --> Round_7_H1.fa, Round_7_H2.fa
24) Run minimap2+SyRI, check coordinates of INV and visually inspect in H1+H2 map. --> Round_8_H1.fa, Round_8_H2.fa (FINAL)
25) Run final BBmap(N50, etc.), BUSCO score on the assemblies.
26) Also check what the non-chromosomal regions contain using RED and BUSCO.


For list of intermediate files generated during the above process in a table, see my_juicebox.assembly_files.txt.
_also on the manuscript, I am referring to the different round number designations because I did not perform quality metric checks at every step on the way. I divided the Rounds of curation by distinctive steps with QC vakues available: Round 1 - visual Hi-C scaffolding with H1 and H2 and H1+H2 maps, Round 2 - attemp to incorporate ONT reads to fill the gaps, further curation with H1+H2 map and visually inspect, Round 3 - extensive manual curation using minimap2+SyRI coordinates to check INVs specificially, as well as checking if long-reads span breakpoints, Round 4 - final manual curation; using -q 0 on 3d-dna pipeline command, allowing mapping quality 0 reads to be visualized on the H1+H2 genome, corrected switch errors mainly to make the Hi-C map look smoother. Used minimap2+SyRI again to check exact INV coordinates only and made all change in H1+H2 map._


## _De novo_ annotation pipeline 
1) Softmask repetitive regions with [RepeatDetector](https://github.com/nextgenusfs/redmask)
2) Align published RNAseq to the softmasked genome using [Hisat2](https://daehwankimlab.github.io/hisat2/) and export as sorted BAM file
3) Annotate genes with [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) using softmasked genome, RNAseq BAM, and Viridiplantae (OrthoDB) protein database
4) Annotate TEs with [EDTA](https://github.com/oushujun/EDTA) on the pre-masked genonme

## Describe polycentric behaviour of chromosomes  
- We used [RepeatOBserver](https://github.com/celphin/RepeatOBserverV1/tree/main) to visualize general patterns of repeat structure and found some unique centromeric repeats distributions
- Search for specific centromeric/telomeric repeat sequences with RepeatOBserver results
- Compare syntenic relationship with Morus genome, that has previously shown polycentric behaviour in an experiment
- Compare synteny in Urticaceae family with [GENESPACE](https://github.com/jtlovell/GENESPACE)


## Describe genes/proteins that are potentially responsible for stinging hairs
- try search the protein sequence against full genome with [miniprot](https://github.com/lh3/miniprot)


