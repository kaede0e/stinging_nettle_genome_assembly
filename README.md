# Stinging nettle genome assembly
Codes/pipeline used in genome assembly report on _Urtica dioica_ssp. _dioica_ (stinging nettle)

This assembly report consists of: 
- _De novo_ genome assembly with PacBio HiFi reads + Hi-C
- Gene annotation with BRAKER3 (planned)
- Repeat annotation with EDTA  (planned) and RepeatOBserver
- Describe polycentric behaviour of nettle chromosomes
- Describe genes/proteins that are potentially responsible for stinging hairs
- Describe fibre synthesis genes presence
- Describe basic phylogenetics with related taxa in Urticaceae family. 

_De novo_ assembly pipeline (exact steps I've taken)
1) Filter HiFi reads by Q20
2) Run Hififasm with HiFi + Hi-C
3) Run Juicer to produce .hic and .assembly file - works as a matrix for Hi-C reads to then map onto. 
4) Use 3d-dna_pipeline to map Hi-C reads to assembly
5) Open the output in Juicebox, manually fix contigs & assign chromosomes
6) Repeat 3-5 on each haplotype separately - Round_1_H1.fa, Round_1_H2.fa
7) Combine H1 and H2 and repeat 3-5 - H1_H2_v1.fa
8) Repeat 3-5 on each haplotype separately again on the original contig-level assembly file - Round_2_H1.fa, Round_2_H2.fa
9) Combine H1 and H2 and repeat 3-5 - H1_H2_v2.fa
10) Perform a draft gene annotation with the Transdecoder pipeline used in lingonberry paper, ensuring that the annotated genes covered most of the chromosomes (BUSCO: 61% but genes were homogeneously distributed across the genome so I used it)
11) Use Anchorwave to align H1 and H2, verify any SVs that are present and whether they are real on Hi-C map
12) Use minimap2 to align H1 and H2 also, similarly do the same to check SVs
13) Use RepeatOBserver to check if some of the "difficult to resolve" regions are overlaping with centromeres or other repetitive regions - in that case, there is not much you can do to resolve the assembly... Also double check that homeologous chromosomes have reasonably similar centromeric positions. 
14) Perform final polishing with Juicebox on H1_H2_v2.fa file and ensure all SVs are real - Round_3_H1.fa, Round_3_H2.fa
15) Assign chromosome numbers that are homeologous between two haplotypes (13 chr each) - Nettle_female_H1_chr_renamed.fa, Nettle_female_H2_chr_renamed.fa
16) Run final Bandage (N50, etc.), BUSCO score, and Merqury (QV, k-mer completeness, etc.) on the assemblies.


