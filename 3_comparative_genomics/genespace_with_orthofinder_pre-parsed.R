###############################################
# -- change paths to those valid on your system
#genomeRepo <- "/home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/reference_data/genomeRepo"
wd <- "/home/kaedeh/scratch/Nettle/GENESPACE/Genespace_workingdir"
path2mcscanx <- "/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/MCScanX"
###############################################
library(devtools)
library(GENESPACE)

# -- before running GENESPACE parse the annotations to fastas with headers that match a gene bed file manually to avvoid all the annoying automated parsing
# -- load your genomes (make sure it matches your parsed bed and peptide.faa file names: genomeA.bed, genomeA.fa AND they have matching number of lines)
# if you want to only plot chromosomes, still use the full list of genes/annotations for this step and GENESPACE can handle that by subsetting the results after running the analysis during Riparian plot step.
genomes2run <- c("Boehmeria_nivea_wild", "Parietaria_judaica", "Urtica_urens", "Urtica_dioica_female_H1_v1")

# -- initalize the run and QC the inputs (this should be quick)
gpar <- init_genespace(
  wd = wd, nCores = 30,
  path2mcscanx = path2mcscanx)

# -- accomplish the run (this part tends to be memory-heavy so give it 90G)
out <- run_genespace(gpar)
save.image("genespace_out.RData")
