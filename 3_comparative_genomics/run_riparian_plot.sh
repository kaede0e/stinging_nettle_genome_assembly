#!/bin/bash
# Riparian plotting! 
# get an interactive job, load all modules and open R for riparian plotting: 
genomeRepo <- "/home/kaedeh/projects/def-gowens/kaedeh/Nettle/raw_data/reference_data/genomeRepo"
wd <- "/home/kaedeh/scratch/Nettle/GENESPACE/Genespace_workingdir_malevsfemale_blkRadius_5"
path2mcscanx <- "/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/MCScanX"
###############################################
library(devtools)
library(GENESPACE)

load("genespace_out_male_female_v3_blkRadius_5.RData")

### all 13 chromosomes, exclude small contigs.
pdf(file = 'riparian_plot_ordered_v3.pdf')
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("skyblue", "darkblue", "purple", "darkred", "salmon"))
ripDat <- plot_riparian(
  gsParam = out, 
  chrLabFun = function(x) gsub("^0", "", gsub("Urtica|dioica|chr|h1|h2|urtica|male|female|contig|_", "", tolower(x))),
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  refGenome = "Urtica_dioica_female_H1_v2", 
  useOrder = FALSE,
  useRegions = FALSE,
  backgroundColor = NULL, 
  chrExpand = 1,
  genomeIDs = c("Urtica_dioica_female_H1_v2", "Urtica_dioica_male_H1_v2", "Urtica_dioica_male_H2_v2", "Urtica_dioica_female_H2_v2"), 
  inversionColor = "#882255",
  forceRecalcBlocks = FALSE,
  scalePlotHeight = 1,
  scalePlotWidth = 5,
  minChrLen2plot = 100000)
dev.off()

### plotting only specified chromosome, with thicker chromosome
pdf(file = 'riparian_plot_ordered_v3_chr8.pdf')
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("skyblue", "darkblue", "purple", "darkred", "salmon"))
roi <- data.frame(
  genome = c("Urtica_dioica_female_H1_v2"), 
  chr = c("Urtica_dioica_female_chr_08"), 
  color = c("grey"))
ripDat <- plot_riparian(
  gsParam = out, 
  chrLabFun = function(x) gsub("^0", "", gsub("Urtica|dioica|chr|h1|h2|urtica|male|female|contig|_", "", tolower(x))),
  braidAlpha = .75,
  chrFill = "darkblue",
  addThemes = ggthemes,
  refGenome = "Urtica_dioica_female_H1_v2", 
  useOrder = FALSE,
  useRegions = FALSE,
  highlightBed = roi,
  backgroundColor = NULL, 
  chrExpand = 2,
  genomeIDs = c("Urtica_dioica_female_H1_v2", "Urtica_dioica_male_H1_v2", "Urtica_dioica_male_H2_v2", "Urtica_dioica_female_H2_v2"), 
  inversionColor = "#882255",
  forceRecalcBlocks = FALSE,
  minChrLen2plot = 100000,
  scalePlotHeight = 1,
  scalePlotWidth = 1)
dev.off()




pdf(file = 'riparian_plot.pdf')
ripDat <- plot_riparian(
  gsParam = out, 
  chrLabFun = function(x) gsub("^0", "", gsub("Salmonberry|redmorph|goldmorph|H1|H2|braker3|v2|anno|transcripts|^lg|_", "", tolower(x))),
  refGenome = "Salmonberry_redmorph_H1_braker3_v2.anno.transcripts", 
  genomeIDs = c("Salmonberry_redmorph_H1_braker3_v2.anno.transcripts", "Salmonberry_redmorph_H2_braker3_v2.anno.transcripts", "Salmonberry_goldmorph_H1_braker3_v2.anno.transcripts", "Salmonberry_goldmorph_H2_braker3_v2.anno.transcripts"), 
  forceRecalcBlocks = FALSE,
  minChrLen2plot = 100,
  scalePlotHeight = 1,
  scalePlotWidth = 2)
dev.off()

# ------ Rubus genome comparative analysis ------- #
# get an interactive job, load all modules and open R for riparian plotting: 
wd <- "/home/kaedeh/scratch/Salmonberry/red_vs_gold/genespace/Genespace_workingdir_Rubus"
path2mcscanx <- "/home/kaedeh/projects/def-gowens/kaedeh/cranberry_genome/bin/MCScanX"
###############################################
library(devtools)
library(GENESPACE)

load("genespace_out_Rubus_v1.RData")

### all 13 chromosomes, exclude small contigs.
pdf(file = 'riparian_plot_Rubus_v1_default.pdf')
#ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
#customPal <- colorRampPalette(
#  c("purple", "darkred", "salmon"))
ripDat <- plot_riparian(
  gsParam = out, 
  chrLabFun = function(x) gsub("^0", "", gsub("chr|h1|h2|scaff|scaffold|contig|_", "", tolower(x))),
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  refGenome = "Rubus_parviflorus_H1", 
  useOrder = FALSE,
  useRegions = FALSE,
  backgroundColor = NULL, 
  chrExpand = 1,
  genomeIDs = c("Rubus_parviflorus_H1", "Rubus_spectabilis_red_H1", "Rubus_Idaeus_JoanJ", "Argentina_anserina", "Fragaria_vesca", "Rosa_rugosa", "Malus_sylvestris", "Pyrus_bretschneideri", "Prunus_persica", "Dryas_octopetala"), 
#  inversionColor = "#882255",
  forceRecalcBlocks = FALSE,
  scalePlotHeight = 1,
  scalePlotWidth = 5,
  minChrLen2plot = 10000)
dev.off()


# ---------------------------------------------------------------------
