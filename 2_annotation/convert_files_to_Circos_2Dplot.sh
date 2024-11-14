#bin/bash 
#convert files to Circos format (see manuscript Figure 1 for which tracks)
#some codes are written in R. 

#a. SyRI output tables to synteny ribbons
cat *synOut.txt | grep "#" | awk '{print $2,$3,$4,$6,$7,$8}' | awk '{print $1,$2,$3,$4"_H2",$5,$6}' | sed s/' '/'    '/g > synOut_Circos.tsv
cat *invOut.txt | grep "#" | awk '{print $2,$3,$4,$6,$7,$8}' | awk '{print $1,$2,$3,$4"_H2",$5,$6}' | sed s/' '/'    '/g > invOut_Circos.tsv
cat *TLOut.txt | grep "#" | awk '{print $2,$3,$4,$6,$7,$8}' | awk '{print $1,$2,$3,$4"_H2",$5,$6}' | sed s/' '/'    '/g > TLOut_Circos.tsv
cat *invTLOut.txt | grep "#" | awk '{print $2,$3,$4,$6,$7,$8}'| awk '{print $1,$2,$3,$4"_H2",$5,$6}' | sed s/' '/'    '/g > invTLOut_Circos.tsv
cat *dupOut.txt | grep "#" | awk '{print $2,$3,$4,$6,$7,$8}' | awk '{print $1,$2,$3,$4"_H2",$5,$6}' | sed s/' '/'    '/g > dupOut_Circos.tsv
cat *invDupOut.txt | grep "#" | awk '{print $2,$3,$4,$6,$7,$8}' | awk '{print $1,$2,$3,$4"_H2",$5,$6}' | sed s/' '/'    '/g > invDupOut_Circos.tsv
### The above script prints _H2 on the end of each haplotype 2 chromosomes since Circos wouldn't allow duplicate chromosome names in their files.  

### Similarly, for the following three tracks, you have to change the first column of the haplotyep 2 "Urtica_dioica_female_chr_01" to "Urtica_dioica_female_chr_01_H2" to be compatible with my Circos plot. 
cat convert_files_to_2Dplot_Circos_histogram_H2.R | \
sed s/"Urtica_dioica_female_chr_01"/"Urtica_dioica_female_chr_01_H2"/g | \
sed s/"Urtica_dioica_female_chr_02"/"Urtica_dioica_female_chr_02_H2"/g | \
sed s/"Urtica_dioica_female_chr_03"/"Urtica_dioica_female_chr_03_H2"/g | \
sed s/"Urtica_dioica_female_chr_04"/"Urtica_dioica_female_chr_04_H2"/g | \
sed s/"Urtica_dioica_female_chr_05"/"Urtica_dioica_female_chr_05_H2"/g | \
sed s/"Urtica_dioica_female_chr_06"/"Urtica_dioica_female_chr_06_H2"/g | \
sed s/"Urtica_dioica_female_chr_07"/"Urtica_dioica_female_chr_07_H2"/g | \
sed s/"Urtica_dioica_female_chr_08"/"Urtica_dioica_female_chr_08_H2"/g | \
sed s/"Urtica_dioica_female_chr_09"/"Urtica_dioica_female_chr_09_H2"/g | \
sed s/"Urtica_dioica_female_chr_10"/"Urtica_dioica_female_chr_10_H2"/g | \
sed s/"Urtica_dioica_female_chr_11"/"Urtica_dioica_female_chr_11_H2"/g | \
sed s/"Urtica_dioica_female_chr_12"/"Urtica_dioica_female_chr_12_H2"/g | \
sed s/"Urtica_dioica_female_chr_13"/"Urtica_dioica_female_chr_13_H2"/g

cat convert_files_to_2Dplot_Circos_histogram_gene.H2.R | \
sed s/"Urtica_dioica_female_chr_01_H2"/"Urtica_dioica_female_chr_01"/g | \
sed s/"Urtica_dioica_female_chr_02_H2"/"Urtica_dioica_female_chr_02"/g | \
sed s/"Urtica_dioica_female_chr_03_H2"/"Urtica_dioica_female_chr_03"/g | \
sed s/"Urtica_dioica_female_chr_04_H2"/"Urtica_dioica_female_chr_04"/g | \
sed s/"Urtica_dioica_female_chr_05_H2"/"Urtica_dioica_female_chr_05"/g | \
sed s/"Urtica_dioica_female_chr_06_H2"/"Urtica_dioica_female_chr_06"/g | \
sed s/"Urtica_dioica_female_chr_07_H2"/"Urtica_dioica_female_chr_07"/g | \
sed s/"Urtica_dioica_female_chr_08_H2"/"Urtica_dioica_female_chr_08"/g | \
sed s/"Urtica_dioica_female_chr_09_H2"/"Urtica_dioica_female_chr_09"/g | \
sed s/"Urtica_dioica_female_chr_10_H2"/"Urtica_dioica_female_chr_10"/g | \
sed s/"Urtica_dioica_female_chr_11_H2"/"Urtica_dioica_female_chr_11"/g | \
sed s/"Urtica_dioica_female_chr_12_H2"/"Urtica_dioica_female_chr_12"/g | \
sed s/"Urtica_dioica_female_chr_13_H2"/"Urtica_dioica_female_chr_13"/g


for files in `ls *.tsv`;
do
	tail -n +2 $files > header_removed_${files}
done

cat header_removed_Gene_per_500kb_window_H2* | \
sed s/"Urtica_dioica_female_chr_01"/"Urtica_dioica_female_chr_01_H2"/g | \
sed s/"Urtica_dioica_female_chr_02"/"Urtica_dioica_female_chr_02_H2"/g | \
sed s/"Urtica_dioica_female_chr_03"/"Urtica_dioica_female_chr_03_H2"/g | \
sed s/"Urtica_dioica_female_chr_04"/"Urtica_dioica_female_chr_04_H2"/g | \
sed s/"Urtica_dioica_female_chr_05"/"Urtica_dioica_female_chr_05_H2"/g | \
sed s/"Urtica_dioica_female_chr_06"/"Urtica_dioica_female_chr_06_H2"/g | \
sed s/"Urtica_dioica_female_chr_07"/"Urtica_dioica_female_chr_07_H2"/g | \
sed s/"Urtica_dioica_female_chr_08"/"Urtica_dioica_female_chr_08_H2"/g | \
sed s/"Urtica_dioica_female_chr_09"/"Urtica_dioica_female_chr_09_H2"/g | \
sed s/"Urtica_dioica_female_chr_10"/"Urtica_dioica_female_chr_10_H2"/g | \
sed s/"Urtica_dioica_female_chr_11"/"Urtica_dioica_female_chr_11_H2"/g | \
sed s/"Urtica_dioica_female_chr_12"/"Urtica_dioica_female_chr_12_H2"/g | \
sed s/"Urtica_dioica_female_chr_13"/"Urtica_dioica_female_chr_13_H2"/g > Gene_per_500kb_window_H2_allChr.tsv

cat header_removed_Gene_per_500kb_window_H1* > Gene_per_500kb_window_H1_allChr.tsv

### The three tracks above is referring to are create by Rscripts as shown below (it's very lengthy and repetitive...) Probably better to code in a different language that I don't know :(

#b. gene annotation (gff3 file from BRAKER3 output) to gene density in histogram, per 500kbp windows
Rscript gene_anno.R 

#---------gene_anno.R --------------#
library(tidyverse)
library(dplyr)
library(ggplot2)

##### 3. BRAKER3 gene density as histogram plot in Circos #####

# Makes sure all windows are in the dataframe, with the 0 hits # 
window_size = 500000 #500kb window is probably best for resolution in Circos for both haplotypes

# H1 ----------------------------------------
empty_windows <- read.table("Nettle_female_H1_karyotype.txt") %>%
  as_tibble(.)%>% 
  rename(CHROM = V1, length = V2) %>% 
  mutate(max_window_size = floor(length/window_size))

## Do this per chromosome and adjust the window number according to "empty_windows" ##
# refer to the dataframe_prep_for_het_plots.R
Gene_H1_Chr13 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_13")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:42) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_13")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr13_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr13 <- Gene_H1_Chr13 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr13_allwindows <- left_join(Chr13_500kb_windows, Gene_per_500kb_window_H1_Chr13) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr13_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr13.tsv")

Gene_H1_Chr12 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_12")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:49) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_12")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr12_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr12 <- Gene_H1_Chr12 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr12_allwindows <- left_join(Chr12_500kb_windows, Gene_per_500kb_window_H1_Chr12) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr12_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr12.tsv")

Gene_H1_Chr11 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_11")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:49) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_11")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr11_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr11 <- Gene_H1_Chr11 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr11_allwindows <- left_join(Chr11_500kb_windows, Gene_per_500kb_window_H1_Chr11) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr11_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr11.tsv")

Gene_H1_Chr10 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_10")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:56) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_10")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr10_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr10 <- Gene_H1_Chr10 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr10_allwindows <- left_join(Chr10_500kb_windows, Gene_per_500kb_window_H1_Chr10) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr10_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr10.tsv")

Gene_H1_Chr9 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_09")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:60) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_09")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr9_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr9 <- Gene_H1_Chr9 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr9_allwindows <- left_join(Chr9_500kb_windows, Gene_per_500kb_window_H1_Chr9) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr9_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr9.tsv")

Gene_H1_Chr8 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_08")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:86) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_08")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr8_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr8 <- Gene_H1_Chr8 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr8_allwindows <- left_join(Chr8_500kb_windows, Gene_per_500kb_window_H1_Chr8) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr8_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr8.tsv")

Gene_H1_Chr7 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_07")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:96) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_07")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr7_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr7 <- Gene_H1_Chr7 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr7_allwindows <- left_join(Chr7_500kb_windows, Gene_per_500kb_window_H1_Chr7) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr7_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr7.tsv")

Gene_H1_Chr6 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_06")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:90) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_06")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr6_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr6 <- Gene_H1_Chr6 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr6_allwindows <- left_join(Chr6_500kb_windows, Gene_per_500kb_window_H1_Chr6) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr6_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr6.tsv")

Gene_H1_Chr5 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_05")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:86) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_05")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr5_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr5 <- Gene_H1_Chr5 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr5_allwindows <- left_join(Chr5_500kb_windows, Gene_per_500kb_window_H1_Chr5) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr5_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr5.tsv")

Gene_H1_Chr4 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_04")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:98) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_04")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr4_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr4 <- Gene_H1_Chr4 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr4_allwindows <- left_join(Chr4_500kb_windows, Gene_per_500kb_window_H1_Chr4) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr4_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr4.tsv")

Gene_H1_Chr3 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_03")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:103) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_03")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr3_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr3 <- Gene_H1_Chr3 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr3_allwindows <- left_join(Chr3_500kb_windows, Gene_per_500kb_window_H1_Chr3) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr3_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr3.tsv")

Gene_H1_Chr2 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_02")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:110) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_02")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr2_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr2 <- Gene_H1_Chr2 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr2_allwindows <- left_join(Chr2_500kb_windows, Gene_per_500kb_window_H1_Chr2) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr2_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr2.tsv")

Gene_H1_Chr1 <- read.table("braker3.Nettle_female_H1_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_01")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:118) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_01")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr1_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H1_Chr1 <- Gene_H1_Chr1 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr1_allwindows <- left_join(Chr1_500kb_windows, Gene_per_500kb_window_H1_Chr1) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr1_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr1.tsv")

# H2 ----------------------------------------
#empty_windows <- read.table("Nettle_female_H2_karyotype.txt") %>%
#  as_tibble(.)%>%
#  rename(CHROM = V1, length = V2) %>%
#  mutate(max_window_size = floor(length/window_size)) # Determine max number of windows per chromosome and change below.

## Do this per chromosome and adjust the window number according to "empty_windows" ##
# refer to the dataframe_prep_for_het_plots.R
Gene_H2_Chr13 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_13")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:43) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_13")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr13_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr13 <- Gene_H2_Chr13 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr13_allwindows <- left_join(Chr13_500kb_windows, Gene_per_500kb_window_H2_Chr13) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr13_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr13.tsv")

Gene_H2_Chr12 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_12")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:48) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_12")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr12_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr12 <- Gene_H2_Chr12 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr12_allwindows <- left_join(Chr12_500kb_windows, Gene_per_500kb_window_H2_Chr12) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr12_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr12.tsv")

Gene_H2_Chr11 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_11")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:49) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_11")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr11_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr11 <- Gene_H2_Chr11 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr11_allwindows <- left_join(Chr11_500kb_windows, Gene_per_500kb_window_H2_Chr11) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr11_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr11.tsv")

Gene_H2_Chr10 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_10")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:54) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_10")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr10_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr10 <- Gene_H2_Chr10 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr10_allwindows <- left_join(Chr10_500kb_windows, Gene_per_500kb_window_H2_Chr10) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr10_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr10.tsv")

Gene_H2_Chr9 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_09")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:60) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_09")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr9_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr9 <- Gene_H2_Chr9 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr9_allwindows <- left_join(Chr9_500kb_windows, Gene_per_500kb_window_H2_Chr9) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr9_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr9.tsv")

Gene_H2_Chr8 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_08")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:72) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_08")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr8_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr8 <- Gene_H2_Chr8 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr8_allwindows <- left_join(Chr8_500kb_windows, Gene_per_500kb_window_H2_Chr8) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr8_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr8.tsv")

Gene_H2_Chr7 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_07")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:98) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_07")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr7_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr7 <- Gene_H2_Chr7 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr7_allwindows <- left_join(Chr7_500kb_windows, Gene_per_500kb_window_H2_Chr7) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr7_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr7.tsv")

Gene_H2_Chr6 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_06")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:92) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_06")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr6_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr6 <- Gene_H2_Chr6 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr6_allwindows <- left_join(Chr6_500kb_windows, Gene_per_500kb_window_H2_Chr6) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr6_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr6.tsv")

Gene_H2_Chr5 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_05")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:85) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_05")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr5_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr5 <- Gene_H2_Chr5 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr5_allwindows <- left_join(Chr5_500kb_windows, Gene_per_500kb_window_H2_Chr5) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr5_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr5.tsv")

Gene_H2_Chr4 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_04")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:97) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_04")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr4_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr4 <- Gene_H2_Chr4 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr4_allwindows <- left_join(Chr4_500kb_windows, Gene_per_500kb_window_H2_Chr4) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr4_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr4.tsv")

Gene_H2_Chr3 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_03")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:100) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_03")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr3_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr3 <- Gene_H2_Chr3 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr3_allwindows <- left_join(Chr3_500kb_windows, Gene_per_500kb_window_H2_Chr3) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr3_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr3.tsv")

Gene_H2_Chr2 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_02")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:113) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_02")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr2_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr2 <- Gene_H2_Chr2 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr2_allwindows <- left_join(Chr2_500kb_windows, Gene_per_500kb_window_H2_Chr2) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr2_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr2.tsv")

Gene_H2_Chr1 <- read.table("braker3.Nettle_female_H2_v1_RNA_ProtsViridiplantaeOrthoDB_gene.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_01")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:119) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_01")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr1_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
Gene_per_500kb_window_H2_Chr1 <- Gene_H2_Chr1 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr1_allwindows <- left_join(Chr1_500kb_windows, Gene_per_500kb_window_H2_Chr1) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr1_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr1.tsv")
#--------------#


#c. TE annotation (bed formatted EDTA output) to TE density in histogram, per 500kbp windows 
Rscript TE_anno.R 

#---------TE_anno.R --------------#
library(tidyverse)
library(dplyr)
library(ggplot2)

##### 2. EDTA TE density as histogram plot in Circos #####

# Makes sure all windows are in the dataframe, with the 0 hits # 
window_size = 500000 #500kb window is probably best for resolution in Circos for both haplotypes

# H1 ----------------------------------------
empty_windows <- read.table("Nettle_female_H1_karyotype.txt") %>%
  as_tibble(.)%>% 
  rename(CHROM = V1, length = V2) %>% 
  mutate(max_window_size = floor(length/window_size))

## Do this per chromosome and adjust the window number according to "empty_windows" ##
# refer to the dataframe_prep_for_het_plots.R
TE_H1_Chr13 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1, 
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "13") %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_13")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:42) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_13")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr13_500kb_windows <- data_frame %>% 
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr13 <- TE_H1_Chr13 %>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr13_allwindows <- left_join(Chr13_500kb_windows, TE_per_500kb_window_H1_Chr13) %>% 
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr13_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr13.tsv")

TE_H1_Chr12 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "12") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_12")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:49) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_12")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr12_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr12 <- TE_H1_Chr12 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr12_allwindows <- left_join(Chr12_500kb_windows, TE_per_500kb_window_H1_Chr12) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr12_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr12.tsv")

TE_H1_Chr11 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "11") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_11")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:49) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_11")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr11_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr11 <- TE_H1_Chr11 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr11_allwindows <- left_join(Chr11_500kb_windows, TE_per_500kb_window_H1_Chr11) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr11_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr11.tsv")

TE_H1_Chr10 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "10") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_10")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:56) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_10")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr10_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr10 <- TE_H1_Chr10 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr10_allwindows <- left_join(Chr10_500kb_windows, TE_per_500kb_window_H1_Chr10) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr10_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr10.tsv")

TE_H1_Chr9 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "09") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_09")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:60) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_09")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr9_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr9 <- TE_H1_Chr9 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr9_allwindows <- left_join(Chr9_500kb_windows, TE_per_500kb_window_H1_Chr9) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr9_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr9.tsv")

TE_H1_Chr8 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "08") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_08")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:86) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_08")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr8_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr8 <- TE_H1_Chr8 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr8_allwindows <- left_join(Chr8_500kb_windows, TE_per_500kb_window_H1_Chr8) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr8_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr8.tsv")

TE_H1_Chr7 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "07") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_07")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:96) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_07")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr7_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr7 <- TE_H1_Chr7 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr7_allwindows <- left_join(Chr7_500kb_windows, TE_per_500kb_window_H1_Chr7) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr7_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr7.tsv")

TE_H1_Chr6 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "06") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_06")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:90) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_06")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr6_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr6 <- TE_H1_Chr6 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr6_allwindows <- left_join(Chr6_500kb_windows, TE_per_500kb_window_H1_Chr6) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr6_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr6.tsv")

TE_H1_Chr5 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "05") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_05")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:86) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_05")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr5_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr5 <- TE_H1_Chr5 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr5_allwindows <- left_join(Chr5_500kb_windows, TE_per_500kb_window_H1_Chr5) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr5_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr5.tsv")

TE_H1_Chr4 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "04") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_04")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:98) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_04")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr4_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr4 <- TE_H1_Chr4 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr4_allwindows <- left_join(Chr4_500kb_windows, TE_per_500kb_window_H1_Chr4) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr4_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr4.tsv")

TE_H1_Chr3 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "03") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_03")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:103) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_03")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr3_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr3 <- TE_H1_Chr3 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr3_allwindows <- left_join(Chr3_500kb_windows, TE_per_500kb_window_H1_Chr3) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr3_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr3.tsv")

TE_H1_Chr2 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "02") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_02")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:110) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_02")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr2_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr2 <- TE_H1_Chr2 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr2_allwindows <- left_join(Chr2_500kb_windows, TE_per_500kb_window_H1_Chr2) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr2_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr2.tsv")

TE_H1_Chr1 <- read.table("Nettle_female_H1.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "01") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_01")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:118) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_01")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr1_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H1_Chr1 <- TE_H1_Chr1 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr1_allwindows <- left_join(Chr1_500kb_windows, TE_per_500kb_window_H1_Chr1) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr1_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr1.tsv")

# H2 ----------------------------------------
empty_windows <- read.table("Nettle_female_H2_karyotype.txt") %>%
  as_tibble(.)%>% 
  rename(CHROM = V1, length = V2) %>% 
  mutate(max_window_size = floor(length/window_size))

## Do this per chromosome and adjust the window number according to "empty_windows" ##
# refer to the dataframe_prep_for_het_plots.R
TE_H2_Chr13 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "13") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_13_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:43) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_13_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr13_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr13 <- TE_H2_Chr13 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr13_allwindows <- left_join(Chr13_500kb_windows, TE_per_500kb_window_H2_Chr13) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr13_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr13.tsv")

TE_H2_Chr12 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "12") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_12_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:48) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_12_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr12_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr12 <- TE_H2_Chr12 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr12_allwindows <- left_join(Chr12_500kb_windows, TE_per_500kb_window_H2_Chr12) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr12_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr12.tsv")

TE_H2_Chr11 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "11") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_11_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:49) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_11_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr11_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr11 <- TE_H2_Chr11 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr11_allwindows <- left_join(Chr11_500kb_windows, TE_per_500kb_window_H2_Chr11) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr11_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr11.tsv")

TE_H2_Chr10 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "10") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_10_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:54) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_10_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr10_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr10 <- TE_H2_Chr10 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr10_allwindows <- left_join(Chr10_500kb_windows, TE_per_500kb_window_H2_Chr10) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr10_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr10.tsv")

TE_H2_Chr9 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "9") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_09_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:60) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_09_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr9_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr9 <- TE_H2_Chr9 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr9_allwindows <- left_join(Chr9_500kb_windows, TE_per_500kb_window_H2_Chr9) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr9_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr9.tsv")

TE_H2_Chr8 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "8") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_08_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:72) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_08_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr8_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr8 <- TE_H2_Chr8 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr8_allwindows <- left_join(Chr8_500kb_windows, TE_per_500kb_window_H2_Chr8) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr8_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr8.tsv")

TE_H2_Chr7 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "7") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_07_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:98) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_07_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr7_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr7 <- TE_H2_Chr7 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr7_allwindows <- left_join(Chr7_500kb_windows, TE_per_500kb_window_H2_Chr7) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr7_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr7.tsv")

TE_H2_Chr6 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "6") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_06_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:92) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_06_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr6_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr6 <- TE_H2_Chr6 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr6_allwindows <- left_join(Chr6_500kb_windows, TE_per_500kb_window_H2_Chr6) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr6_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr6.tsv")

TE_H2_Chr5 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "5") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_05_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:85) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_05_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr5_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr5 <- TE_H2_Chr5 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr5_allwindows <- left_join(Chr5_500kb_windows, TE_per_500kb_window_H2_Chr5) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr5_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr5.tsv")

TE_H2_Chr4 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "4") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_04_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:97) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_04_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr4_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr4 <- TE_H2_Chr4 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr4_allwindows <- left_join(Chr4_500kb_windows, TE_per_500kb_window_H2_Chr4) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr4_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr4.tsv")

TE_H2_Chr3 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "3") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_03_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:100) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_03_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr3_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr3 <- TE_H2_Chr3 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr3_allwindows <- left_join(Chr3_500kb_windows, TE_per_500kb_window_H2_Chr3) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr3_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr3.tsv")

TE_H2_Chr2 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "2") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_02_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:113) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_02_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr2_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr2 <- TE_H2_Chr2 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr2_allwindows <- left_join(Chr2_500kb_windows, TE_per_500kb_window_H2_Chr2) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr2_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr2.tsv")

TE_H2_Chr1 <- read.table("Nettle_female_H2.chr.fa.mod.EDTA.TEanno.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHR = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHR == "1") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_01_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:119) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_01_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Chr1_500kb_windows <- data_frame %>%
  mutate(window = window_1 - window_size,
         window_end = window + (window_size - 1)) %>%
  select(CHROM, window, window_end)
TE_per_500kb_window_H2_Chr1 <- TE_H2_Chr1 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H2_Chr1_allwindows <- left_join(Chr1_500kb_windows, TE_per_500kb_window_H2_Chr1) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H2_Chr1_allwindows, "Circos_input/TE_per_500kb_window_H2_Chr1.tsv")
#--------------#

#d. Shannon diversity score calculated in each 5kbp windows across the genome to Mean shannon diversity score over 250 kbp window in line plot
Rscript Shannon_line.R 

#---------Shannon_line.R --------------#
library(tidyverse)
library(dplyr)
library(ggplot2)

##### 1. RepeatOBserver Shannon diversity plot for line plot in Circos #####

# Makes sure all windows are in the dataframe, with the 0 hits # 
window_size = 250000 #250kb window is probably best for resolution in Circos for both haplotypes

# H1 ----------------------------------------
empty_windows <- read.table("Nettle_female_H1_karyotype.txt") %>%
  as_tibble(.)%>% 
  rename(CHROM = V1, length = V2) %>% 
  mutate(max_window_size = floor(length/window_size))

## Do this per chromosome and adjust the window number according to "empty_windows" ##
# refer to the dataframe_prep_for_het_plots.R
RO_H1_Chr13 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr13_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_13")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:84) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_13")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr13 <- RO_H1_Chr13 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_13")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr13, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr13.tsv")

RO_H1_Chr12 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr12_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_12")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:97) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_12")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr12 <- RO_H1_Chr12 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_12")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr12, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr12.tsv")

RO_H1_Chr11 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr11_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_11")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:97) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_11")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr11 <- RO_H1_Chr11 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_11")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr11, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr11.tsv")

RO_H1_Chr10 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr10_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_10")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:111) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_10")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr10 <- RO_H1_Chr10 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_10")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr10, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr10.tsv")

RO_H1_Chr9 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr9_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_09")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:120) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_09")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr9 <- RO_H1_Chr9 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_09")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr9, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr9.tsv")

RO_H1_Chr8 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr8_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_08")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:172) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_08")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr8 <- RO_H1_Chr8 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_08")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr8, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr8.tsv")

RO_H1_Chr7 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr7_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_07")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:191) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_07")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr7 <- RO_H1_Chr7 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_07")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr7, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr7.tsv")

RO_H1_Chr6 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr6_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_06")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:179) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_06")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr6 <- RO_H1_Chr6 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_06")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr6, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr6.tsv")

RO_H1_Chr5 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr5_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_05")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:172) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_05")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr5 <- RO_H1_Chr5 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_05")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr5, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr5.tsv")

RO_H1_Chr4 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr4_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_04")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:196) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_04")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr4 <- RO_H1_Chr4 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_04")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr4, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr4.tsv")

RO_H1_Chr3 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr3_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_03")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:206) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_03")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr3 <- RO_H1_Chr3 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_03")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr3, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr3.tsv")

RO_H1_Chr2 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr2_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_02")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:219) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_02")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr2 <- RO_H1_Chr2 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_02")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr2, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr2.tsv")

RO_H1_Chr1 <- read.table("RepeatOBserver_results/NettleFemaleR5_H1-AT_Chr1_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_01")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:236) {   #the max. number = empty_windows$max_window_size+1                   
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_01")    
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H1_Chr1 <- RO_H1_Chr1 %>% 
  filter(CHROM == "Urtica_dioica_female_chr_01")%>% 
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H1_Chr1, "Circos_input/Shannon_mean_per_250kb_window_H1_Chr1.tsv")

# H2 ----------------------------------------
empty_window_H2 <- read.table("Nettle_female_H2_karyotype.txt") %>%
  as_tibble(.)%>% 
  rename(CHROM = V1, length = V2) %>% 
  mutate(max_window_size = floor(length/window_size))

RO_H2_Chr13 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr13_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_13_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:85) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_13_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr13 <- RO_H2_Chr13 %>%
  filter(CHROM == "Urtica_dioica_female_chr_13_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr13, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr13.tsv")

RO_H2_Chr12 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr12_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_12_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:95) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_12_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr12 <- RO_H2_Chr12 %>%
  filter(CHROM == "Urtica_dioica_female_chr_12_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr12, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr12.tsv")

RO_H2_Chr11 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr11_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_11_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:97) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_11_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr11 <- RO_H2_Chr11 %>%
  filter(CHROM == "Urtica_dioica_female_chr_11_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr11, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr11.tsv")

RO_H2_Chr10 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr10_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_10_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:107) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_10_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr10 <- RO_H2_Chr10 %>%
  filter(CHROM == "Urtica_dioica_female_chr_10_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr10, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr10.tsv")

RO_H2_Chr9 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr9_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_09_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:119) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_09_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr9 <- RO_H2_Chr9 %>%
  filter(CHROM == "Urtica_dioica_female_chr_09_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr9, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr9.tsv")

RO_H2_Chr8 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr8_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_08_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:143) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_08_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr8 <- RO_H2_Chr8 %>%
  filter(CHROM == "Urtica_dioica_female_chr_08_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr8, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr8.tsv")

RO_H2_Chr7 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr7_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_07_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:195) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_07_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr7 <- RO_H2_Chr7 %>%
  filter(CHROM == "Urtica_dioica_female_chr_07_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr7, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr7.tsv")

RO_H2_Chr6 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr6_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_06_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:183) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_06_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr6 <- RO_H2_Chr6 %>%
  filter(CHROM == "Urtica_dioica_female_chr_06_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr6, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr6.tsv")

RO_H2_Chr5 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr5_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_05_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:169) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_05_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr5 <- RO_H2_Chr5 %>%
  filter(CHROM == "Urtica_dioica_female_chr_05_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr5, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr5.tsv")

RO_H2_Chr4 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr4_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_04_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:193) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_04_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr4 <- RO_H2_Chr4 %>%
  filter(CHROM == "Urtica_dioica_female_chr_04_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr4, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr4.tsv")

RO_H2_Chr3 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr3_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_03_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:200) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_03_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr3 <- RO_H2_Chr3 %>%
  filter(CHROM == "Urtica_dioica_female_chr_03_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr3, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr3.tsv")

RO_H2_Chr2 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr2_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_02_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:227) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_02_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr2 <- RO_H2_Chr2 %>%
  filter(CHROM == "Urtica_dioica_female_chr_02_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr2, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr2.tsv")

RO_H2_Chr1 <- read.table("RepeatOBserver_results/NettleFemaleR5_H2-AT_Chr1_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1,
         Shannon = V2) %>%
  mutate(CHROM = "Urtica_dioica_female_chr_01_H2")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:237) {   #the max. number = empty_windows$max_window_size+1
  # creating a vector to append to
  # data frame
  vec <- c(as.list((i)*window_size), "Urtica_dioica_female_chr_01_H2")
  # assigning this vector to ith row
  data_frame[i, ] <- vec
}
Shannon_mean_per_250kb_window_H2_Chr1 <- RO_H2_Chr1 %>%
  filter(CHROM == "Urtica_dioica_female_chr_01_H2")%>%
  mutate(window = floor(POS/window_size)*window_size,
         window_middle = window + (0.5*window_size)) %>%
  group_by(CHROM, window, window_middle) %>%
  summarize(mean_Shannon = mean(Shannon, na.rm = TRUE))
write_tsv(Shannon_mean_per_250kb_window_H2_Chr1, "Circos_input/Shannon_mean_per_250kb_window_H2_Chr1.tsv")
#--------------#

#e. Centromeric range prediction by RepeatOBserver
### Formatted the "NettleFemaleR8_H1-AT_total_possible_range.txt" so that the output table looks like: 
#chr  (tab)  #start  (tab)  #end

