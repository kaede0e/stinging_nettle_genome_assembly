library(tidyverse)
library(dplyr)
library(ggplot2)

##### 3. BRAKER3 gene density as histogram plot in Circos #####

# Makes sure all windows are in the dataframe, with the 0 hits #
window_size = 500000 #500kb window is probably best for resolution in Circos for both haplotypes

# H1 ----------------------------------------
#empty_windows <- read.table("Nettle_female_H1_karyotype.txt") %>%
#  as_tibble(.)%>%
#  rename(CHROM = V1, length = V2) %>%
#  mutate(max_window_size = floor(length/window_size)) # Determine max number of windows per chromosome and change below.

## Do this per chromosome and adjust the window number according to "empty_windows" ##
# refer to the dataframe_prep_for_het_plots.R
Gene_H1_Chr13 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H1_Chr12 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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
Gene_per_500kb_window_H1_Chr12 <- Gene_H1_Chr12 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr12_allwindows <- left_join(Chr12_500kb_windows, Gene_per_500kb_window_H1_Chr12) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr12_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr12.tsv")

Gene_H1_Chr11 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H1_Chr10 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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
Gene_per_500kb_window_H1_Chr10 <- Gene_H1_Chr10 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H1_Chr10_allwindows <- left_join(Chr10_500kb_windows, Gene_per_500kb_window_H1_Chr10) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H1_Chr10_allwindows, "Circos_input/Gene_per_500kb_window_H1_Chr10.tsv")

Gene_H1_Chr9 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H1_Chr8 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H1_Chr7 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_07")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:95) {   #the max. number = empty_windows$max_window_size+1
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

Gene_H1_Chr6 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_06")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:88) {   #the max. number = empty_windows$max_window_size+1
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

Gene_H1_Chr5 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H1_Chr4 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H1_Chr3 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H1_Chr2 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_02")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:109) {   #the max. number = empty_windows$max_window_size+1
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

Gene_H1_Chr1 <- read.table("braker3.Nettle_female_H1_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_01")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:116) {   #the max. number = empty_windows$max_window_size+1
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
Gene_H2_Chr13 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H2_Chr12 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H2_Chr11 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H2_Chr10 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_10")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:53) {   #the max. number = empty_windows$max_window_size+1
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

Gene_H2_Chr9 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H2_Chr8 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_08")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:71) {   #the max. number = empty_windows$max_window_size+1
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

Gene_H2_Chr7 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H2_Chr6 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_06")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:91) {   #the max. number = empty_windows$max_window_size+1
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

Gene_H2_Chr5 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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
Gene_per_500kb_window_H2_Chr5 <- Gene_H2_Chr5 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr5_allwindows <- left_join(Chr5_500kb_windows, Gene_per_500kb_window_H2_Chr5) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr5_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr5.tsv")

Gene_H2_Chr4 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
  select(V1,V2,V3) %>%
  rename(CHROM = V1,
         POS = V2,
         POS_END = V3) %>%
  filter(CHROM == "Urtica_dioica_female_chr_04")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>%
  as.tibble()
for(i in 1:96) {   #the max. number = empty_windows$max_window_size+1
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

Gene_H2_Chr3 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H2_Chr2 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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

Gene_H2_Chr1 <- read.table("braker3.Nettle_female_H2_v2_RNA_ProtsViridiplantaeOrthoDB.bed", header = FALSE) %>%
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
Gene_per_500kb_window_H2_Chr1 <- Gene_H2_Chr1 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_Genes = n())
Gene_per_500kb_window_H2_Chr1_allwindows <- left_join(Chr1_500kb_windows, Gene_per_500kb_window_H2_Chr1) %>%
  replace(is.na(.), 0)
write_tsv(Gene_per_500kb_window_H2_Chr1_allwindows, "Circos_input/Gene_per_500kb_window_H2_Chr1.tsv")
# afterwards, you have to change the first column "Urtica_dioica_female_chr_01" to "Urtica_dioica_female_chr_01_H2" to be compatible with my Circos plot. 

