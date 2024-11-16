library(tidyverse)
library(dplyr)
library(ggplot2)

##### 2. EDTA TE density as histogram plot in Circos #####

# Makes sure all windows are in the dataframe, with the 0 hits # 
window_size = 500000 #500kb window is probably best for resolution in Circos for both haplotypes

# H1 ----------------------------------------
empty_windows <- read.table("Nettle_female_H1v2_karyotype.txt") %>%
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
  filter(CHR == "9") %>%
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
  filter(CHR == "8") %>%
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
  filter(CHR == "7") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_07")
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
  filter(CHR == "6") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_06")
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
  filter(CHR == "5") %>%
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
  filter(CHR == "4") %>%
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
  filter(CHR == "3") %>%
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
  filter(CHR == "2") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_02")
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
  filter(CHR == "1") %>%
  mutate(CHROM = "Urtica_dioica_female_chr_01")
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
TE_per_500kb_window_H1_Chr1 <- TE_H1_Chr1 %>%
  mutate(window = floor(POS/window_size)*window_size,
         window_end = window + (window_size - 1)) %>%
  group_by(CHROM, window, window_end) %>%
  summarize(n_TEs = n())
TE_per_500kb_window_H1_Chr1_allwindows <- left_join(Chr1_500kb_windows, TE_per_500kb_window_H1_Chr1) %>%
  replace(is.na(.), 0)
write_tsv(TE_per_500kb_window_H1_Chr1_allwindows, "Circos_input/TE_per_500kb_window_H1_Chr1.tsv")

# H2 ----------------------------------------
empty_windows_H2 <- read.table("Nettle_female_H2v2_karyotype.txt") %>%
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
for(i in 1:53) {   #the max. number = empty_windows$max_window_size+1
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
for(i in 1:71) {   #the max. number = empty_windows$max_window_size+1
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
for(i in 1:91) {   #the max. number = empty_windows$max_window_size+1
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
for(i in 1:86) {   #the max. number = empty_windows$max_window_size+1
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
for(i in 1:96) {   #the max. number = empty_windows$max_window_size+1
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
for(i in 1:118) {   #the max. number = empty_windows$max_window_size+1
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

