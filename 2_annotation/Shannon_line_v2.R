library(tidyverse)
library(dplyr)
library(ggplot2)

##### 1. RepeatOBserver Shannon diversity plot for line plot in Circos #####

# Makes sure all windows are in the dataframe, with the 0 hits # 
window_size = 250000 #250kb window is probably best for resolution in Circos for both haplotypes

# H1 ----------------------------------------
empty_windows <- read.table("Nettle_female_H1v2_karyotype.txt") %>%
  as_tibble(.)%>% 
  rename(CHROM = V1, length = V2) %>% 
  mutate(max_window_size = floor(length/window_size))

## Do this per chromosome and adjust the window number according to "empty_windows" ##
# refer to the dataframe_prep_for_het_plots.R
RO_H1_Chr13 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr13_Shannon_div.txt", header = FALSE) %>%
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

RO_H1_Chr12 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr12_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_12")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:95) {   #the max. number = empty_windows$max_window_size+1                   
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

RO_H1_Chr11 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr11_Shannon_div.txt", header = FALSE) %>%
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

RO_H1_Chr10 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr10_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_10")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:108) {   #the max. number = empty_windows$max_window_size+1                   
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

RO_H1_Chr9 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr9_Shannon_div.txt", header = FALSE) %>%
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

RO_H1_Chr8 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr8_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_08")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:171) {   #the max. number = empty_windows$max_window_size+1                   
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

RO_H1_Chr7 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr7_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_07")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:190) {   #the max. number = empty_windows$max_window_size+1                   
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

RO_H1_Chr6 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr6_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_06")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:176) {   #the max. number = empty_windows$max_window_size+1                   
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

RO_H1_Chr5 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr5_Shannon_div.txt", header = FALSE) %>%
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

RO_H1_Chr4 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr4_Shannon_div.txt", header = FALSE) %>%
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

RO_H1_Chr3 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr3_Shannon_div.txt", header = FALSE) %>%
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

RO_H1_Chr2 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr2_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_02")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:218) {   #the max. number = empty_windows$max_window_size+1                   
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

RO_H1_Chr1 <- read.table("RepeatOBserver_results/NettleFemaleR8_H1-AT_Chr1_Shannon_div.txt", header = FALSE) %>%
  rename(POS = V1, 
         Shannon = V2) %>% 
  mutate(CHROM = "Urtica_dioica_female_chr_01")
data_frame = data.frame(
  window_1 = numeric(), CHROM = character(),stringsAsFactors = FALSE) %>% 
  as.tibble()
for(i in 1:232) {   #the max. number = empty_windows$max_window_size+1                   
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

