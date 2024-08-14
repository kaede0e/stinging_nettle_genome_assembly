#this script is a modification from a Mojtaba script
library(tidyverse)
library(tidysq)
library(data.table)
args = commandArgs(trailingOnly = TRUE)

MIX_ASSEM <- args[1]
PREFIX <- args[2]
MIX_FASTA <- args[3]
SAVE_DIR <- args[4]

HAP12_assembly <- read.table(MIX_ASSEM, fill = TRUE, col.names = paste("V", 1:200, sep = ""))
read_fasta(MIX_FASTA, alphabet = "dna_ext") -> HAP12_FASTA

#fragmented contigs
HAP12_assembly %>%
  filter(grepl("^>",V1)) %>%
  select(fragment_ID = V1,
         original_order = V2,
         length = V3)   %>%
  separate(fragment_ID,
           into = c("contig","lable1","lable2"),
           sep = ":::",
           remove=F) %>%
  filter(grepl("fragment",lable1)) ->  FRAGMENTS

#Run if some of the contigs are fragmented
if (nrow(FRAGMENTS) > 0 ) {

#calculate the coordination of the fragmented contigs in the original contig
FRAGMENTS %>%
  separate(lable1,
           into = c("Fragment","Fragment_order")) %>%
  group_by(contig) %>%
  arrange(as.numeric(Fragment_order)) %>%
  mutate(end=cumsum(length)) %>%
  ungroup() %>%
  #mutate(end = ifelse(Fragment_order == 1, length, lag(length) + length)) %>%
  mutate(start = (end-length) + 1) %>%
  mutate(fragment_ID = gsub(">","",fragment_ID)) %>%
  mutate(contig = gsub(">","",contig)) %>%
  select(contig,
         fragment_ID,
         original_order,
         start,
         end,
         length) -> coord_fragments

#extract hw sequence of the fragmented contigs from the original contig
FRAGMENTED_CONTIGS <- NULL
for (CONTIG in 1:nrow(coord_fragments)) {
  HAP12_FASTA %>%
    filter(name == as.character(coord_fragments[CONTIG,1])) %>%
    pull(sq) %>%
    bite(as.numeric(coord_fragments[CONTIG,4]):as.numeric(coord_fragments[CONTIG,5])) %>%
    enframe() %>%
    mutate(name = as.character(coord_fragments[CONTIG,2])) %>%
    select(sq = value,
           name) %>%
    bind_rows(.,
          FRAGMENTED_CONTIGS) -> FRAGMENTED_CONTIGS
}

#remove the the original contig and replace the fragments
HAP12_FASTA %>%
  filter(!name %in%  pull(distinct(coord_fragments,contig))) %>%
  bind_rows(.,
            FRAGMENTED_CONTIGS) -> HAP12_FASTA_FRAG

      rm(HAP12_FASTA, #comente lo de rm
        coord_fragments,
        FRAGMENTED_CONTIGS,
        CONTIG)
      #HAP12_FASTA_FRAG -> HAP12_FASTA_FRAG #aquí está la diferencia , debería ser HAP12_FASTA_FRAG -> HAP12_FASTA así qu elo voy a corregir debajo
      HAP12_FASTA_FRAG -> HAP12_FASTA
      rm(HAP12_FASTA_FRAG)
}else {
  print("NO fragmented contig")
}

rm(FRAGMENTS)
#acá voy
#Extract chromosome info for assembly file



#no fragments
HAP12_assembly %>%
  filter(!grepl("^>",V1)) %>%
  mutate(CHR = row_number()) %>%
  gather(contig_order_in_chr, order_sign, V1:V200) %>% #cambie V1 a V200
  filter(!is.na(order_sign)) %>%
  mutate(contig_order_in_chr=as.numeric(gsub("V","",contig_order_in_chr))) %>%
  mutate(order = abs(as.numeric(order_sign))) %>%
  arrange(CHR,
          contig_order_in_chr)  %>%
  full_join(.,
            select(filter(HAP12_assembly,grepl("^>",V1)),fragment_ID = V1,original_order = V2),
            by=c("order"="original_order")) %>%
  mutate(fragment_ID = gsub(">","",fragment_ID)) %>%
  mutate(HAP = ifelse(grepl("h1",fragment_ID),"H1","H2")) %>%
  filter(!grepl("fragment",  fragment_ID)) %>%
  mutate(Chromosome = paste0(HAP, "_HiC_scaffold_", CHR)) %>%
  mutate(orientation = ifelse(order_sign >= 0, "F", "R")) %>%
  select(Chromosome, contig_order_in_chr, fragment_ID, orientation)  -> mid_data

#fragments no debris
HAP12_assembly %>%
  filter(!grepl("^>",V1)) %>%
  mutate(CHR = row_number()) %>%
  gather(contig_order_in_chr, order_sign, V1:V200) %>% #cambie V1 a V200
  filter(!is.na(order_sign)) %>%
  mutate(contig_order_in_chr=as.numeric(gsub("V","",contig_order_in_chr))) %>%
  mutate(order = abs(as.numeric(order_sign))) %>%
  arrange(CHR,
          contig_order_in_chr)  %>%
  full_join(.,
            select(filter(HAP12_assembly,grepl("^>",V1)),fragment_ID = V1,original_order = V2, size= V3),
            by=c("order"="original_order")) %>%
  mutate(fragment_ID = gsub(">","",fragment_ID)) %>%
  mutate(HAP = ifelse(grepl("h1",fragment_ID),"H1","H2")) %>%
  filter(grepl("fragment", fragment_ID))  %>%
  filter(!grepl("debris",  fragment_ID))  -> mid_data2

#debris
HAP12_assembly %>%
  filter(!grepl("^>",V1)) %>%
  mutate(CHR = row_number()) %>%
  gather(contig_order_in_chr, order_sign, V1:V200) %>% #cambie V1 a V200
  filter(!is.na(order_sign)) %>%
  mutate(contig_order_in_chr=as.numeric(gsub("V","",contig_order_in_chr))) %>%
  mutate(order = abs(as.numeric(order_sign))) %>%
  arrange(CHR,
          contig_order_in_chr)  %>%
  full_join(.,
            select(filter(HAP12_assembly,grepl("^>",V1)),fragment_ID = V1,original_order = V2),
            by=c("order"="original_order")) %>%
  mutate(fragment_ID = gsub(">","",fragment_ID)) %>%
  mutate(HAP = ifelse(grepl("h1",fragment_ID),"H1","H2")) %>%
  filter(grepl("debris", fragment_ID)) %>%
   mutate(Chromosome = paste0(HAP, "_HiC_scaffold_", CHR)) %>%
  mutate(orientation = ifelse(order_sign >= 0, "F", "R")) %>%
  select(Chromosome, contig_order_in_chr, fragment_ID, orientation) -> debris


#create table for fragments
mid_data2 %>% group_by(CHR) %>% filter(size == max(size)) %>% ungroup() %>% mutate(Chromosome = paste0(HAP, "_HiC_scaffold_", CHR)) -> mid_data2_2

chromosome_df <- mid_data2_2 %>% select(CHR, Chromosome) %>% distinct()

mid_data2  %>% left_join(chromosome_df, by = "CHR") %>% mutate(orientation = ifelse(order_sign >= 0, "F", "R"))  %>%  select(Chromosome, contig_order_in_chr, fragment_ID, orientation) -> mid_data2


#merge tables
#bind_rows(mid_data, mid_data2, debris) -> CONTIG_ORDER_CHR
if( nrow(mid_data2) > 0) {
bind_rows(mid_data, mid_data2, debris) -> CONTIG_ORDER_CHR
} else {
mid_data -> CONTIG_ORDER_CHR
}



 CONTIG_ORDER_CHR %>% mutate(numeric_part = as.integer(sub(".*_scaffold_(\\d+)", "\\1", Chromosome)) ) %>% mutate(HAP = str_extract(Chromosome, "H[12]")) %>%  arrange(HAP,  numeric_part, contig_order_in_chr,  numeric_part, contig_order_in_chr ) %>% select(Chromosome, contig_order_in_chr, fragment_ID, orientation, numeric_part, HAP) -> CONTIG_ORDER_CHR

CONTIG_ORDER_CHR %>% group_by(HAP)  %>% arrange(HAP, numeric_part) %>%  mutate(NEW_CHR = as.integer(factor(numeric_part))) %>%  ungroup() %>% mutate(Chromosome = paste(HAP,"HiC_scaffold",NEW_CHR,sep = "_")) %>% select(Chromosome, contig_order_in_chr, fragment_ID, orientation) ->  CONTIG_ORDER_CHR

#estoy ya se ve bien, solo comparar con lo que tenía antes y ver si es igual en cuanto a numero de rows y así

#mid_data %>%
 # group_by(CHR,HAP) %>%
 # tally() %>%
 # ungroup() %>%
 # group_by(CHR) %>%
 # filter(n==max(n)) %>%
 # ungroup() %>%
 # select(CHR,HAP_new = HAP) -> CHRHAP

#mid_data %>%
 # full_join(.,
  #          CHRHAP) %>%
 # group_by(HAP_new) %>%
 #   arrange(HAP_new,CHR) %>%
 #   mutate(NEW_CHR = as.integer(factor(CHR))) %>%
 # ungroup() %>%
 # mutate(Chromosome = paste(HAP_new,"HiC_scaffold",NEW_CHR,sep = "_")) %>%
 # mutate(orientation = ifelse(order_sign<0,"R","F")) %>%
 # select(Chromosome,
  #       contig_order_in_chr,
  #       fragment_ID,
  #       orientation) -> CONTIG_ORDER_CHR

HAP12_FASTA %>%  #pero lo había borrado antes,por eso no salía
  full_join(.,
            CONTIG_ORDER_CHR,by=c("name" = "fragment_ID")) -> REV_FOR_CONTIG


 REV_FOR_CONTIG %>%
  filter(orientation == "R") %>%
  select(-orientation) -> REV_CONTIG


REV_FOR_CONTIG %>%
  filter(orientation == "F") %>%
  select(-orientation) -> FOR_CONTIG


REV_CONTIG %>% mutate(sq = reverse(complement(sq))) -> REV_COMP_CONTIGS


bind_rows(FOR_CONTIG,
          REV_COMP_CONTIGS) -> CORRECTED_CONTIGS


CORRECTED_CONTIGS %>%
  distinct(Chromosome) %>%
  pull(Chromosome) -> CHROMOSOME

NS <- strrep("N", 500)

FINAL_FASTA <- NULL
for (CHR in 1:length(CHROMOSOME)) { CORRECTED_CONTIGS %>% filter(Chromosome == as.character(CHROMOSOME[CHR])) %>%  arrange(contig_order_in_chr) %>% pull(sq) %>% paste(  . ,sq(c(rep(NS,  length(.) - 1), ""), "dna_ext"))  %>%  collapse() %>% enframe() %>% mutate(name = as.character(CHROMOSOME[CHR]))  %>% select(sq = value,  name) %>% bind_rows(., FINAL_FASTA) -> FINAL_FASTA }

#ordenaring contigs
FINAL_FASTA %>%
  mutate(length = get_sq_lengths(sq)) %>%
  arrange(desc(length)) %>%
  select(-length)-> FINAL_FASTA

write_fasta(pull(filter(FINAL_FASTA,grepl("^H1",name)),sq),pull(filter(FINAL_FASTA,grepl("^H1",name)),name),paste0(SAVE_DIR,"/",PREFIX,"_hap1.reviewed.chr_assembled.fasta"))

write_fasta(pull(filter(FINAL_FASTA,grepl("^H2",name)),sq),pull(filter(FINAL_FASTA,grepl("^H2",name)),name),paste0(SAVE_DIR,"/",PREFIX,"_hap2.reviewed.chr_assembled.fasta"))
