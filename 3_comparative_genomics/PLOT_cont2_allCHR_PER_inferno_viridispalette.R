library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(viridis)

#functions
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

filterMum <- function(df, minl=1000, flanks=1e4){
    coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
        summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
        ungroup %>% arrange(desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
    merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
        mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

diagMum <- function(df){
    ## Find best qid order
    rid.o = df %>% group_by(qid, rid) %>% summarize(base=sum(abs(qe-qs)),
                                                    rs=weighted.mean(rs, abs(qe-qs))) %>%
        ungroup %>% arrange(desc(base)) %>% group_by(qid) %>% do(head(., 1)) %>%
        ungroup %>% arrange(desc(rid), desc(rs)) %>%
        mutate(qid=factor(qid, levels=unique(qid)))
    ## Find best qid strand
    major.strand = df %>% group_by(qid) %>%
        summarize(major.strand=ifelse(sum(sign(qe-qs)*abs(qe-qs))>0, '+', '-'),
                  maxQ=max(c(qe, qs)))
    merge(df, major.strand) %>% mutate(qs=ifelse(major.strand=='-', maxQ-qs, qs),
                                       qe=ifelse(major.strand=='-', maxQ-qe, qe),
                                       qid=factor(qid, levels=levels(rid.o$qid)))
}


        args <- commandArgs(trailingOnly = TRUE)

                f <-  args[1] #delta
#                chr <- args[2]
                #bed <- args[3] #trans
        #       bed2 <-  args[2] #inv
        #       bed3 <- args[3] # centromers
                #bed4 <- args[5] #chr centromers
                #st <- args[3]
                #end <- args[4]

         #      centromers <-  read.table(bed3,  header=F)
          #     centromers <- centromers[order(centromers[,1], centromers[,2]), ]
           #    centromers$start <- pmin(centromers$V2, centromers$V3)
            #   centromers$end <- pmax(centromers$V2, centromers$V3)
               #centromers <- subset(centromers, centromers$V1== bed4)
             #  centromers_s <- centromers$start
               #TRANS <- read.table(bed,  header=F)
               #TRANS_L <- TRANS$V7
               #TRANS_R <- TRANS$V8

              # INV <- read.table(bed2,  header=F)
               #INV <- subset(INV, INV$V7 == chr)
               #INV$lenght <- INV$V2 - INV$V1
               #INV_long <- subset(INV, INV$lenght > 100000)

               #INV_L <- INV_long$V1
               #INV_R <- INV_long$V2
               #INV_L <- INV$V7
               #INV_R <- INV$V8

                 name <- gsub(".txt", "", f)

                # name2 <-  paste(name, chr, "Ha412", ".jpg" , sep="_")
                 #mumgp = readDelta(f)
                 mumgp <- read.table(f, header=F)
                 colnames(mumgp) <- c("rs", "re", "qs", "qe", "sim", "qid", "rid", "strand")
                 mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                 if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$sim, mumgp$qid, mumgp$rid, mumgp$strand))
                 #mumgp$qid <- paste(mumgp$qid, "Arg", sep="_")

                 #mumgp <- read.table(f, header=F)
                 #colnames(mumgp) <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
                 #mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                  #  if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$error, mumgp$qid, mumgp$rid, mumgp$strand))

 #                mumgp <- subset(mumgp, mumgp$rid == chr) #& mumgp$rs > st & mumgp$re < end )
                # mumgp <- subset(mumgp, mumgp$sim > 95) #I removed thisfor anchorwave
                 mumgp.filt = filterMum(mumgp, minl=1e4)
                 mumgp.filt.diag = diagMum(mumgp.filt)

# If you want to modify y-axis order as you wish,
  desired_order <- c("F_H1_chr_01", "F_H2_chr_01", "M_H1_chr_01", "M_H2_chr_01", "F_H1_chr_02", "F_H2_chr_02", "M_H1_chr_02", "M_H2_chr_02", "F_H1_chr_03", "F_H2_chr_03", "M_H1_chr_03", "M_H2_chr_03", "F_H1_chr_04", "F_H2_chr_04", "M_H1_chr_04", "M_H2_chr_04", "F_H1_chr_05", "F_H2_chr_05", "M_H1_chr_05", "M_H2_chr_05", "F_H1_chr_06", "F_H2_chr_06", "M_H1_chr_06", "M_H2_chr_06", "F_H1_chr_07", "F_H2_chr_07", "M_H1_chr_07", "M_H2_chr_07", "F_H1_chr_08", "F_H2_chr_08", "M_H1_chr_08", "M_H2_chr_08", "F_H1_chr_09", "F_H2_chr_09", "M_H1_chr_09", "M_H2_chr_09", "F_H1_chr_10", "F_H2_chr_10", "M_H1_chr_10", "M_H2_chr_10", "F_H1_chr_11", "F_H2_chr_11", "M_H1_chr_11", "M_H2_chr_11", "F_H1_chr_12", "F_H2_chr_12", "M_H1_chr_12", "M_H2_chr_12", "F_H1_chr_13", "F_H2_chr_13", "M_H1_chr_13", "M_H2_chr_13")
  mumgp.filt.diag$qid <- factor(mumgp.filt.diag$qid, levels = desired_order)

              P1 <- ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=sim)) +
  geom_segment(show.legend=FALSE, size=3) +
                      geom_point(alpha=0.09) + theme_bw() +
  scale_y_discrete(drop = FALSE) +
  facet_grid(qid~rid, scales='free', space='free', switch='both') +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(strip.text.y=element_text(angle=180, size=10),
        strip.text.x=element_text(size=10),
        strip.background=element_blank(),
        legend.position=c(1,-.03), legend.justification=c(1,1),
        legend.direction='horizontal',
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing=unit(0, 'cm')) +
  xlab('reference sequence') + ylab('assembly')  +  scale_colour_gradient2(low = inferno(1), high = inferno(50), midpoint = 97.5)

jpeg( paste(name, ".jpg", sep ="_"),  width=4000, height= 4000)
print(P1)
dev.off()
