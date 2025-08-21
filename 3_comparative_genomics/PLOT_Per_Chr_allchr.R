library(dplyr)
library(magrittr)
library(GenomicRanges)
library(knitr)
library(ggplot2)
library(tidyr)
library(tidyverse)

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
               chr <- args[2]
                #bed <- args[3] #trans
        #        bed2 <-  args[3] #inv
        #       bed3 <- args[4] # centromers
        #       bed4 <- args[5] #chr centromers
                #st <- args[3]
                #end <- args[4]

         #      centromers <-  read.table(bed3,  header=F)
          #     centromers <- centromers[order(centromers[,1], centromers[,2]), ]
          #     centromers$start <- pmin(centromers$V2, centromers$V3)
          #     centromers$end <- pmax(centromers$V2, centromers$V3)
          #     centromers <- subset(centromers, centromers$V1== bed4)
          #     centromers_s <- centromers$start
               #TRANS <- read.table(bed,  header=F)
               #TRANS_L <- TRANS$V7
               #TRANS_R <- TRANS$V8

         #      INV <- read.table(bed2,  header=F)
         #      INV <- subset(INV, INV$V7 == chr)
         #      INV$lenght <- INV$V2 - INV$V1
         #      INV_long <- subset(INV, INV$lenght > 100000)

         #      INV_L <- INV_long$V1
         #      INV_R <- INV_long$V2
               #INV_L <- INV$V7
               #INV_R <- INV$V8

                 name <- gsub(".txt", "", f)

                # name2 <-  paste(name, chr, "Ha412", ".jpg" , sep="_")
                 #mumgp = readDelta(f)
                 mumgp <- read.table(f, header=F)
                 colnames(mumgp) <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
                 mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                 if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$error, mumgp$qid, mumgp$rid, mumgp$strand))
                 #mumgp$qid <- paste(mumgp$qid, "Arg", sep="_")

                 #mumgp <- read.table(f, header=F)
                 #colnames(mumgp) <- c("rs", "re", "qs", "qe", "error", "qid", "rid", "strand")
                 #mumgp[c("qs", "qe")] <- t(mapply(\(a, b, c, d, e, f, g, h){
                  #  if(h == "-") c(d,c) else c(c,d) }, mumgp$rs, mumgp$re, mumgp$qs, mumgp$qe, mumgp$error, mumgp$qid, mumgp$rid, mumgp$strand))

                 mumgp <- subset(mumgp, mumgp$rid == chr) #& mumgp$rs > st & mumgp$re < end )
                #FILTER PER CHR

                 mumgp$H1_H1_name <- paste(sub("^(H\\d+).*", "\\1", mumgp$rid), round(mumgp$rs/1000000, digits = 2), round(mumgp$re/1000000, digits = 2), sub("^(H\\d+).*", "\\1", mumgp$qid), round(mumgp$qs/1000000, digits = 2), round(mumgp$qe/1000000, digits = 2), sep = "-")

                 mumgp.filt = filterMum(mumgp, minl=1e4)
                 mumgp.filt.diag = diagMum(mumgp.filt)

 #  mumgp.filt.diag$H1_H1_name <- paste(sub("^(H\\d+).*", "\\1", mumgp.filt.diag$rid), round(mumgp.filt.diag$rs/1000000, digits = 1), round(mumgp.filt.diag$re/1000000, digits = 1), sub("^(H\\d+).*", "\\1", mumgp.filt.diag$qid), round(mumgp.filt.diag$qs/1000000, digits = 1), round(mumgp.filt.diag$qe/1000000, digits = 1), sep = "-")


        #add vertical lines to the transposition region
            #  data_vline <- data.frame(group = c("Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr6", "Chr6", "Chr6", "Chr8", "Chr9", "Chr11", "Chr12", "Chr14", "Chr15", "Chr16", "Chr17"), vline
#=c(NA, NA, NA, NA, 125000000, 130000000, 130000000, 155000000, NA, NA, NA, NA, NA, NA, NA, NA))

              P1 <- ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment(show.legend=FALSE, size=3) + geom_point(alpha=0.09) + theme_bw() +
   geom_text(aes(x = re, y = qe, label = H1_H1_name), size = 3) +
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
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')
#+
        # geom_vline(xintercept=INV_R, color = "green", size=1) +
        # geom_vline(xintercept=INV_L, color = "green", size=1) +
#        geom_vline(xintercept=centromers_s, color = "purple", size=1, alpha = .05) +
        # annotate("rect", ymin = 0, ymax =  max( mumgp.filt.diag$qe), xmin = INV_L, xmax = INV_R, alpha = 0.2, fill = "green")
        # geom_vline(xintercept=end, color = "purple", size=1, alpha = .05) +
        #geom_vline(data = mumgp.filt.diag %>% filter(rid == "Ha412HOChr03"), aes(xintercept = 28000000), linetype="dotted", size=2) +
        #geom_vline(data = mumgp.filt.diag %>% filter(rid == "Ha412HOChr03"), aes(xintercept = 42000000), linetype="dotted", size=2) +
        #geom_vline(data = mumgp.filt.diag %>% filter(rid == "Ha412HOChr03"), aes(xintercept = 89000000), linetype="dotted", size=2)
       # geom_hline(yintercept=TRANS_L, color = "green", size=1,  alpha = .05) +
       # geom_hline(yintercept=TRANS_R, color = "green", size=1,  alpha = .05) +
       # geom_hline(yintercept=INV_R, color = "red", size=1, alpha = .2) +
       # geom_hline(yintercept=INV_L, color = "red", size=1, alpha = .2) +
       # annotate("rect", xmin = 0, xmax =  max( mumgp.filt.diag$re), ymin = INV_L, ymax = INV_R, alpha = 0.2, fill = "red")

jpeg( paste(name, chr, ".jpg", sep ="_"),  width=2000, height= 2000)
print(P1)
dev.off()
