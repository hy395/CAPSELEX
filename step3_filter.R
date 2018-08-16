###################
# probe filtering #
###################
library(Biostrings)
library(data.table)
library(parallel)
source("../final/scripts/utils/all_functions.R")

files <- list.files("raw/rdata/experiments", pattern="rds", full.names=T)
zeros1 <- list.files("../data/selex_initial/raw_seqs", full.names = T)
zeros2 <- list.files("../final/raw/rdata/rdata_cycleZero", full.names=T)

outdir <- "filter/"
dir.create(outdir)
dir.create(paste0(outdir, "seqs_2k/"))
dir.create(paste0(outdir, "seqs_neg/"))

#######################################
# remove experiments
# 1. remove experiments with inconsistency enrichment b/w cycles
# 2. remove experiments without cycle 3 (no enrichment)
# 3. remove experiments with < 500 probes occuringing > 2 times (no enrichment)
# 4. remove experiments that: most common 8mer in the final 2000 probes occuring < 100 times (no enrichment)
# remove probes
# 1. remove probes with dustScore > 2 (low complexity)
# 2. top 2k positive probes
#######################################
corrs <-readRDS("raw/corr_bw_cycles_z/spearman_corr.rds")

for (i in 1:length(files)){
  print(paste0("processing #",i,": ", files[i]))
  tmp_cor <- corrs[files[i]]
  if (is.na(tmp_cor) | tmp_cor<0.7) {
    print("failed 0.7 cycle correlation filter, skip")
    next
  }
  
  # decide experiment info
  tf <- gsub("raw/rdata/experiments/|_.*","",sub("_",".",files[i]))
  barcode <- gsub(".*_|.rds","",files[i])
  probes <- readRDS(files[i])
  
  if (length(grep("_3",names(probes)))!=1) {
    print("don't have cycle3, skip")
  }
  
  select <- as.character(probes[[grep("_3",names(probes))]])
  select <- unname(select)
  toremove <- grep("N",select)
  if (length(toremove)!=0){select <- select[-toremove]}
  select <- DNAStringSet(select)
  probeset <- c(select, reverseComplement(select))
  if (sum(tabdt(as.character(probeset))$Freq>2) < 500) {
    print("don't have enough enriched sequences, skip")
    next
  }
  
  # initial cycle
  if (length(grep(barcode, zeros1))==1) { # cycle0 comes with CAP-SELEX
    initial <- readDNAStringSet(zeros1[grep(barcode, zeros1)], format="fastq")
  } else if (length(grep(barcode, zeros2))==1) { # cycle0 comes with HT-SELEX
    initial <- readRDS(zeros2[grep(barcode, zeros2)])
  } else { # otherwise, use cycle1 as cycle0, skip if cycle1 doesn't exist
    cycle1 <- gsub("_[0-9]_","_1_",names(probes)[1])
    if (!cycle1%in%names(probes)) {
      print("don't have cycle0 or cycle1, skip")
      next
      }
    initial <- probes[[cycle1]]
  }
  toremove <- grep("N", initial)
  if (length(toremove)!=0){initial <- initial[-toremove]}
  initial <- unique(c(initial, reverseComplement(initial)))
  
  # compute z-score
  a <- probes2kmers_longk(probeset, initial, 5)
  z <- enrichment_over_cycle0(a[1,,drop=F], a[2,])
  combined <- rbind(a,z=z[1,])
  combined <- combined[,order(combined[3,], decreasing=T)]
  
  # filter probe
  a <- combined[,combined[1,]>=2] # probes have at least 2 counts
  a <- a[,dustScore(DNAStringSet(colnames(a))) < 2] # probes dust core < 2
  
  # save positive probes
  a <- a[,1:min(10000, ncol(a))]
  pos <- keep1Strand(colnames(a))
  towrite_2k <- DNAStringSet(pos)[1:min(2000, length(pos))]
  
  # final filtering: if most enriched kmer occurs < 100 times. Remove that experiment.
  tmp <- probes2kmers(c(towrite_2k, reverseComplement(towrite_2k)), 8)
  if(max(tmp) < 100) {
    print("final sequences don't enrich for core sequence, skip")
    next
  }
  
  # save 2k positive probes
  writeXStringSet(towrite_2k, paste0(outdir, "seqs_2k/",gsub(".rds","_pos.fasta",gsub(".*/","",files[i]))), format="fasta")
  
  # save 10k negative probes
  neg <- DNAStringSet(sample(initial[!initial%in%probeset],5000))
  writeXStringSet(neg, paste0(outdir, "seqs_neg/",gsub(".rds","_neg.fasta",gsub(".*/","",files[i]))), format="fasta")
}

# #########################
# # check reproducibility #
# #########################
# library(ggplot2)
# library(ggrepel)
# library(grid)
# library(gridExtra)
# source("../final/scripts/utils/all_functions.R")
# 
# files <- list.files(paste0(outdir, "seqs_2k"), full.names=T, pattern="pos")
# tf <- gsub(paste0(outdir, "seqs_2k/|_.*"),"",files)
# vocab <- colnames(oligonucleotideFrequency(DNAStringSet(paste0(  rep("A", 8), collapse="")), width=8))
# unique_kmers <- keep1Strand(vocab)
# p_list <- list()
# consistency <- data.frame(x=rep("",length(tf)), y=rep("",length(tf)), cor=rep(0,length(tf)), stringsAsFactors = F)
# for (i in 1:(length(tf)-1)) {
#   tmp <- tf[i]
#   if (tolower(tmp)!=tolower(tf[i+1])) next
# 
#   a <- readDNAStringSet(files[i])
#   b <- readDNAStringSet(files[i+1])
#   c1 <- probes2kmers(c(a, reverseComplement(a)),8)
#   c2 <- probes2kmers(c(b, reverseComplement(b)),8)
#   x_name <- gsub(paste0(outdir, "seqs_2k/|_pos.fasta"),"",files[i])
#   y_name <- gsub(paste0(outdir, "seqs_2k/|_pos.fasta"),"",files[i+1])
# 
#   # compute consistency
#   x_top <- names(sort(c1, decreasing=T)[1:100])
#   y_top <- names(sort(c2, decreasing=T)[1:100])
#   consistency[i,1] <- x_name
#   consistency[i,2] <- y_name
#   consistency[i,3] <- round(length(intersect(x_top, y_top)) / length(union(x_top, y_top)), 3)
# 
#   # plot
#   labels <- union(names(sort(c1, decreasing=T)[1:20]), names(sort(c2, decreasing=T)[1:20]))
#   labels <- labels[labels%in%unique_kmers]
#   toplot <- data.frame(c1,c2)[unique_kmers,]
#   toplot$kmer <- rownames(toplot)
#   toplot$color <- factor(ifelse(toplot$kmer%in%labels,1,0))
#   p <- ggplot(toplot, aes(x=c1,y=c2, colour=color)) +
#     scale_color_manual(values=c("black","red")) +
#     geom_point() + theme(legend.position="none") +
#     geom_text_repel(data = toplot[labels,], aes(label = kmer, colour=factor(0))) +
#     labs(x=x_name, y=y_name, title=paste0(tmp,"_vs_", tf[i+1], ": ", consistency[i,3]))
#   p_list[[i]] <- p
# 
#   print(i)
# }
# 
# # draw
# dir.create(paste0(outdir, "replicates/"))
# p_list <- p_list[!sapply(p_list, is.null)]
# consistency <- consistency[consistency[,1]!="",]
# for (i in 1:5) {
#   start <- (i-1)*9 + 1
#   end <- i*9
#   png(paste0(outdir, "replicates/plot_",i,".png"), 2000, 1600)
#   do.call("grid.arrange", c(p_list[start:min(end, length(p_list))], ncol=3, nrow=3))
#   dev.off()
#   print(i)
# }
# saveRDS(consistency, paste0(outdir, "replicates/replicates_consistency.rds"))
# 
# #######################
# # compare consistency #
# #######################
# library(ggplot2)
# a <- readRDS("raw/replicates_2000probes/consistency.rds")
# b <- readRDS(paste0(outdir, "replicates/replicates_consistency.rds"))
# 
# pdf(paste0(outdir, "replicates/QC.pdf"), 12, 12)
# boxplot(list(before=a[,3], after=b[,3]), xlab="Quality Control", ylab="consistency")
# stripchart(list(before=a[,3], after=b[,3]), vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
# p <- wilcox.test(a[,3], b[,3], alternative="less")$p.value
# text(1.5, 0.7, labels = paste0("p < ", round(p,3)))
# dev.off()
