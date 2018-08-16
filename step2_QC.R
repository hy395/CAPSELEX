library(Biostrings)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
source("../final/scripts/utils/all_functions.R")

##############
# 8mer stats #
##############
# 1. compute 8mer count, consider reverse complement. save to file.
# 2. compute z-score. save to file.
# 3. compute fold-change.
files <- list.files("raw/rdata/experiments", pattern="rds", full.names=T)
zeros1 <- list.files("../data/selex_initial/raw_seqs", full.names = T)
zeros2 <- list.files("../final/raw/rdata/rdata_cycleZero", full.names=T)

list_cnts <- list()
list_z <- list()
# didn't compute FC. 
vocab <- colnames(oligonucleotideFrequency(DNAStringSet(paste0(  rep("A", 8), collapse="")), width=8))
for (i in 1:length(files)) {
  barcode <- gsub(".*_","",gsub("raw/rdata/experiments/|.rds","",files[i]))
  
  # compute 8mer count
  probes <- readRDS(files[i])
  probes <- mclapply(probes, function(x) { # remove reads with N
    Ns <- grep("N",x)
    if (length(Ns)!=0){x <- x[-Ns]}
    return(x)
  }, mc.cores=length(probes), mc.preschedule=F)
  tmp <- lapply(probes, function(x) { # count
    probes2kmers(c(x,reverseComplement(x)), kmer.len=8, cores=8)
  })
  cnts <- do.call(rbind, tmp)
  list_cnts[[i]] <- cnts
  
  # find cycle0 data
  if (length(grep(barcode, zeros1))==1) { # cycle0 comes with CAP-SELEX
    initial <- readDNAStringSet(zeros1[grep(barcode, zeros1)], format="fastq")
  } else if (length(grep(barcode, zeros2))==1) { # cycle0 comes with HT-SELEX
    initial <- readRDS(zeros2[grep(barcode, zeros2)])
  } else { # otherwise, use cycle1 as cycle0, skip if cycle1 doesn't exist
    cycle1 <- gsub("_[0-9]_","_1_",names(probes)[1])
    if (!cycle1%in%names(probes)) next
    initial <- probes[[cycle1]]
  }
  Ns <- grep("N",initial); if (length(Ns)!=0) {initial <- initial[-Ns]} # remove N
  reads <- unique(c(initial, reverseComplement(initial)))
  mm <- mm_train(reads,5)
  p <- mm_pred(vocab, mm)
  
  list_z[[i]] <- enrichment_over_cycle0(cnts, p) # calculate z-score
  print(sprintf("processing %s out of %s experiments...",i , length(files)))
}
names(list_cnts) <- names(list_z) <- files
saveRDS(list_cnts, "raw/rdata/cnt_8mer.rds")
saveRDS(list_z, "raw/rdata/z_8mer.rds")

########
# plot #
########

##########################
# b/w cycles correlation #
##########################
# by fold-change
dir.create("raw/corr_bw_cycles_z")
z_list <- readRDS("raw/rdata/z_8mer.rds")
corrs <- vector()

for (i in 1:length(z_list)) {
  # criteria for open a new figure
  if (i%%25 == 1) {
    png(paste0("raw/corr_bw_cycles_z/results_",i,".png"), 1600, 1400)
    par(mfrow=c(5,5))
  }
  # plot
  cycle2 <- grep("_2",rownames(z_list[[i]]))
  cycle3 <- grep("_3",rownames(z_list[[i]]))
  if (length(cycle2)==1&length(cycle3)==1) {
    tmp <- t(z_list[[i]][c(cycle2,cycle3),])
    select <- tmp[,1] > 10 & tmp[,2] > 10
    corrs[i] <- cor(tmp[select,1], tmp[select,2], method="spearman")
    plot(tmp, main=paste0(sum(select)," overlap kmers with correlation ",round(corrs[i], 2)), pch=20)
  } else {
    plot(1:10, 1:10)
  }
  # criteria for closing a figure
  if (i%%25 == 0| i==length(z_list)) {
    dev.off()
  }
  print(i)
}
dev.off()
names(corrs) <- names(z_list)
png("raw/corr_bw_cycles_z/spearman_corr.png",800,800)
hist(corrs, breaks=100)
dev.off()
saveRDS(corrs, file="raw/corr_bw_cycles_z/spearman_corr.rds")

##############################
# b/w replicates correlation #
##############################
z_list <- readRDS("raw/rdata/z_8mer.rds")
tf <- gsub("raw/rdata/experiments/|_.*","",sub("_",".",names(z_list)))
z_list <- z_list[order(tf)]
tf <- tf[order(tf)]
vocab <- colnames(oligonucleotideFrequency(DNAStringSet(paste0(  rep("A", 8), collapse="")), width=8))
unique_kmers <- keep1Strand(vocab)
p_list <- list()
correlation <- data.frame(x=rep("",length(tf)), y=rep("",length(tf)), cor=rep(0,length(tf)), stringsAsFactors = F)
for (i in 1:(length(tf)-1)) {
  tmp <- tf[i]
  if (tolower(tmp)!=tolower(tf[i+1])) next
  
  x <- z_list[[i]][grep("_3", rownames(z_list[[i]])),]
  y <- z_list[[i+1]][grep("_3", rownames(z_list[[i+1]])),]
  
  x_name <- rownames(z_list[[i]])[grep("_3", rownames(z_list[[i]]))]
  y_name <- rownames(z_list[[i+1]])[grep("_3", rownames(z_list[[i+1]]))]
  labels <- union(names(sort(x, decreasing=T)[1:20]), names(sort(y, decreasing=T)[1:20]))
  labels <- labels[labels%in%unique_kmers]
  
  toplot <- data.frame(x,y)[unique_kmers,]
  toplot$kmer <- rownames(toplot)
  toplot$color <- factor(ifelse(toplot$kmer%in%labels,1,0))
  
  top <- toplot[toplot[,1] > 10 | toplot[,2] > 10,]
  correlation[i,1] <- x_name
  correlation[i,2] <- y_name
  correlation[i,3] <- cor(top[,1],top[,2], method="spearman")
  
  p <- ggplot(toplot, aes(x=x,y=y, colour=color)) +
    scale_color_manual(values=c("black","red")) +
    geom_point() + theme(legend.position="none") +
    geom_text_repel(data = toplot[labels,], aes(label = kmer, colour=factor(0))) +
    labs(x=x_name, y=y_name, title=paste0(tmp,"_vs_", tf[i+1],": ",round(correlation[i,3],3)))
  
  p_list[[i]] <- p
  print(i)
}
# draw
dir.create("raw/corr_bw_replicates_z")
p_list <- p_list[!sapply(p_list, is.null)]
for (i in 1:5) {
  start <- (i-1)*9 + 1
  end <- i*9
  png(paste0("raw/corr_bw_replicates_z/plot_",i,".png"), 2000, 1600)
  do.call("grid.arrange", c(p_list[start:min(end, length(p_list))], ncol=3, nrow=3))
  dev.off()
  print(i)
}
correlation <- correlation[!correlation[,1]=="",]
saveRDS(correlation, "raw/corr_bw_replicates_FC/spearman_cors.rds")

########################################
# b/w replicates top probe consistency #
########################################
files <- list.files("raw/rdata/experiments", full.names=T)
tf <- gsub("raw/rdata/experiments/|_.*","",sub("_",".",files))
files <- files[order(tf)]
tf <- tf[order(tf)]
vocab <- colnames(oligonucleotideFrequency(DNAStringSet(paste0(  rep("A", 8), collapse="")), width=8))
unique_kmers <- keep1Strand(vocab)
p_list <- list()
consistency <- data.frame(x=rep("",length(tf)), y=rep("",length(tf)), cor=rep(0,length(tf)), stringsAsFactors = F)
for (i in 1:(length(tf)-1)) {
  print(i)
  
  tmp <- tf[i]
  if (tolower(tmp)!=tolower(tf[i+1])) next
  
  # read probes
  a <- readRDS(files[i])
  b <- readRDS(files[i+1])
  a <- a[[grep("_3", names(a))]]
  b <- b[[grep("_3", names(b))]]
  
  # take the top 2000 enriched probes from each
  aaa <- tabdt(as.character(c(a, reverseComplement(a))))
  bbb <- tabdt(as.character(c(b, reverseComplement(b))))
  a_select <- keep1Strand(aaa$x[order(aaa$Freq, decreasing=T)[1:10000]])
  b_select <- keep1Strand(bbb$x[order(bbb$Freq, decreasing=T)[1:10000]])
  
  a_select <- DNAStringSet(a_select[1:min(2000, length(a_select))])
  b_select <- DNAStringSet(b_select[1:min(2000, length(b_select))])
  
  c1 <- colSums(oligonucleotideFrequency(c(a_select, reverseComplement(a_select)),8))
  c2 <- colSums(oligonucleotideFrequency(c(b_select, reverseComplement(b_select)),8))
  
  x_name <- gsub("raw/rdata/experiments/|.rds","",files[i])
  y_name <- gsub("raw/rdata/experiments/|.rds","",files[i+1])
  
  # compute consistency
  x_top <- names(sort(c1, decreasing=T)[1:100])
  y_top <- names(sort(c2, decreasing=T)[1:100])
  consistency[i,1] <- x_name
  consistency[i,2] <- y_name
  consistency[i,3] <- round(length(intersect(x_top, y_top)) / length(union(x_top, y_top)), 3)
  
  # plot
  labels <- union(names(sort(c1, decreasing=T)[1:20]), names(sort(c2, decreasing=T)[1:20]))
  labels <- labels[labels%in%unique_kmers]
  toplot <- data.frame(c1,c2)[unique_kmers,]
  toplot$kmer <- rownames(toplot)
  toplot$color <- factor(ifelse(toplot$kmer%in%labels,1,0))
  toplot <- data.frame(c1,c2)[unique_kmers,]
  toplot$kmer <- rownames(toplot)
  toplot$color <- factor(ifelse(toplot$kmer%in%labels,1,0))
  p <- ggplot(toplot, aes(x=c1,y=c2, colour=color)) +
    scale_color_manual(values=c("black","red")) +
    geom_point() + theme(legend.position="none") +
    geom_text_repel(data = toplot[labels,], aes(label = kmer, colour=factor(0))) +
    labs(x=x_name, y=y_name, title=paste0(tmp,"_vs_", tf[i+1], ": ", consistency[i,3]))
  p_list[[i]] <- p
}
consistency <- consistency[consistency[,1]!="",]
dir.create("raw/replicates_2000probes")
saveRDS(consistency, "raw/replicates_2000probes/consistency.rds")

# draw
p_list <- p_list[!sapply(p_list, is.null)]
for (i in 1:5) {
  start <- (i-1)*9 + 1
  end <- i*9
  png(paste0("raw/replicates_2000probes/plot_",i,".png"), 2000, 1600)
  do.call("grid.arrange", c(p_list[start:min(end, length(p_list))], ncol=3, nrow=3))
  dev.off()
  print(i)
}
