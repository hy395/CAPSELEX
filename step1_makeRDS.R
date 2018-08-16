library(Biostrings)
source("../final/scripts/utils/all_functions.R")

###################################
# exps pass QC by Jolma. Table S2 #
###################################
# note: barcode for MEIS2_HOXA13 is wrong. A T is missing at the beginning. I changed it in the txt file.
qc_pass <- read.table("raw/Table_S2.txt", stringsAsFactors = F, sep="\t")
qc_pass <- unique(qc_pass[seq(from=1, to=nrow(qc_pass), by=5),])
qc_pass$pattern <- ""
for (i in 1:nrow(qc_pass)) {
  tmp <- as.character(qc_pass[i, ])
  pattern <- ifelse(length(grep("_",tmp[1])!=0), paste(tmp[1], "[0-9]", tmp[5], tmp[4], sep="_"), paste(tmp[1], "htSELEX", "[0-9]", tmp[5], tmp[4], sep="_"))
  qc_pass[i, "pattern"] <- pattern
}
qc_pass$exp_name <- gsub("_\\[0-9\\]","",qc_pass$pattern)

# 9 experiments in the table but not on ENA are removed
inputdir <- "../data/CAP-SELEX/raw"
tmp <- which(sapply(qc_pass$pattern, function(x) {
  length(grep(x, list.files(inputdir, recursive=T)))==0
}))
qc_pass <- qc_pass[-tmp, ]
write.table(qc_pass, file="raw/Table_S2_trim.txt")

##############
# create rds #
##############
dir.create("raw/rdata/experiments/")
outdir <- "raw/rdata/experiments/"
files <- list.files(inputdir, pattern="fastq.gz", recursive = T, full.names = T)
depth_table <- matrix(0, nrow(qc_pass), 5, dimnames = list(qc_pass$exp_name, paste0("cycle", 1:5)))
for (i in 1:nrow(qc_pass)) {
  pattern <- qc_pass[i, "pattern"]
  exp_name <- qc_pass[i, "exp_name"]
  exps <- files[grep(pattern, files)]
  
  reads <- mclapply(exps, function(x) readDNAStringSet(x, format="fastq"), mc.cores=length(exps))
  names(reads) <- gsub(".*/","",exps)
  saveRDS(reads, paste0(outdir, exp_name, ".rds"))  
  
  tmp_depth <- sapply(reads, length)
  cycles <- as.numeric(sapply(strsplit(names(reads), "_"), function(x) x[3]))
  depth_table[i, cycles] <- tmp_depth
  
  print(i)
}
write.table(depth_table, file="raw/depth.txt")
