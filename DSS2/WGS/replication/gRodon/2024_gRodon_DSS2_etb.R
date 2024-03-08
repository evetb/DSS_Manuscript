#running gRodon on each bin from DSS#2

#### Libraries ####

library(gRodon)
library(Biostrings)

getwd() # Get working directory
setwd("") # Set working directory

#### Getting filename variables ####

F2_files <- list.files(path = "/F2")
F2_ffn <- grep(".*fasta", F2_files, value=T)

F2_cds <- list.files(path = "/F2_gff")
F2_cds_files <- grep("CDS_names+", F2_cds, value=T)

M2_files <- list.files(path = "/M2")
M2_ffn <- grep(".*fasta", M2_files, value=T)

M2_cds <- list.files(path = "/M2_gff")
M2_cds_files <- grep("CDS_names+", M2_cds, value=T)

#### gRodon on all .fasta files, per cage ####

# strsplit uses regex;
# since the period character is used 
# as a wildcard in regex, we need to
# escape from it to split on period, 
# using backslash \ 
# then we have to escape the \ as well, 
# so we end up adding two \ and then 
# the period, as seen below

# using paste, I am combining the first two parts of 
# the bin name back together

for (i in F2_ffn){
  print(i)
  BIN <- strsplit(i, "\\.")[[1]][1:2]
  BIN <- paste(BIN[1], BIN[2], sep = ".")
  print(BIN)
}

# Using a for loop to properly get 
# each file based on its name,
# as well as running through gRodon vignette steps
# using partial mode
# and properly labelling the output file
# https://microbialgamut.com/gRodon-vignette

setwd("/F2")

for (i in F2_ffn){
  BIN <- strsplit(i, "\\.")[[1]][1:2]
  BIN <- paste(BIN[1], BIN[2], sep = ".")
  genes <- readDNAStringSet(i)
  CDS_IDs <- readLines((paste("CDS_names_",BIN,".fa.txt",sep="")))
  gene_IDs <- gsub(" .*","",names(genes))
  genes <- genes[gene_IDs %in% CDS_IDs]
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
  capture.output(predictGrowth(genes, 
                               highly_expressed, mode = "partial"), 
                 file=((paste(BIN,".tsv",sep=""))), append=FALSE, type=c("output"), split=TRUE)
}

setwd("/M2")

for (i in M2_ffn){
  BIN <- strsplit(i, "\\.")[[1]][1:2]
  BIN <- paste(BIN[1], BIN[2], sep = ".")
  genes <- readDNAStringSet(i)
  CDS_IDs <- readLines((paste("CDS_names_",BIN,".fa.txt",sep="")))
  gene_IDs <- gsub(" .*","",names(genes))
  genes <- genes[gene_IDs %in% CDS_IDs]
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
  capture.output(predictGrowth(genes, 
                               highly_expressed, mode = "partial"), 
                 file=((paste(BIN,".tsv",sep=""))), append=FALSE, type=c("output"), split=TRUE)
}
