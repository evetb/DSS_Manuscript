# running gRodon on the daily assemblies 
# of the DSS #2 mouse experiment

#### Libraries #####

library(gRodon)
library(Biostrings)
library(tictoc)

#### Getting filename variables ####

F2_ffns <- list.files(path = "")
F2_ffn_files <- grep(".*fasta", F2_ffns, value=T)

F2_cds <- list.files(path = "")
F2_cds_files <- grep("CDS_names+", F2_cds, value=T)

M2_ffns <- list.files(path = "")
M2_ffn_files <- grep(".*fasta", M2_ffns, value=T)

M2_cds <- list.files(path = "")
M2_cds_files <- grep("CDS_names+", M2_cds, value=T)

#### gRodon on all .fasta files, per cage ####

###### MMv2 ####

setwd("")

# Using a for loop to properly get 
# each file based on its name,
# as well as running through gRodon vignette steps
# using metagenome mode version 2 (MMv2)
# and properly labelling the output file
# https://microbialgamut.com/gRodon-vignette

for (i in F2_ffn_files){
  BIN <- strsplit(i, "\\.")[[1]][1:2]
  BIN <- BIN[1]
  genes <- readDNAStringSet(paste("F2_daily_ffn/", i, sep=""))
  CDS_IDs <- readLines((paste("F2_daily_gff/CDS_names_",BIN,".txt",sep="")))
  gene_IDs <- gsub(" .*","",names(genes))
  genes <- genes[gene_IDs %in% CDS_IDs]
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
  read_depths <- read.delim(paste("doc/", BIN, ".tsv", sep=""),
                            stringsAsFactors = FALSE)
  depths <- read_depths$meandepth
  names(depths) <- read_depths$X.rname
  depth_of_coverage <- depths[gsub(" .*", "", names(genes))]
  capture.output(predictGrowth(genes, highly_expressed, mode = "metagenome_v2", 
                               depth_of_coverage = depth_of_coverage, n_le = 1000), 
                 file=((paste("out/v2/",BIN,".tsv",sep=""))), append=FALSE, type=c("output"), split=TRUE)
}

setwd("")

for (i in M2_ffn_files){
  BIN <- strsplit(i, "\\.")[[1]][1:2]
  BIN <- BIN[1]
  genes <- readDNAStringSet(paste("M2_daily_ffn/", i, sep=""))
  CDS_IDs <- readLines((paste("M2_daily_gff/CDS_names_",BIN,".txt",sep="")))
  gene_IDs <- gsub(" .*","",names(genes))
  genes <- genes[gene_IDs %in% CDS_IDs]
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
  read_depths <- read.delim(paste("doc/", BIN, ".tsv", sep=""),
                            stringsAsFactors = FALSE)
  depths <- read_depths$meandepth
  names(depths) <- read_depths$X.rname
  depth_of_coverage <- depths[gsub(" .*", "", names(genes))]
  capture.output(predictGrowth(genes, highly_expressed, mode = "metagenome_v2", 
                               depth_of_coverage = depth_of_coverage, n_le = 1000), 
                 file=((paste("out/v2/", BIN,".tsv",sep=""))), append=FALSE, type=c("output"), split=TRUE)
}