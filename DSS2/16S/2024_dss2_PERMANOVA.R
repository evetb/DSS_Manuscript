# PERMANOVA analyses on beta-diversity matrices from DSS #2

# Libraries ------------------------------------------------------------------------------

library(phyloseq)
library(tidyverse)
library(pairwiseAdonis)
#library(RColorBrewer)
#library(data.table)
library(vegan) #for alpha diversity, betadisp
#library(cowplot)
#library(gridExtra)
#library(ape)
#library(biomformat)
library(qiime2R)
#library(rbiom)
#library(stats)
#library(nlme)
#library(compositions)
library(writexl)

# QIIME2R ---------------------------------------------------------------------------------------------

# Using qiime2R to create phyloseq objects from experimental data

#### DSS2 ####

dss2_physeq <- qza_to_phyloseq(features = "nnn-filtered-table.qza",
                               tree = "rooted-tree.qza",
                               taxonomy = "taxonomy.qza",
                               metadata = "DSS2_metadata2.tsv")

dss2_meta <- read.table("DSS2_metadata2.tsv", header = TRUE)

setwd("E:/DESKTOP-UFR1DAD/Documents/McGill/DSS_Manuscript/DSS2/16S/diversity")

# DSS #2 experiment ----------------------------------------------------------------------------------

# checking number of reads per sample
# using sample_sums() on phyloseq object
# sorting small to large using sort
# https://www.nicholas-ollberding.com/post/introduction-to-phyloseq/

sort(sample_sums(dss2_physeq))

# shows that M2W8 only has 2 reads

# quick histogram plot of samples

hist(sample_sums(dss2_physeq), main="Histogram: Read Counts", xlab="Total Reads", 
     border="black", col="#228833", breaks=20)

#### Subsetting ####

# subsetting DSS #2 by cage

F1 <- subset_samples(dss2_physeq, cage == 'F1')
F2 <- subset_samples(dss2_physeq, cage == 'F2')
M1 <- subset_samples(dss2_physeq, cage == 'M1')
M2 <- subset_samples(dss2_physeq, cage == 'M2')

# Subsetting DSS #2 by sorted fraction

Epos <- subset_samples(dss2_physeq, sf == 'EdU+')
Eneg <- subset_samples(dss2_physeq, sf == 'EdU-')
Whole <- subset_samples(dss2_physeq, sf == 'Whole')

# plotting again for sample distribution

hist(sample_sums(Epos), main="Histogram: Read Counts - Epos", xlab="Total Reads", 
     border="black", col="#228833", breaks=20)
hist(sample_sums(Eneg), main="Histogram: Read Counts - Eneg", xlab="Total Reads", 
     border="black", col="#228833", breaks=20)
hist(sample_sums(Whole), main="Histogram: Read Counts - Whole", xlab="Total Reads", 
     border="black", col="#228833", breaks=20)

# you can see that the Whole samples tend to have a lot less reads
# max Epos/Eneg around 50,000 vs max Whole around 35,000
# though the samples >30,000 for Epos & Eneg are outliers

# excluding sample M2W8 since it has sample size = 2
Whole <- subset_samples(dss2_physeq, sf == 'Whole' & samples != 'M2W8') 

# plotting again
hist(sample_sums(Whole), main="Histogram: Read Counts - Whole", xlab="Total Reads", 
     border="black", col="#228833", breaks=20)

# Subsetting sorted fraction by sex

F.Epos <- subset_samples(Epos, sex == 'Female')
M.Epos <- subset_samples(Epos, sex == 'Male')

F.Whole <- subset_samples(Whole, sex == 'Female')
M.Whole <- subset_samples(Whole, sex == 'Male')

# Subsetting DSS #2 by health state

# PreDSS <- subset_samples(dss2_physeq, hs == 'Pre-DSS')
# DSS <- subset_samples(dss2_physeq, hs == 'DSS')
# Recovery <- subset_samples(dss2_physeq, hs == 'Recovery')
# FollowUp <- subset_samples(dss2_physeq, hs == 'Follow-up')

# subsetting each cage into its respective sorted fractions

# F1Epos <- subset_samples(F1, sf == 'EdU+')
# F1Eneg <- subset_samples(F1, sf == 'EdU-')
# F1W <- subset_samples(F1, sf == 'Whole')
# 
# F2Epos <- subset_samples(F2, sf == 'EdU+')
# F2Eneg <- subset_samples(F2, sf == 'EdU-')
# F2W <- subset_samples(F2, sf == 'Whole')
# 
# M1Epos <- subset_samples(M1, sf == 'EdU+')
# M1Eneg <- subset_samples(M1, sf == 'EdU-')
# M1W <- subset_samples(M1, sf == 'Whole')
# 
# M2Epos <- subset_samples(M2, sf == 'EdU+')
# M2Eneg <- subset_samples(M2, sf == 'EdU-')
# M2W <- subset_samples(M2, sf == 'Whole')

# Rarefaction curves --------------------------------------------------------------------------------

# from vegan
# as per https://micca.readthedocs.io/en/latest/phyloseq.html

# new format, have to add as() to convert otu object to matrix
# https://github.com/joey711/phyloseq/issues/1641

Epos_otu_mat <- as(t(otu_table(Epos)), "matrix")
Whole_otu_mat <- as(t(otu_table(Whole)), "matrix")

rarecurve(Epos_otu_mat, step = 50, cex = 0.5) # 

rarecurve(Whole_otu_mat, step = 50, cex = 0.5) # 

# PERMDISP ----------------------------------------------------------------------------------------

# calculating dispersion & testing for significance
# using vegan's beta_disper() function
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/
# https://rfunctions.blogspot.com/2019/03/betadisper-and-adonis-homogeneity-of.html

## functions ----------------------------------------------------------
# Key:
## RMMC = repeated measures, multiple comparisons
## hs = health state
## bc = Bray-Curtis
## wu = weighted Unifrac
## uu = unweighted Unifrac
## jc = Jaccard

Permadisp_RMMC_hs_bc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["hs"]], 
                       #bias.adjust = TRUE,
                       type = "centroid")
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(dist(phy_bd$distances),
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

Permadisp_RMMC_hs_wu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = TRUE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["hs"]], 
                       #bias.adjust = TRUE,
                       type = "centroid")
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(dist(phy_bd$distances),
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

Permadisp_RMMC_hs_uu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = FALSE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["hs"]], 
                       #bias.adjust = TRUE,
                       type = "centroid")
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(dist(phy_bd$distances),
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

Permadisp_RMMC_hs_jc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'jaccard', weighted = FALSE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["hs"]], 
                       #bias.adjust = TRUE,
                       type = "centroid")
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(dist(phy_bd$distances),
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

## exporting to excel --------------------------------------------------

### Epos ----------------------------------------------------------------

write_xlsx(as.data.frame(Permadisp_RMMC_hs_bc(Epos)), "Epos_Permadisp_hs_bc.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_wu(Epos)), "Epos_Permadisp_hs_wu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_uu(Epos)), "Epos_Permadisp_hs_uu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_jc(Epos)), "Epos_Permadisp_hs_jc.xlsx")

### Whole ----------------------------------------------------------------

write_xlsx(as.data.frame(Permadisp_RMMC_hs_bc(Whole)), "Whole_Permadisp_hs_bc.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_wu(Whole)), "Whole_Permadisp_hs_wu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_uu(Whole)), "Whole_Permadisp_hs_uu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_jc(Whole)), "Whole_Permadisp_hs_jc.xlsx")

# PERMANOVA ------------------------------------------------------------------------------------------

# Creating a PERMANOVA function (code adapted from Taguer, 2021)

# To rarefy your counts 
# (equivalent to p-sampling depth in Qiime2 - 
# qiime diversity core-metrics-phylogenetic),
# Add the following line for each phyloseq object, 
# which rarefies all samples 
# to a depth of 95% of the read count 
# of the sample with the lowest number of reads: 

# rarefy_even_depth(physeq_object, rngseed=1, sample.size=0.95*min(sample_sums(expt1)), replace=F)

# There's an explanation in this tutorial here: 
# https://micca.readthedocs.io/en/latest/phyloseq.html#import-data-and-preparation

## DSS #2 PERMANOVA Functions ----------------------------------------------------------------

### On EdU+ fraction ------------------------------------------------------------------------

#### Rarefaction ---------------------------------------------------------------------
Epos_rare <- rarefy_even_depth(Epos, rngseed=1, sample.size=0.95*min(sample_sums(Epos)), replace=F)

#### PERMANOVA significance tests --------------------------------------------------

# creating Permanova function, 
# exporting the table of results (data frame)
# only tells global significance 
# and not post-hoc multiple comparisons

# here doing Bray-Curtis

Permanova_hs <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  adonis2(rbiom_bc ~ hs, data = as(sample_data(phylo_exp), "data.frame"))
}

Permanova_hs(Epos)

## post-hoc multiple comparisons using pairwise.adonis2()
## from https://github.com/pmartinezarbizu/pairwiseAdonis

PermaMC_hs_bc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["hs"]]), p.adjust.m = 'BH') 
}

PermaMC_hs_bc(Epos)

##### per day ---------------------------------------------

PermaMC_day_bc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]]), p.adjust.m = 'BH') 
}

PermaMC_day_wu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = TRUE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]]), p.adjust.m = 'BH') 
}

PermaMC_day_uu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = FALSE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]]), p.adjust.m = 'BH') 
}

PermaMC_day_jc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'jaccard', weighted = FALSE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]]), p.adjust.m = 'BH') 
}

# NOTE: set.seed() guarantees repeatable results
# otherwise your p-values change each time you run the function
# https://stackoverflow.com/questions/56938448/why-does-adonis-from-vegan-returns-a-different-p-value-every-time-it-is

# exporting to excel

###### Epos ----------------------------------------------------

write_xlsx(as.data.frame(PermaMC_day_bc(Epos)), "Epos_PermaMC_day_bc.xlsx")
write_xlsx(as.data.frame(PermaMC_day_wu(Epos)), "Epos_PermaMC_day_wu.xlsx")
write_xlsx(as.data.frame(PermaMC_day_uu(Epos)), "Epos_PermaMC_day_uu.xlsx")
write_xlsx(as.data.frame(PermaMC_day_jc(Epos)), "Epos_PermaMC_day_jc.xlsx")

###### Whole ----------------------------------------------------

write_xlsx(as.data.frame(PermaMC_day_bc(Whole)), "Whole_PermaMC_day_bc.xlsx")
write_xlsx(as.data.frame(PermaMC_day_wu(Whole)), "Whole_PermaMC_day_wu.xlsx")
write_xlsx(as.data.frame(PermaMC_day_uu(Whole)), "Whole_PermaMC_day_uu.xlsx")
write_xlsx(as.data.frame(PermaMC_day_jc(Whole)), "Whole_PermaMC_day_jc.xlsx")

##### Accounting for Repeated Measures ------------------
# from: 
# https://thebiobucket.blogspot.com/2011/04/repeat-measure-adonis-lately-i-had-to.html

###### Functions ----------------------------------------------

# repeated measures in combination with pairwise comparisons

## Key:
## MC = multiple comparisons
## hs = health state
## bc = Bray-Curtis
## wu = weighted Unifrac
## uu = unweighted Unifrac
## jc = Jaccard

PermaMC_hs_bc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH') 
  }
  return(fit.rand.hs.mc)
}

PermaMC_hs_wu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = TRUE, tree = tree)
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH') 
  }
  return(fit.rand.hs.mc)
}

PermaMC_hs_uu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = FALSE, tree = tree)
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH') 
  }
  return(fit.rand.hs.mc)
}

PermaMC_hs_jc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'jaccard', weighted = FALSE, tree = tree)
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["hs"]][idx], 
                                      p.adjust.m = 'BH') 
  }
  return(fit.rand.hs.mc)
}

###### Calculations -----------------------------------------

PermaMC_hs_bc_Epos <- PermaMC_hs_bc(Epos) 
PermaMC_hs_bc_Whole <- PermaMC_hs_bc(Whole) 

PermaMC_hs_wu_Epos <- PermaMC_hs_wu(Epos)
PermaMC_hs_wu_Whole <- PermaMC_hs_wu(Whole)

PermaMC_hs_uu_Epos <- PermaMC_hs_uu(Epos)
PermaMC_hs_uu_Whole <- PermaMC_hs_uu(Whole)

PermaMC_hs_jc_Epos <- PermaMC_hs_jc(Epos)
PermaMC_hs_jc_Whole <- PermaMC_hs_jc(Whole)

# exporting to excel

write_xlsx(PermaMC_hs_bc_Epos, "PermaMC_hs_bc_Epos.xlsx")
write_xlsx(PermaMC_hs_bc_Whole, "PermaMC_hs_bc_Whole.xlsx")

write_xlsx(PermaMC_hs_wu_Epos, "PermaMC_hs_wu_Epos.xlsx")
write_xlsx(PermaMC_hs_wu_Whole, "PermaMC_hs_wu_Whole.xlsx")

write_xlsx(PermaMC_hs_uu_Epos, "PermaMC_hs_uu_Epos.xlsx")
write_xlsx(PermaMC_hs_uu_Whole, "PermaMC_hs_uu_Whole.xlsx")

write_xlsx(PermaMC_hs_jc_Epos, "PermaMC_hs_jc_Epos.xlsx")
write_xlsx(PermaMC_hs_jc_Whole, "PermaMC_hs_jc_Whole.xlsx")

####### Per sorted fraction (sf) & sex ------------------------------------------

# subsetting (as written above)

F.Epos <- subset_samples(Epos, sex == 'Female')
M.Epos <- subset_samples(Epos, sex == 'Male')

F.Whole <- subset_samples(Whole, sex == 'Female')
M.Whole <- subset_samples(Whole, sex == 'Male')

# calculations

## Females
PermaMC_hs_bc_FEpos <- PermaMC_hs_bc(F.Epos) 
PermaMC_hs_bc_FWhole <- PermaMC_hs_bc(F.Whole) 

PermaMC_hs_wu_FEpos <- PermaMC_hs_wu(F.Epos)
PermaMC_hs_wu_FWhole <- PermaMC_hs_wu(F.Whole)

PermaMC_hs_uu_FEpos <- PermaMC_hs_uu(F.Epos)
PermaMC_hs_uu_FWhole <- PermaMC_hs_uu(F.Whole)

PermaMC_hs_jc_FEpos <- PermaMC_hs_jc(F.Epos)
PermaMC_hs_jc_FWhole <- PermaMC_hs_jc(F.Whole)

## Males

PermaMC_hs_bc_MEpos <- PermaMC_hs_bc(M.Epos) 
PermaMC_hs_bc_MWhole <- PermaMC_hs_bc(M.Whole) 

PermaMC_hs_wu_MEpos <- PermaMC_hs_wu(M.Epos)
PermaMC_hs_wu_MWhole <- PermaMC_hs_wu(M.Whole)

PermaMC_hs_uu_MEpos <- PermaMC_hs_uu(M.Epos)
PermaMC_hs_uu_MWhole <- PermaMC_hs_uu(M.Whole)

PermaMC_hs_jc_MEpos <- PermaMC_hs_jc(M.Epos)
PermaMC_hs_jc_MWhole <- PermaMC_hs_jc(M.Whole)

# exporting to excel

## Females
write_xlsx(PermaMC_hs_bc_FEpos, "PermaMC_hs_bc_FEpos.xlsx")
write_xlsx(PermaMC_hs_bc_FWhole, "PermaMC_hs_bc_FWhole.xlsx")

write_xlsx(PermaMC_hs_wu_FEpos, "PermaMC_hs_wu_FEpos.xlsx")
write_xlsx(PermaMC_hs_wu_FWhole, "PermaMC_hs_wu_FWhole.xlsx")

write_xlsx(PermaMC_hs_uu_FEpos, "PermaMC_hs_uu_FEpos.xlsx")
write_xlsx(PermaMC_hs_uu_FWhole, "PermaMC_hs_uu_FWhole.xlsx")

write_xlsx(PermaMC_hs_jc_FEpos, "PermaMC_hs_jc_FEpos.xlsx")
write_xlsx(PermaMC_hs_jc_FWhole, "PermaMC_hs_jc_FWhole.xlsx")

## Males

write_xlsx(PermaMC_hs_bc_MEpos, "PermaMC_hs_bc_MEpos.xlsx")
write_xlsx(PermaMC_hs_bc_MWhole, "PermaMC_hs_bc_MWhole.xlsx")

write_xlsx(PermaMC_hs_wu_MEpos, "PermaMC_hs_wu_MEpos.xlsx")
write_xlsx(PermaMC_hs_wu_MWhole, "PermaMC_hs_wu_MWhole.xlsx")

write_xlsx(PermaMC_hs_uu_MEpos, "PermaMC_hs_uu_MEpos.xlsx")
write_xlsx(PermaMC_hs_uu_MWhole, "PermaMC_hs_uu_MWhole.xlsx")

write_xlsx(PermaMC_hs_jc_MEpos, "PermaMC_hs_jc_MEpos.xlsx")
write_xlsx(PermaMC_hs_jc_MWhole, "PermaMC_hs_jc_MWhole.xlsx")

#### Ordination (PCoA) -------------------------------------------------------------

# see the following for plot_ordination:
# https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_ordination

# creating PCoA plot from the Epos phyloseq object
# (and using forcats to fct_relevel the health states)

##### Per sorted fraction -----------------------------------

## Bray-Curtis

PCoA_hs_epos_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  #print(rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ BC") + 
    geom_point(aes(shape = sex), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=24)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
    #stat_ellipse(aes(fill = hs) geom = "polygon", alpha = 0.15)
}


## Weighted Unifrac

PCoA_hs_epos_wu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ WU") + 
    geom_point(aes(shape = sex), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
}

## Unweighted Unifrac

PCoA_hs_epos_uu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ UU") + 
    geom_point(aes(shape = sex), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
}

## Jaccard

PCoA_hs_epos_jc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ cells") + 
    geom_point((aes(shape = sex)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    # theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
    #       axis.text = element_text(size=20),
    #       legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    theme(plot.title = element_text(hjust=0.5, size = 45), axis.title=element_text(size=33), 
          axis.text.x = element_text(size=24),axis.text.y=element_text(size=24), 
          legend.text=element_text(size=33), legend.title=element_text(size=33)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
}


PCoA_hs_epos_bc(Epos) 
PCoA_hs_epos_wu(Epos) 
PCoA_hs_epos_uu(Epos)
PCoA_hs_epos_jc(Epos)

# documentation: https://rdrr.io/github/cmmr/rbiom/man/beta.div.html

### On Whole community ------------------------------------------------------------------------

#### Rarefaction ---------------------------------------------------------------------------
# W_rare <- rarefy_even_depth(Whole, rngseed=1, sample.size=0.95*min(sample_sums(Whole)), replace=F)

#### PERMANOVA significance tests ------------------------------------------------------

# globally on health states (hs)
Perm_W_rare <- Permanova_hs(W_rare)
Perm_W_rare

# w/post-hoc multiple comparisons btwn hs

PermaMC_hs(W_rare)
W_PermaMC_hs.df <- as.data.frame(PermaMC_hs(W_rare))

write_csv(W_PermaMC_hs.df, "Whole_PERMANOVA_hs.csv")

# w/post-hoc multiple comparisons btwn days

PermaMC_day(W_rare)
W_PermaMC_day.df <- as.data.frame(PermaMC_day(W_rare))

#### Ordination (PCoA) --------------------------------------------------------------------

# Bray-Curtis

PCoA_hs_w_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole BC") + 
    geom_point(aes(shape = sex), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
}

## Weighted Unifrac

PCoA_hs_w_wu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole WU") + 
    geom_point(aes(shape = sex), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
}

## Unweighted Unifrac

PCoA_hs_w_uu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole UU") + 
    geom_point(aes(shape = sex), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
}

## Jaccard 

PCoA_hs_w_jc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole community") + 
    geom_point((aes(shape = sex)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 45), axis.title=element_text(size=33), 
          axis.text.x = element_text(size=24),axis.text.y=element_text(size=24), 
          legend.text=element_text(size=33), legend.title=element_text(size=33)) +
    # theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
    #       axis.text = element_text(size=20),
    #       legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State", shape = "Sex") +
    stat_ellipse(linewidth = 2)
}

PCoA_hs_w_bc(Whole) 
PCoA_hs_w_wu(Whole) 
PCoA_hs_w_uu(Whole) 
PCoA_hs_w_jc(Whole) 

## Per sorted fraction & sex ---------------------------------------------

F.Epos
F.Whole
M.Epos
M.Whole

### Ordination ------------------------------------------------------------

#### Epos ---------------------------------------------------------------

##### Female -------------------------------------------------------------

## Bray-Curtis

PCoA_epos_f_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  #print(rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Female BC") + 
    stat_ellipse(linewidth = 2) +
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") 
}

## Weighted Unifrac

PCoA_epos_f_wu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Female WU") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Unweighted Unifrac

PCoA_epos_f_uu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Female UU") + 
    geom_point(aes(shape = cage), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Jaccard

PCoA_epos_f_jc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Female JC") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

PCoA_epos_f_bc(F.Epos)
PCoA_epos_f_wu(F.Epos)
PCoA_epos_f_uu(F.Epos)
PCoA_epos_f_jc(F.Epos)

##### Male ----------------------------------------------------------------

## Bray-Curtis

PCoA_epos_m_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  #print(rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Male BC") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
  #stat_ellipse(aes(fill = hs) geom = "polygon", alpha = 0.15)
}

## Weighted Unifrac

PCoA_epos_m_wu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Male WU") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Unweighted Unifrac

PCoA_epos_m_uu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Male UU") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Jaccard

PCoA_epos_m_jc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "EdU+ Male JC") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20),
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

PCoA_epos_m_bc(M.Epos)
PCoA_epos_m_wu(M.Epos)
PCoA_epos_m_uu(M.Epos)
PCoA_epos_m_jc(M.Epos)

#### Whole --------------------------------------------------------------

##### Female -----------------------------------------------------

## Bray-Curtis

PCoA_w_f_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Female BC") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Weighted Unifrac

PCoA_w_f_wu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Female WU") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Unweighted Unifrac

PCoA_w_f_uu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Female UU") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20),
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Jaccard

PCoA_w_f_jc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Female JC") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20),
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

PCoA_w_f_bc(F.Whole)
PCoA_w_f_wu(F.Whole)
PCoA_w_f_uu(F.Whole)
PCoA_w_f_jc(F.Whole)

##### Male -----------------------------------------------------

## Bray-Curtis

PCoA_w_m_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Male BC") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20), 
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Weighted Unifrac

PCoA_w_m_wu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Male WU") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20),
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Unweighted Unifrac

PCoA_w_m_uu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Male UU") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20),
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

## Jaccard

PCoA_w_m_jc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["hs"]] <- fct_relevel(phylo_exp@sam_data[["hs"]], 
                                            c("Baseline", "Pre-symptomatic", "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, 
                  color= "hs", # omitting this gets rid of the ellipsoid lines
                  title = "Whole Male JC") + 
    geom_point((aes(shape = cage)), size=7) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    theme_bw() +  #theme(text = element_text(size = 18)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    theme(plot.title = element_text(hjust=0.5, size = 36), axis.title=element_text(size=26), 
          axis.text = element_text(size=20),
          legend.text=element_text(size=20), legend.title=element_text(size=26)) +
    labs(color = "Health State",
         shape = "Cage") +
    stat_ellipse(linewidth = 2)
}

PCoA_w_m_bc(M.Whole)
PCoA_w_m_wu(M.Whole)
PCoA_w_m_uu(M.Whole)
PCoA_w_m_jc(M.Whole)