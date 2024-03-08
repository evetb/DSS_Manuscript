# PERMANOVA analysis of microbial community beta-diversity
# Accounting for repeated measures and multiple comparisons

# Libraries ---------------------------------------------

library(phyloseq)
library(tidyverse)
library(qiime2R)
library(pairwiseAdonis)

# Data import & manipulation ----------------------------

# Using qiime2R to create phyloseq objects from my experimental data

dss2_physeq <- qza_to_phyloseq(features = "nnn-filtered-table.qza",
                               tree = "rooted-tree.qza",
                               taxonomy = "taxonomy.qza",
                               metadata = "DSS2_metadata.tsv")

dss2_meta <- read.table("DSS2_metadata.tsv", header = TRUE)

# Subsetting my dataset to only include the EdU+ bacterial communities
Epos <- subset_samples(dss2_physeq, sf == 'EdU+')

# rarefying my data to ensure even, random sampling
Epos_rare <- rarefy_even_depth(Epos, rngseed=1, sample.size=0.95*min(sample_sums(Epos)), replace=F)

# Performing PERMANOVA -----------------------------------

# Accounting for repeated measures from: 
# https://thebiobucket.blogspot.com/2011/04/repeat-measure-adonis-lately-i-had-to.html

# Taking my taxa counts and phylogenetic tree from a phyloseq object, Epos_rare
# and using these to compute Bray-Curtis distances
counts <- otu_table(Epos_rare)
tree <- phy_tree(Epos_rare)
rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)

### computing the true R2-value:
print(fit <- adonis2(rbiom_bc ~ Epos_rare@sam_data[["day"]], permutations=1))
print(fit <- adonis2(rbiom_bc ~ Epos_rare@sam_data[["hs"]], permutations=1))

### number of permutations
B <- 1999

### set up a "permControl" object:
### we turn off mirroring as time should only flow in one direction
ctrl <- how(blocks = Epos_rare@sam_data[["cage"]], within = Within(type = "series", mirror = FALSE))

### Number of observations:
nobs <- nrow(Epos_rare@sam_data)

### check permutation (...rows represent the sample id):
### ..they are ok!
### within in each repeated sample (= sites) timepoints are shuffled,
### with keeping the sequence intact (e.g., for site 1: 1,2,3 - 2,3,1 - 3,2,1)
shuffle(nobs, control = ctrl)

### loop:
### in adonis(...) you need to put permutations = 1, otherwise
### adonis will not run

# comparing between days
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand.day <- adonis2(rbiom_bc ~ Epos_rare@sam_data[["day"]][idx],permutations = 1)
}

# comparing between health states (hs)
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand.hs <- adonis2(rbiom_bc ~ Epos_rare@sam_data[["hs"]][idx],permutations = 1)
}

# doing the above (repeated measures) 
# in combination with multiple comparisons
# using pairwise.adonis()
# this function is comparing between 
# the different health states ("hs")
# controlling for cage of origin

PermaMC_hs_rm <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  nobs <- nrow(phylo_exp@sam_data)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs <- pairwise.adonis(rbiom_bc, 
                                   phylo_exp@sam_data[["hs"]][idx], p.adjust.m = 'BH') 
  }
  return(fit.rand.hs)
}

PermaMC_hs_rm(Epos_rare)
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1  Colitis vs Pre-DSS  1 0.1960076 0.6834888 0.02561467   0.805      0.974    
# 2 Colitis vs Recovery  1 0.1476282 0.4885365 0.01844332   0.974      0.974    
# 3      Colitis vs DSS  1 0.2101471 0.7363940 0.02754276   0.743      0.974    
# 4 Pre-DSS vs Recovery  1 0.2159282 0.7099962 0.03126360   0.816      0.974    
# 5      Pre-DSS vs DSS  1 0.1714879 0.6032831 0.02669006   0.910      0.974    
# 6     Recovery vs DSS  1 0.2190225 0.7241156 0.03186551   0.821      0.974

PermaMC_hs_rm(W_rare)
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1      DSS vs Colitis  1 0.2031065 0.9372580 0.03613559   0.467     0.6936    
# 2      DSS vs Pre-DSS  1 0.1425037 0.6804292 0.03000072   0.834     0.8340    
# 3     DSS vs Recovery  1 0.3308622 1.6784346 0.07088453   0.075     0.2000    
# 4  Colitis vs Pre-DSS  1 0.1786481 0.8583200 0.03319318   0.578     0.6936    
# 5 Colitis vs Recovery  1 0.3187504 1.6155042 0.06069786   0.090     0.2000    
# 6 Pre-DSS vs Recovery  1 0.2999243 1.6005235 0.06781729   0.100     0.2000 