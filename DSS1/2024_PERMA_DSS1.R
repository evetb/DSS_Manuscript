# PERMANOVA analyses on beta-diversity matrices from DSS #1

# Libraries -----------------------------------------------------------------------------

library(phyloseq)
library(tidyverse)
library(pairwiseAdonis)
library(vegan)
library(cowplot)
library(gridExtra)
library(ape)
library(biomformat)
library(qiime2R)
library(rbiom)
library(stats)
library(nlme)
library(compositions)
library(ggforce)
library(writexl) #for going from R to excel

# QIIME2R ---------------------------------------------------------------------------------------------

# Using qiime2R to create phyloseq objects from experimental data

#### DSS1 ####

dss1_physeq <- qza_to_phyloseq(features = "table.qza",
                               tree = "rooted-tree.qza",
                               taxonomy = "taxonomy.qza",
                               metadata = "DSS1_metadata2.tsv")

dss1_meta <- read.table("DSS1_metadata2.tsv", header = TRUE)

# DSS #1 experiment ----------------------------------------------------------------------------------

#### Subsetting ####

# subsetting DSS #1 by cage

F1 <- subset_samples(dss1_physeq, cage == 'F1')
F2 <- subset_samples(dss1_physeq, cage == 'F2')
M1 <- subset_samples(dss1_physeq, cage == 'M1')
M2 <- subset_samples(dss1_physeq, cage == 'M2')

# subsetting DSS #1 by sex

Females <- subset_samples(dss1_physeq, sex == 'Female')
Males <- subset_samples(dss1_physeq, sex == 'Male')

# Rarefaction curves --------------------------------------------------------------------------------

# from vegan
# as per https://micca.readthedocs.io/en/latest/phyloseq.html

# new format, have to add as() to convert otu object to matrix
# https://github.com/joey711/phyloseq/issues/1641

otu_mat <- as(t(otu_table(dss1_physeq)), "matrix")

rarecurve(otu_mat, step = 50, cex = 0.5)

rarecurve(t(otu_table(F1)), step = 50, cex = 0.5)
rarecurve(t(otu_table(F2)), step = 50, cex = 0.5)
rarecurve(t(otu_table(M1)), step = 50, cex = 0.5)
rarecurve(t(otu_table(M2)), step = 50, cex = 0.5)

# PERMDISP ----------------------------------------------------------------------------------------

# calculating dispersion & testing for significance
# using vegan's beta_disper() function
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/
# https://rfunctions.blogspot.com/2019/03/betadisper-and-adonis-homogeneity-of.html

## functions ----------------------------------------------------------

# RMMC stands for repeated measures, multiple comparisons

Permadisp_RMMC_hs_bc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["health.state"]], 
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
                                      phylo_exp@sam_data[["health.state"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

Permadisp_RMMC_hs_wu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = TRUE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["health.state"]], 
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
                                      phylo_exp@sam_data[["health.state"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

Permadisp_RMMC_hs_uu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = FALSE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["health.state"]], 
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
                                      phylo_exp@sam_data[["health.state"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

Permadisp_RMMC_hs_jc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'jaccard', weighted = FALSE, tree = tree)
  phy_bd <- betadisper(rbiom_bc, group = phylo_exp@sam_data[["health.state"]], 
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
                                      phylo_exp@sam_data[["health.state"]][idx], 
                                      p.adjust.m = 'BH')
  }
  return(fit.rand.hs.mc)
}

## exporting to excel --------------------------------------------------

### All samples ----------------------------------------------------------------

write_xlsx(as.data.frame(Permadisp_RMMC_hs_bc(dss1_physeq)), "dss1_Permadisp_hs_bc.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_wu(dss1_physeq)), "dss1_Permadisp_hs_wu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_uu(dss1_physeq)), "dss1_Permadisp_hs_uu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_jc(dss1_physeq)), "dss1_Permadisp_hs_jc.xlsx")

### Females ----------------------------------------------------------------

write_xlsx(as.data.frame(Permadisp_RMMC_hs_bc(Females)), "Females_Permadisp_hs_bc.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_wu(Females)), "Females_Permadisp_hs_wu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_uu(Females)), "Females_Permadisp_hs_uu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_jc(Females)), "Females_Permadisp_hs_jc.xlsx")

### Males ----------------------------------------------------------------

write_xlsx(as.data.frame(Permadisp_RMMC_hs_bc(Males)), "Males_Permadisp_hs_bc.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_wu(Males)), "Males_Permadisp_hs_wu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_uu(Males)), "Males_Permadisp_hs_uu.xlsx")
write_xlsx(as.data.frame(Permadisp_RMMC_hs_jc(Males)), "Males_Permadisp_hs_jc.xlsx")

# PERMANOVA ------------------------------------------------------------------------------------------

# Creating a PERMANOVA function (code adapted from Taguer, 2021)

# To rarefy your counts (equivalent to p-sampling depth in Qiime2 
# - qiime diversity core-metrics-phylogenetic),
# Add the following line for each phyloseq object, which rarefies 
# all samples to a depth of 95% of the read count 
# of the sample with the lowest number of reads: 

# rarefy_even_depth(phylo, rngseed=1, sample.size=0.95*min(sample_sums(expt1)), replace=F)

# where phylo is the phyloseq object.
# There's an explanation in this tutorial here: 
# https://micca.readthedocs.io/en/latest/phyloseq.html#import-data-and-preparation

## DSS #1 PERMANOVA Functions ----------------------------------------------------------------

### On all cages ------------------------------------------------------------------------

#### Rarefaction ---------------------------------------------------------------------
dss1_rare <- rarefy_even_depth(dss1_physeq, rngseed=1, 
                               sample.size=0.95*min(sample_sums(dss1_physeq)), replace=F)


#### PERMANOVA significance tests --------------------------------------------------

##### Overall PERMANOVA --------------------------------

# creating Permanova function, exporting the table of results (data frame)
# only tells global significance and not post-hoc multiple comparisons

Permanova_hs_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  rbiom_bc = rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  adonis2(rbiom_bc ~ health.state, data = as(sample_data(phylo_exp), "data.frame"))
}

Permanova_hs_bc(dss1_rare)
Permanova_hs_bc(dss1_physeq)

##### Multiple comparisons ----------------------------------------

## post-hoc multiple comparisons using pairwise.adonis2()
## from https://github.com/pmartinezarbizu/pairwiseAdonis

# I input the distance matrix generated by rbiom 
# & specified I wanted comparisons across health states (hs)
# making sure that these are considered as factors & not characters,
# by using as.factor()

# per health state

PermaMC_hs_bc <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  rbiom_bc = rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["health.state"]]))
  }

PermaMC_hs_wu <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  rbiom_bc = rbiom::beta.div(counts, 'unifrac', weighted = TRUE, tree = tree)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["health.state"]]))
}

# per day

PermaMC_day_bc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]])) 
}

PermaMC_day_wu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = TRUE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]])) 
}

PermaMC_day_uu <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'unifrac', weighted = FALSE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]])) 
}

PermaMC_day_jc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'jaccard', weighted = FALSE, tree = tree)
  set.seed(123)
  pairwise.adonis(rbiom_bc, as.factor(phylo_exp@sam_data[["day"]])) 
}

PermaMC_hs_bc(dss1_rare)
PermaMC_hs_wu(dss1_rare)

PermaMC_day_bc(dss1_physeq)
PermaMC_day_wu(dss1_physeq) 
PermaMC_day_uu(dss1_physeq)
PermaMC_day_jc(dss1_physeq)

# saving results of PermaMC as a data frame
dss1_PermaMC_hs_bc.df <- as.data.frame(PermaMC_hs_bc(dss1_rare))
dss1_PermaMC_hs_wu.df <- as.data.frame(PermaMC_hs_wu(dss1_rare))

# saving the data frame as a .csv file
write_csv(dss1_PermaMC_hs_bc.df, "dss1_PERMANOVA_hs_bc.csv")
write_csv(dss1_PermaMC_hs_wu.df, "dss1_PERMANOVA_hs_wu.csv")

# p.adjust is 1 for all (not significant)
write_csv(as.data.frame(PermaMC_day_bc(dss1_physeq)), "dss1_PermaMC_day_bc.csv")
write_csv(as.data.frame(PermaMC_day_wu(dss1_physeq)), "dss1_PermaMC_day_wu.csv")
write_csv(as.data.frame(PermaMC_day_uu(dss1_physeq)), "dss1_PermaMC_day_uu.csv")
write_csv(as.data.frame(PermaMC_day_jc(dss1_physeq)), "dss1_PermaMC_day_jc.csv")

##### Repeated measures -------------------------------------------------------

###### Individual testing -----------------------------------------------------

# from: 
# https://thebiobucket.blogspot.com/2011/04/repeat-measure-adonis-lately-i-had-to.html

# transforming count data to matrix
# with samples as rows
# because i need nobs (number of observations) of counts and ctrl 
# to have equal row numbers to do shuffle
# https://github.com/joey711/phyloseq/issues/613

counts <- as(otu_table(dss1_physeq), "matrix")
if(taxa_are_rows(dss1_physeq)){counts <- t(counts)}

tree <- phy_tree(dss1_physeq)
rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)

### number of permutations
B <- 1999

### set up a "permControl" object:
### we turn off mirroring as time should only flow in one direction
### (note that permControl was deprecated and now you use how() instead)
### https://fromthebottomoftheheap.net/2013/12/17/permute-0.8-0-released/

ctrl <- how(blocks = dss1_physeq@sam_data[["cage"]], 
                    within = Within(type = "series", mirror = FALSE))

### Number of observations:
nobs <- nrow(counts) #52 rows 

### check permutation (...rows represent the sample id)
### within each repeated sample (= cage), time points are shuffled,
### with keeping the sequence intact (e.g., for site 1: 1,2,3 - 2,3,1 - 3,2,1)
shuffle(nobs, control = ctrl)

### loop:
### in adonis(...) you need to put permutations = 1, 
### otherwise adonis will not run

# day - nothing significant
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand.day <- adonis2(rbiom_bc ~ dss1_physeq@sam_data[["day"]][idx],permutations = 1)
}

# health state (hs) - nothing significant
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand.hs <- adonis2(rbiom_bc ~ dss1_physeq@sam_data[["health.state"]][idx],permutations = 1)
}

###### Functions --------------------------------------------

# doing the above in combination with pairwise comparisons
# using pairwise.adonis2 w/new syntax (see ANOVA_DISP.R) doesn't work with [idx]
# it should be accounting for cage anyway i think w/the ctrl variable
# so i can use the normal pairwise.adonis

PermaMC_hs_bc <- function(phylo_exp) {
  counts <- otu_table(phylo_exp)
  tree <- phy_tree(phylo_exp)
  rbiom_bc <- rbiom::beta.div(counts, 'bray-curtis', weighted = TRUE, tree = tree)
  B <- 1999
  ctrl <- how(blocks = phylo_exp@sam_data[["cage"]], 
              within = Within(type = "series", mirror = FALSE))
  counts <- as(otu_table(phylo_exp), "matrix")
  if(taxa_are_rows(phylo_exp)){counts <- t(counts)}
  # wondering if the above counts stuff is even necessary
  # since the below is just counting rows
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["health.state"]][idx], 
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
  # wondering if the above counts stuff is even necessary
  # since the below is just counting rows
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["health.state"]][idx], 
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
  # wondering if the above counts stuff is even necessary
  # since the below is just counting rows
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["health.state"]][idx], 
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
  # wondering if the above counts stuff is even necessary
  # since the below is just counting rows
  nobs <- nrow(counts)
  set.seed(123)
  for(i in 2:(B+1)){
    idx <- shuffle(nobs, control = ctrl)
    fit.rand.hs.mc <- pairwise.adonis(rbiom_bc,
                                      phylo_exp@sam_data[["health.state"]][idx], 
                                      p.adjust.m = 'BH') 
  }
  return(fit.rand.hs.mc)
}

###### Calculations ------------------------------------------------------------

dss1_PermaMC_hs_bc.df <- PermaMC_hs_bc(dss1_physeq)
dss1_PermaMC_hs_wu.df <- PermaMC_hs_wu(dss1_physeq)
dss1_PermaMC_hs_uu.df <- PermaMC_hs_uu(dss1_physeq)
dss1_PermaMC_hs_jc.df <- PermaMC_hs_jc(dss1_physeq)

# exporting to excel

write_xlsx(dss1_PermaMC_hs_bc.df, "dss1_PermaMC_hs_bc.df.xlsx")
write_xlsx(dss1_PermaMC_hs_wu.df, "dss1_PermaMC_hs_wu.df.xlsx")
write_xlsx(dss1_PermaMC_hs_uu.df, "dss1_PermaMC_hs_uu.df.xlsx")
write_xlsx(dss1_PermaMC_hs_jc.df, "dss1_PermaMC_hs_jc.df.xlsx")

#### Ordination (PCoA) -------------------------------------------------------------

# see the following for plot_ordination:
# https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_ordination
# also https://cran.r-project.org/web/packages/rbiom/rbiom.pdf

# creating PCoA plot from the Epos phyloseq object
# (and using forcats in tidyverse to fct_relevel the health states)

# documentation: https://rdrr.io/github/cmmr/rbiom/man/beta.div.html

# numbers to specify shapes in scale_shape_manual from:
# https://www.datanovia.com/en/blog/ggplot-point-shapes-best-tips/

PCoA_hs_bc_all <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                            c("Baseline", "Pre-symptomatic", 
                                              "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state",
                  title = "Bray-Curtis") + 
    geom_point((aes(shape = sex)), size=5) + theme_bw() + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    # theme(text = element_text(size = 12)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    scale_shape_manual(values=c(16,18)) +
    labs(color = "Health State", shape = "Sex") + 
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_wu_all <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_wu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_wu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", 
                  title = "Weighted Unifrac") + 
    geom_point((aes(shape = sex)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    #scale_shape_manual(values=c(19,5,20,1)) +
    labs(color = "Health State", shape = "Sex") + 
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_uu_all <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_uu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_uu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state",
                  title = "Unweighted Unifrac") + 
    geom_point((aes(shape = sex)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) +
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    scale_shape_manual(values=c(16,18)) +
    labs(color = "Health State", shape = "Sex") + 
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_jac_all <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_jac.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_jac.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state",
                  title = "Jaccard") + 
    geom_point((aes(shape = sex)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    scale_shape_manual(values=c(16,18)) +
    labs(color = "Health State", shape = "Sex") + 
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_bc_all(dss1_rare)
PCoA_hs_wu_all(dss1_rare)
PCoA_hs_uu_all(dss1_rare)
PCoA_hs_jac_all(dss1_rare)

PCoA_hs_bc_all(dss1_physeq)
PCoA_hs_wu_all(dss1_physeq)
PCoA_hs_uu_all(dss1_physeq)
PCoA_hs_jac_all(dss1_physeq)

### Per sex  ------------------------------------------------------------------------

#### Rarefaction ---------------------------------------------------------------------------

Fem_rare <- rarefy_even_depth(Females, rngseed=1, 
                              sample.size=0.95*min(sample_sums(Females)), 
                              replace=F)

Male_rare <- rarefy_even_depth(Males, rngseed=1, 
                               sample.size=0.95*min(sample_sums(Males)), 
                               replace=F)

#### Repeated Measures -----------------------------------------------------------------------

# using the same formulas as above

##### Calculations -------------------------------------

dss1_PermaMC_hs_bc_f.df <- PermaMC_hs_bc(Females)
dss1_PermaMC_hs_wu_f.df <- PermaMC_hs_wu(Females)
dss1_PermaMC_hs_uu_f.df <- PermaMC_hs_uu(Females)
dss1_PermaMC_hs_jc_f.df <- PermaMC_hs_jc(Females)

dss1_PermaMC_hs_bc_m.df <- PermaMC_hs_bc(Males)
dss1_PermaMC_hs_wu_m.df <- PermaMC_hs_wu(Males)
dss1_PermaMC_hs_uu_m.df <- PermaMC_hs_uu(Males)
dss1_PermaMC_hs_jc_m.df <- PermaMC_hs_jc(Males)

# exporting to excel

write_xlsx(dss1_PermaMC_hs_bc_f.df, "dss1_PermaMC_hs_bc_f.df.xlsx")
write_xlsx(dss1_PermaMC_hs_wu_f.df, "dss1_PermaMC_hs_wu_f.df.xlsx")
write_xlsx(dss1_PermaMC_hs_uu_f.df, "dss1_PermaMC_hs_uu_f.df.xlsx")
write_xlsx(dss1_PermaMC_hs_jc_f.df, "dss1_PermaMC_hs_jc_f.df.xlsx")

write_xlsx(dss1_PermaMC_hs_bc_m.df, "dss1_PermaMC_hs_bc_m.df.xlsx")
write_xlsx(dss1_PermaMC_hs_wu_m.df, "dss1_PermaMC_hs_wu_m.df.xlsx")
write_xlsx(dss1_PermaMC_hs_uu_m.df, "dss1_PermaMC_hs_uu_m.df.xlsx")
write_xlsx(dss1_PermaMC_hs_jc_m.df, "dss1_PermaMC_hs_jc_m.df.xlsx")

#### Ordination (PCoA) --------------------------------------------------------------------

# shape by cages

PCoA_hs_bc_fm <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Female mice - BC") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    #theme(text = element_text(size = 12)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") + 
    stat_ellipse(linewidth = 2) +
    #ggforce::geom_mark_ellipse(size = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_wu_fm <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_wu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_wu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Female mice - WU") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") + 
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_uu_fm <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_uu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_uu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Female mice - UU") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") +
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_jac_fm <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_jac.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_jac.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Female mice - JC") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") + 
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_bc_fm(Fem_rare) 
PCoA_hs_wu_fm(Fem_rare) 
PCoA_hs_uu_fm(Fem_rare) 
PCoA_hs_jac_fm(Fem_rare)

PCoA_hs_bc_fm(Females) 
PCoA_hs_wu_fm(Females) 
PCoA_hs_uu_fm(Females) 
PCoA_hs_jac_fm(Females) 

PCoA_hs_bc_m <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Male mice - BC") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") + 
    guides(color = guide_legend(order = 2), 
           shape = guide_legend(order = 1)) +
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_wu_m <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_wu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_wu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Male mice - WU") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") + 
    guides(color = guide_legend(order = 2), 
           shape = guide_legend(order = 1)) +
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_uu_m <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_uu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_uu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Male mice - UU") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") +
    guides(color = guide_legend(order = 2), 
           shape = guide_legend(order = 1)) +
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_jac_m <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_jac.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_jac.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Male mice - JC") + 
    geom_point((aes(shape = cage)), size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    #geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State", shape = "Cage") + 
    guides(color = guide_legend(order = 2), 
           shape = guide_legend(order = 1)) +
    stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}


PCoA_hs_bc_m(Male_rare) # amazing sep btwn all vs Rec
PCoA_hs_wu_m(Male_rare) # clear sep btwn Baseline & Rec
PCoA_hs_uu_m(Male_rare) # clear sep btwn all but pre-symp
PCoA_hs_jac_m(Male_rare) # clear sep btwn all but pre-symp

PCoA_hs_bc_m(Males) # separation between Baseline/Presymp/Symp vs Recovery
PCoA_hs_wu_m(Males) # separation between Baseline & Recovery (days 1-13/14 on left; rest on right)
PCoA_hs_uu_m(Males) # incredible separation between all but pre-symp
PCoA_hs_jac_m(Males) # incredible separation between all but pre-symp

### Per cage ------------------------------------------------------------------------------

#### Rarefaction ---------------------------------------------------------------------------
F1_rare <- rarefy_even_depth(F1, rngseed=1, sample.size=0.95*min(sample_sums(F1)), replace=F)
F2_rare <- rarefy_even_depth(F2, rngseed=1, sample.size=0.95*min(sample_sums(F2)), replace=F)
M1_rare <- rarefy_even_depth(M1, rngseed=1, sample.size=0.95*min(sample_sums(M1)), replace=F)
M2_rare <- rarefy_even_depth(M2, rngseed=1, sample.size=0.95*min(sample_sums(M2)), replace=F)

#### Ordination (PCoA) --------------------------------------------------------------------

PCoA_hs_bc_fm_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage F2 - Bray-Curtis") + 
    geom_point(size=5) + theme_bw() +  
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    #theme(text = element_text(size = 12)) + 
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") + 
    #stat_ellipse(linewidth = 2) +
    #ggforce::geom_mark_ellipse(size = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_wu_fm_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_wu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_wu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage F2 - WUnifrac") + 
    geom_point(size=5) + theme_bw() +  theme(text = element_text(size = 12)) +
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") + 
    #stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_uu_fm_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_uu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_uu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage F2 - Un-Unifrac") + 
    geom_point(size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") +
    #stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_jac_fm_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_jac.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_jac.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage F2 - Jaccard") + 
    geom_point(size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") + 
    #stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_bc_fm_cage(F1_rare)  
PCoA_hs_wu_fm_cage(F1_rare)  
PCoA_hs_uu_fm_cage(F1_rare)  
PCoA_hs_jac_fm_cage(F1_rare) 

PCoA_hs_bc_fm_cage(F2_rare)  
PCoA_hs_wu_fm_cage(F2_rare)  
PCoA_hs_uu_fm_cage(F2_rare)  
PCoA_hs_jac_fm_cage(F2_rare) 

PCoA_hs_bc_fm_cage(F1) 
PCoA_hs_wu_fm_cage(F1) 
PCoA_hs_uu_fm_cage(F1) 
PCoA_hs_jac_fm_cage(F1) 

PCoA_hs_bc_fm_cage(F2) 
PCoA_hs_wu_fm_cage(F2) 
PCoA_hs_uu_fm_cage(F2) 
PCoA_hs_jac_fm_cage(F2) 

PCoA_hs_bc_m_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_bc.rlog = rbiom::beta.div(counts, 'bray-curtis', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_bc.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage M2 - Bray-Curtis") + 
    geom_point(size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") + 
    #stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_wu_m_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_wu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=TRUE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_wu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage M2 - WUnifrac") + 
    geom_point(size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") + 
    #stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_uu_m_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_uu.rlog = rbiom::beta.div(counts, 'unifrac', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_uu.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage M2 - Un-Unifrac") + 
    geom_point(size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") +
    #stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_jac_m_cage <- function(phylo_exp) {
  counts = otu_table(phylo_exp)
  tree = phy_tree(phylo_exp)
  phylo_exp@sam_data[["health.state"]] <- fct_relevel(phylo_exp@sam_data[["health.state"]], 
                                                      c("Baseline", "Pre-symptomatic", 
                                                        "Symptomatic", "Recovery"))
  rbiom_jac.rlog = rbiom::beta.div(counts, 'jaccard', weighted=FALSE, tree=tree)
  pcoa.rlog = ordinate(phylo_exp, method="PCoA", distance=rbiom_jac.rlog)
  plot_ordination(phylo_exp, pcoa.rlog, color="health.state", title = "Cage M2 - Jaccard") + 
    geom_point(size=5) + theme_bw() +  theme(text = element_text(size = 12)) + 
    geom_text(aes(label = day), colour = "black", hjust=1, vjust=1) +
    scale_color_manual(values=c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
    labs(color = "Health State") + 
    #stat_ellipse(linewidth = 2) +
    theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=22), 
          axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
          legend.text=element_text(size=22), legend.title=element_text(size=22))
}

PCoA_hs_bc_m_cage(M1) 
PCoA_hs_wu_m_cage(M1) 
PCoA_hs_uu_m_cage(M1) 
PCoA_hs_jac_m_cage(M1) 

PCoA_hs_bc_m_cage(M2) 
PCoA_hs_wu_m_cage(M2) 
PCoA_hs_uu_m_cage(M2) 
PCoA_hs_jac_m_cage(M2) 
