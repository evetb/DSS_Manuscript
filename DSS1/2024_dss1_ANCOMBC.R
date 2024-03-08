# Running ANCOM-BC2 with DSS #1 data

# Libraries ------------

library(qiime2R)
library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(writexl)
library(microViz) # to fix taxa names
library(lme4)

# QIIME2R ------------------

# Creating phyloseq object from QIIME2 output

dss1_phy <- qza_to_phyloseq(features = "table.qza",
                            tree = "rooted-tree.qza",
                            taxonomy = "taxonomy.qza",
                            metadata = "DSS1_metadata2.tsv")

# Data manipulation -----------------------------------

## manipulating tax_table w/in phyloseq object -----------------------

# using tax_fix() 
# see https://david-barnett.github.io/microViz/reference/tax_fix.html

dss1_phy_fix <- tax_fix(dss1_phy, min_length = 2, unknowns = NA)

# checking results - worked!
unique(dss1_phy_fix@tax_table[,"Genus"])
unique(dss1_phy_fix@tax_table[,"Species"])

# trying tax_fix() for "uncultured" (Genus), "unidentified", 
# "uncultured_bacterium", "uncultured_organism", 
# "uncultured_prokaryote", "uncultured_marine", 
# "uncultured_rumen" (Species)

dss1_phy_fix <- tax_fix(dss1_phy, 
                        min_length = 2, 
                        unknowns = 
                          c("uncultured", "unidentified", "uncultured_bacterium",
                            "uncultured_organism", "uncultured_prokaryote", "uncultured_marine",
                            "uncultured_rumen"))

# checking how many times certain names come up
table(dss1_phy_fix@tax_table[,"Species"])

# removing pre-symptomatic samples
# due to not enough samples for comparisons
# ps_filter() from microViz
# https://david-barnett.github.io/microViz/reference/ps_filter.html

dss1_phy_fix_filtered <- ps_filter(dss1_phy_fix, health.state != "Pre-symptomatic")

## making a column called sample ---------------------------------
## which will be numeric
## one value per cage

dss1_phy_fix_filtered@sam_data[["sample"]] <-  with(dss1_phy_fix_filtered@sam_data, 
                                      ifelse(cage == "F1", 1, 
                                             ifelse(cage == "F2", 2, 
                                                    ifelse(cage == "M1", 3, 4))))

# ANCOMBC2: all cages -----------------------------------------------

## Phylum --------------------

# rand_formula syntax:
# https://github.com/FrederickHuangLin/ANCOMBC/discussions/161

dss1_phylum_out <- ancombc2(dss1_phy_fix, 
                         assay_name = "counts",
                         tax_level = "Phylum",
                         fix_formula = "health.state",
                         rand_formula = "(1|cage)",
                         p_adj_method = "bonferroni",
                         pseudo = 10^-6,
                         pseudo_sens = TRUE,
                         prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                         group = "health.state",
                         struc_zero = TRUE,
                         alpha = 0.05, n_cl = 1, verbose = TRUE,
                         global = TRUE,
                         dunnet = TRUE,
                         iter_control = list(tol = 1e-2, max_iter = 20, 
                                             verbose = TRUE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         lme_control = lme4::lmerControl(),
                         mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

## Class ---------------------------

dss1_class_out <- ancombc2(dss1_phy_fix, 
                         assay_name = "counts",
                         tax_level = "Class",
                         fix_formula = "health.state",
                         rand_formula = "(1|cage)",
                         p_adj_method = "bonferroni",
                         pseudo = 0,
                         pseudo_sens = TRUE,
                         prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                         group = "health.state",
                         struc_zero = TRUE,
                         alpha = 0.05, n_cl = 1, verbose = TRUE,
                         global = TRUE,
                         dunnet = TRUE,
                         iter_control = list(tol = 1e-2, max_iter = 20, 
                                             verbose = TRUE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         lme_control = lme4::lmerControl(),
                         mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

## Order ---------------------------

dss1_order_out <- ancombc2(dss1_phy_fix, 
                           assay_name = "counts",
                           tax_level = "Order",
                           fix_formula = "health.state",
                           rand_formula = "(1|cage)",
                           p_adj_method = "bonferroni",
                           pseudo = 0,
                           pseudo_sens = TRUE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "health.state",
                           struc_zero = TRUE,
                           alpha = 0.05, n_cl = 1, verbose = TRUE,
                           global = TRUE,
                           dunnet = TRUE,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

## Family --------------------------

dss1_family_out <- ancombc2(dss1_phy_fix, 
                           assay_name = "counts",
                           tax_level = "Family",
                           fix_formula = "health.state",
                           rand_formula = "(1|cage)",
                           p_adj_method = "bonferroni",
                           pseudo = 0,
                           pseudo_sens = TRUE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "health.state",
                           struc_zero = TRUE,
                           alpha = 0.05, n_cl = 1, verbose = TRUE,
                           global = TRUE,
                           dunnet = TRUE,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

## Genus ---------------------------

dss1_genus_out <- ancombc2(dss1_phy_fix, 
                           assay_name = "counts",
                           tax_level = "Genus",
                           fix_formula = "health.state",
                           rand_formula = "(1|cage)",
                           p_adj_method = "bonferroni",
                           pseudo = 0,
                           pseudo_sens = FALSE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "health.state",
                           struc_zero = TRUE,
                           alpha = 0.05, n_cl = 1, verbose = TRUE,
                           global = TRUE,
                           dunnet = TRUE,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

## ASV -------------------------------

dss1_asv_out <- ancombc2(dss1_phy_fix, 
                                  assay_name = "counts",
                                  fix_formula = "health.state",
                                  rand_formula = "(1|cage)",
                                  p_adj_method = "bonferroni",
                                  pseudo = 0,
                                  pseudo_sens = TRUE,
                                  prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                  group = "health.state",
                                  struc_zero = TRUE,
                                  alpha = 0.05, n_cl = 1, verbose = TRUE,
                                  global = TRUE,
                                  dunnet = TRUE,
                                  iter_control = list(tol = 1e-2, max_iter = 1, 
                                                      verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 1),
                                  lme_control = lme4::lmerControl(),
                                  mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 1))

dss1_asv_out_filtered <- ancombc2(dss1_phy_fix_filtered, 
                           assay_name = "counts",
                           fix_formula = "health.state",
                           rand_formula = "(1|cage)",
                           p_adj_method = "bonferroni",
                           pseudo = 0,
                           pseudo_sens = TRUE,
                           prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                           group = "health.state",
                           struc_zero = TRUE,
                           alpha = 0.05, n_cl = 1, verbose = TRUE,
                           global = TRUE,
                           dunnet = TRUE,
                           iter_control = list(tol = 1e-2, max_iter = 1, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 1),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 1))


## Exporting ANCOMBC output ----------------------------

# From list to excel data frame-type object
# only exporting the dunnett's results
# since that's what I'm interested in

write_xlsx(dss1_phylum_out[["res_dunn"]], "dss1_phylum_out_dunn.xlsx")
write_xlsx(dss1_class_out[["res_dunn"]], "dss1_class_out_dunn.xlsx")
write_xlsx(dss1_order_out[["res_dunn"]], "dss1_order_out_dunn.xlsx")
write_xlsx(dss1_family_out[["res_dunn"]], "dss1_family_out_dunn.xlsx")
write_xlsx(dss1_genus_out[["res_dunn"]], "dss1_genus_out_dunn.xlsx")
write_xlsx(dss1_asv_out[["res_dunn"]], "dss1_asv_out_dunn.xlsx")
