# running ANCOM-BC2 w/DSS #2 data

# Libraries ------------

library(qiime2R)
library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(writexl)
library(microViz)

# QIIME2R ------------------

# Creating phyloseq object from QIIME2 output

dss2_phy <- qza_to_phyloseq(features = "nnn-filtered-table.qza",
                            tree = "rooted-tree.qza",
                            taxonomy = "taxonomy.qza",
                            metadata = "DSS2_metadata2.tsv")

# Data manipulation -----------------------------------

# checking unique Genus names in dss2_phy
# as per https://microbiome.github.io/tutorials/cleaning_taxonomy_table.html

sort(unique(tax_table(dss2_phy)[,"Genus"])) # 229 taxa

table(tax_table(dss2_phy)[,"Genus"]) # get number of appearances per name
# "uncultured" appears 598 times

## manipulating tax_table w/in phyloseq object -----------------------

# using tax_fix() 
# see https://david-barnett.github.io/microViz/reference/tax_fix.html

dss2_phy_fix <- tax_fix(dss2_phy, min_length = 2, unknowns = NA)

# checking results - worked!
unique(dss2_phy_fix@tax_table[,"Genus"]) # 251 taxa
unique(dss2_phy_fix@tax_table[,"Species"])

# checking how many times certain names come up
table(dss2_phy_fix@tax_table[,"Genus"]) # originally 598 cases of uncultured
table(dss2_phy_fix@tax_table[,"Species"])

# trying tax_fix() for "uncultured" (Genus), "unidentified", "uncultured_bacterium",
# "uncultured_organism", "uncultured_prokaryote", "uncultured_marine", "uncultured_rumen" (Species)

dss2_phy_fix <- tax_fix(dss2_phy, 
                        min_length = 2, 
                        unknowns = c("uncultured", "uncultured Genus", 
                                     "unidentified", "uncultured_bacterium",
                                     "uncultured_organism"))

# still have some "uncultured Genus", but that might be unavoidable

# ANCOMBC-2: per sorted fraction (sf) including the pre-symptomatic (PS) group --------------------------------------

# subsetting per sf
# excluding sample M2W8 due to contamination

Epos_phy_PS <- subset_samples(dss2_phy_fix, sf == 'EdU+')
Eneg_phy_PS <- subset_samples(dss2_phy_fix, sf == 'EdU-')
Whole_phy_PS <- dss2_phy_fix %>% 
  subset_samples(sf == 'Whole') %>% 
  subset_samples(samples != "M2W8")

## Calculations ---------------

### Phylum ------------------------

#### EdU+ ------------------------------

dss2_Epos_phylum_PS <- ancombc2(Epos_phy_PS, 
                                 assay_name = "counts",
                                 tax_level = "Phylum",
                                 fix_formula = "hs",
                                 rand_formula = "(1|cage)",
                                 p_adj_method = "bonferroni",
                                 pseudo = 0,
                                 pseudo_sens = TRUE,
                                 prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                 group = "hs",
                                 struc_zero = TRUE,
                                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                                 global = TRUE,
                                 dunnet = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

#### Whole --------------------------------

dss2_Whole_phylum_PS <- ancombc2(Whole_phy_PS, 
                                  assay_name = "counts",
                                  tax_level = "Phylum",
                                  fix_formula = "hs",
                                  rand_formula = "(1|cage)",
                                  p_adj_method = "bonferroni",
                                  pseudo = 0,
                                  pseudo_sens = TRUE,
                                  prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                  group = "hs",
                                  struc_zero = TRUE,
                                  alpha = 0.05, n_cl = 1, verbose = TRUE,
                                  global = TRUE,
                                  dunnet = TRUE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                                      verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  lme_control = lme4::lmerControl(),
                                  mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

### Class ---------------------------------

#### EdU+ ------------------------------

dss2_Epos_class_PS <- ancombc2(Epos_phy_PS, 
                                assay_name = "counts",
                                tax_level = "Class",
                                fix_formula = "hs",
                                rand_formula = "(1|cage)",
                                p_adj_method = "bonferroni",
                                pseudo = 0,
                                pseudo_sens = TRUE,
                                prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                group = "hs",
                                struc_zero = TRUE,
                                alpha = 0.05, n_cl = 1, verbose = TRUE,
                                global = TRUE,
                                dunnet = TRUE,
                                iter_control = list(tol = 1e-2, max_iter = 20, 
                                                    verbose = TRUE),
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(),
                                mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

#### Whole --------------------------------

dss2_Whole_class_PS <- ancombc2(Whole_phy_PS, 
                                 assay_name = "counts",
                                 tax_level = "Class",
                                 fix_formula = "hs",
                                 rand_formula = "(1|cage)",
                                 p_adj_method = "bonferroni",
                                 pseudo = 0,
                                 pseudo_sens = TRUE,
                                 prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                 group = "hs",
                                 struc_zero = TRUE,
                                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                                 global = TRUE,
                                 dunnet = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

### Order ---------------------------------

#### EdU+ ------------------------------

dss2_Epos_order_PS <- ancombc2(Epos_phy_PS, 
                                assay_name = "counts",
                                tax_level = "Order",
                                fix_formula = "hs",
                                rand_formula = "(1|cage)",
                                p_adj_method = "bonferroni",
                                pseudo = 0,
                                pseudo_sens = TRUE,
                                prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                group = "hs",
                                struc_zero = TRUE,
                                alpha = 0.05, n_cl = 1, verbose = TRUE,
                                global = TRUE,
                                dunnet = TRUE,
                                iter_control = list(tol = 1e-2, max_iter = 20, 
                                                    verbose = TRUE),
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(),
                                mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

#### Whole --------------------------------

dss2_Whole_order_PS <- ancombc2(Whole_phy_PS, 
                                 assay_name = "counts",
                                 tax_level = "Order",
                                 fix_formula = "hs",
                                 rand_formula = "(1|cage)",
                                 p_adj_method = "bonferroni",
                                 pseudo = 0,
                                 pseudo_sens = TRUE,
                                 prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                 group = "hs",
                                 struc_zero = TRUE,
                                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                                 global = TRUE,
                                 dunnet = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

### Family ---------------------------------

#### EdU+ ------------------------------

dss2_Epos_family_PS <- ancombc2(Epos_phy_PS, 
                                 assay_name = "counts",
                                 tax_level = "Family",
                                 fix_formula = "hs",
                                 rand_formula = "(1|cage)",
                                 p_adj_method = "bonferroni",
                                 pseudo = 0,
                                 pseudo_sens = TRUE,
                                 prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                 group = "hs",
                                 struc_zero = TRUE,
                                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                                 global = TRUE,
                                 dunnet = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

#### Whole --------------------------------

dss2_Whole_family_PS <- ancombc2(Whole_phy_PS, 
                                  assay_name = "counts",
                                  tax_level = "Family",
                                  fix_formula = "hs",
                                  rand_formula = "(1|cage)",
                                  p_adj_method = "bonferroni",
                                  pseudo = 0,
                                  pseudo_sens = TRUE,
                                  prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                  group = "hs",
                                  struc_zero = TRUE,
                                  alpha = 0.05, n_cl = 1, verbose = TRUE,
                                  global = TRUE,
                                  dunnet = TRUE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                                      verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  lme_control = lme4::lmerControl(),
                                  mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

### Genus ----------------------------------

#### EdU+ ------------------------------

dss2_Epos_genus_PS <- ancombc2(Epos_phy_PS, 
                                assay_name = "counts",
                                tax_level = "Genus",
                                fix_formula = "hs",
                                rand_formula = "(1|cage)",
                                p_adj_method = "bonferroni",
                                pseudo = 0,
                                pseudo_sens = TRUE,
                                prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                group = "hs",
                                struc_zero = TRUE,
                                alpha = 0.05, n_cl = 1, verbose = TRUE,
                                global = TRUE,
                                dunnet = TRUE,
                                iter_control = list(tol = 1e-2, max_iter = 20, 
                                                    verbose = TRUE),
                                em_control = list(tol = 1e-5, max_iter = 100),
                                lme_control = lme4::lmerControl(),
                                mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

#### Whole --------------------------------

dss2_Whole_genus_PS <- ancombc2(Whole_phy_PS, 
                                 assay_name = "counts",
                                 tax_level = "Genus",
                                 fix_formula = "hs",
                                 rand_formula = "(1|cage)",
                                 p_adj_method = "bonferroni",
                                 pseudo = 0,
                                 pseudo_sens = TRUE,
                                 prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                 group = "hs",
                                 struc_zero = TRUE,
                                 alpha = 0.05, n_cl = 1, verbose = TRUE,
                                 global = TRUE,
                                 dunnet = TRUE,
                                 iter_control = list(tol = 1e-2, max_iter = 20, 
                                                     verbose = TRUE),
                                 em_control = list(tol = 1e-5, max_iter = 100),
                                 lme_control = lme4::lmerControl(),
                                 mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

### Species -------------------------------------

#### EdU+ ------------------------------

dss2_Epos_species_PS <- ancombc2(Epos_phy_PS, 
                                  assay_name = "counts",
                                  tax_level = "Species",
                                  fix_formula = "hs",
                                  rand_formula = "(1|cage)",
                                  p_adj_method = "bonferroni",
                                  pseudo = 0,
                                  pseudo_sens = TRUE,
                                  prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                  group = "hs",
                                  struc_zero = TRUE,
                                  alpha = 0.05, n_cl = 1, verbose = TRUE,
                                  global = TRUE,
                                  dunnet = TRUE,
                                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                                      verbose = TRUE),
                                  em_control = list(tol = 1e-5, max_iter = 100),
                                  lme_control = lme4::lmerControl(),
                                  mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

#### Whole ------------------------------

dss2_Whole_species_PS <- ancombc2(Whole_phy_PS, 
                                   assay_name = "counts",
                                   tax_level = "Species",
                                   fix_formula = "hs",
                                   rand_formula = "(1|cage)",
                                   p_adj_method = "bonferroni",
                                   pseudo = 0,
                                   pseudo_sens = TRUE,
                                   prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                                   group = "hs",
                                   struc_zero = TRUE,
                                   alpha = 0.05, n_cl = 1, verbose = TRUE,
                                   global = TRUE,
                                   dunnet = TRUE,
                                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                                       verbose = TRUE),
                                   em_control = list(tol = 1e-5, max_iter = 100),
                                   lme_control = lme4::lmerControl(),
                                   mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

### ASV -------------------------------------

#### EdU+ ------------------------------

dss2_Epos_asv_PS <- ancombc2(Epos_phy_PS, 
                              assay_name = "counts",
                              #tax_level = "Phylum",
                              fix_formula = "hs",
                              rand_formula = "(1|cage)",
                              p_adj_method = "bonferroni",
                              pseudo = 0,
                              pseudo_sens = TRUE,
                              prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                              group = "hs",
                              struc_zero = TRUE,
                              alpha = 0.05, n_cl = 1, verbose = TRUE,
                              global = TRUE,
                              dunnet = TRUE,
                              iter_control = list(tol = 1e-2, max_iter = 20, 
                                                  verbose = TRUE),
                              em_control = list(tol = 1e-5, max_iter = 100),
                              lme_control = lme4::lmerControl(),
                              mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

#### Whole --------------------------------

dss2_Whole_asv_PS <- ancombc2(Whole_phy_PS, 
                               assay_name = "counts",
                               #tax_level = "Phylum",
                               fix_formula = "hs",
                               rand_formula = "(1|cage)",
                               p_adj_method = "bonferroni",
                               pseudo = 0,
                               pseudo_sens = TRUE,
                               prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                               group = "hs",
                               struc_zero = TRUE,
                               alpha = 0.05, n_cl = 1, verbose = TRUE,
                               global = TRUE,
                               dunnet = TRUE,
                               iter_control = list(tol = 1e-2, max_iter = 20, 
                                                   verbose = TRUE),
                               em_control = list(tol = 1e-5, max_iter = 100),
                               lme_control = lme4::lmerControl(),
                               mdfdr_control = list(fwer_ctrl_method = "bonferroni", B = 100))

## Exporting ANCOMBC2 results --------------------------------

# only exporting Dunnett results
# using writexl package

write_xlsx(dss2_Epos_phylum_PS[["res_dunn"]], "dss2_Epos_phylum_PS.xlsx")
write_xlsx(dss2_Whole_phylum_PS[["res_dunn"]], "dss2_Whole_phylum_PS.xlsx")

write_xlsx(dss2_Epos_class_PS[["res_dunn"]], "dss2_Epos_class_PS.xlsx")
write_xlsx(dss2_Whole_class_PS[["res_dunn"]], "dss2_Whole_class_PS.xlsx")

write_xlsx(dss2_Epos_order_PS[["res_dunn"]], "dss2_Epos_order_PS.xlsx")
write_xlsx(dss2_Whole_order_PS[["res_dunn"]], "dss2_Whole_order_PS.xlsx")

write_xlsx(dss2_Epos_family_PS[["res_dunn"]], "dss2_Epos_family_PS.xlsx")
write_xlsx(dss2_Whole_family_PS[["res_dunn"]], "dss2_Whole_family_PS.xlsx")

write_xlsx(dss2_Epos_genus_PS[["res_dunn"]], "dss2_Epos_genus_PS.xlsx")
write_xlsx(dss2_Whole_genus_PS[["res_dunn"]], "dss2_Whole_genus_PS.xlsx")

write_xlsx(dss2_Whole_species_PS[["res_dunn"]], "dss2_Whole_species_PS.xlsx")

write_xlsx(dss2_Epos_asv_PS[["res_dunn"]], "dss2_Epos_asv_PS.xlsx")
write_xlsx(dss2_Whole_asv_PS[["res_dunn"]], "dss2_Whole_asv_PS.xlsx")