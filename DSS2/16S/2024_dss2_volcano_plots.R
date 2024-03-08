# ANCOMBC2 volcano plots for DSS #2 experiment

# Libraries -----------------------------------------

library(tidyverse)
library(ggrepel)
library(writexl)

# Data manipulation ---------------------------------

# data frames for ANCOMBC2 Dunnett

## Epos ------------------------------------------------

dss2_Epos_phylum_PS_dunn.df <- dss2_Epos_phylum_PS[["res_dunn"]]
dss2_Epos_class_PS_dunn.df <- dss2_Epos_class_PS[["res_dunn"]]
dss2_Epos_order_PS_dunn.df <- dss2_Epos_order_PS[["res_dunn"]]
dss2_Epos_family_PS_dunn.df <- dss2_Epos_family_PS[["res_dunn"]]
dss2_Epos_genus_PS_dunn.df <- dss2_Epos_genus_PS[["res_dunn"]]
dss2_Epos_species_PS_dunn.df <- dss2_Epos_species_PS[["res_dunn"]]
dss2_Epos_asv_PS_dunn.df <- dss2_Epos_asv_PS[["res_dunn"]]

dss2_Epos_phylum_PS_dunn.df
dss2_Epos_class_PS_dunn.df
dss2_Epos_order_PS_dunn.df
dss2_Epos_family_PS_dunn.df
dss2_Epos_genus_PS_dunn.df
dss2_Epos_species_PS_dunn.df
dss2_Epos_asv_PS_dunn.df

## Whole ---------------------------------------------

dss2_Whole_phylum_PS_dunn.df <- dss2_Whole_phylum_PS[["res_dunn"]]
dss2_Whole_class_PS_dunn.df <- dss2_Whole_class_PS[["res_dunn"]]
dss2_Whole_order_PS_dunn.df <- dss2_Whole_order_PS[["res_dunn"]]
dss2_Whole_family_PS_dunn.df <- dss2_Whole_family_PS[["res_dunn"]]
dss2_Whole_genus_PS_dunn.df <- dss2_Whole_genus_PS[["res_dunn"]]
dss2_Whole_species_PS_dunn.df <- dss2_Whole_species_PS[["res_dunn"]]
dss2_Whole_asv_PS_dunn.df <- dss2_Whole_asv_PS[["res_dunn"]]

dss2_Whole_phylum_PS_dunn.df
dss2_Whole_class_PS_dunn.df
dss2_Whole_order_PS_dunn.df
dss2_Whole_family_PS_dunn.df
dss2_Whole_genus_PS_dunn.df
dss2_Whole_species_PS_dunn.df
dss2_Whole_asv_PS_dunn.df

## Baseline vs Pre-symptomatic -------------------------------------

### Epos ----------------------------------------------------------

# add a column of NAs

dss2_Epos_phylum_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Epos_class_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Epos_order_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Epos_family_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Epos_genus_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Epos_species_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Epos_asv_PS_dunn.df$diffexpressed_PS <- "NO"

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss2_Epos_phylum_PS_dunn.df$diffexpressed_PS[dss2_Epos_phylum_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Epos_phylum_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Epos_class_PS_dunn.df$diffexpressed_PS[dss2_Epos_class_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Epos_class_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Epos_order_PS_dunn.df$diffexpressed_PS[dss2_Epos_order_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Epos_order_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Epos_family_PS_dunn.df$diffexpressed_PS[dss2_Epos_family_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Epos_family_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Epos_genus_PS_dunn.df$diffexpressed_PS[dss2_Epos_genus_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Epos_genus_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Epos_species_PS_dunn.df$diffexpressed_PS[dss2_Epos_species_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Epos_species_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Epos_asv_PS_dunn.df$diffexpressed_PS[dss2_Epos_asv_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Epos_asv_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"

# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss2_Epos_phylum_PS_dunn.df$diffexpressed_PS[dss2_Epos_phylum_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Epos_phylum_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Epos_class_PS_dunn.df$diffexpressed_PS[dss2_Epos_class_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Epos_class_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Epos_order_PS_dunn.df$diffexpressed_PS[dss2_Epos_order_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Epos_order_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Epos_family_PS_dunn.df$diffexpressed_PS[dss2_Epos_family_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Epos_family_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Epos_genus_PS_dunn.df$diffexpressed_PS[dss2_Epos_genus_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Epos_genus_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Epos_species_PS_dunn.df$diffexpressed_PS[dss2_Epos_species_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Epos_species_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Epos_asv_PS_dunn.df$diffexpressed_PS[dss2_Epos_asv_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Epos_asv_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss2_Epos_phylum_PS_dunn.df$label_PS <- NA
dss2_Epos_phylum_PS_dunn.df$label_PS[dss2_Epos_phylum_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Epos_phylum_PS_dunn.df$taxon[dss2_Epos_phylum_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Epos_class_PS_dunn.df$label_PS <- NA
dss2_Epos_class_PS_dunn.df$label_PS[dss2_Epos_class_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Epos_class_PS_dunn.df$taxon[dss2_Epos_class_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Epos_order_PS_dunn.df$label_PS <- NA
dss2_Epos_order_PS_dunn.df$label_PS[dss2_Epos_order_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Epos_order_PS_dunn.df$taxon[dss2_Epos_order_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Epos_family_PS_dunn.df$label_PS <- NA
dss2_Epos_family_PS_dunn.df$label_PS[dss2_Epos_family_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Epos_family_PS_dunn.df$taxon[dss2_Epos_family_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Epos_genus_PS_dunn.df$label_PS <- NA
dss2_Epos_genus_PS_dunn.df$label_PS[dss2_Epos_genus_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Epos_genus_PS_dunn.df$taxon[dss2_Epos_genus_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Epos_species_PS_dunn.df$label_PS <- NA
dss2_Epos_species_PS_dunn.df$label_PS[dss2_Epos_species_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Epos_species_PS_dunn.df$taxon[dss2_Epos_species_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Epos_asv_PS_dunn.df$label_PS <- NA
dss2_Epos_asv_PS_dunn.df$label_PS[dss2_Epos_asv_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Epos_asv_PS_dunn.df$taxon[dss2_Epos_asv_PS_dunn.df$diffexpressed_PS != "NO"]

### Whole ----------------------------------------------------------

# add a column of NAs

dss2_Whole_phylum_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Whole_class_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Whole_order_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Whole_family_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Whole_genus_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Whole_species_PS_dunn.df$diffexpressed_PS <- "NO"
dss2_Whole_asv_PS_dunn.df$diffexpressed_PS <- "NO"

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss2_Whole_phylum_PS_dunn.df$diffexpressed_PS[dss2_Whole_phylum_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Whole_phylum_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Whole_class_PS_dunn.df$diffexpressed_PS[dss2_Whole_class_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Whole_class_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Whole_order_PS_dunn.df$diffexpressed_PS[dss2_Whole_order_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Whole_order_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Whole_family_PS_dunn.df$diffexpressed_PS[dss2_Whole_family_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Whole_family_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Whole_genus_PS_dunn.df$diffexpressed_PS[dss2_Whole_genus_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Whole_genus_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Whole_species_PS_dunn.df$diffexpressed_PS[dss2_Whole_species_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Whole_species_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"
dss2_Whole_asv_PS_dunn.df$diffexpressed_PS[dss2_Whole_asv_PS_dunn.df$`lfc_hsPre-symptomatic` > 0.6 & dss2_Whole_asv_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "UP"

# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss2_Whole_phylum_PS_dunn.df$diffexpressed_PS[dss2_Whole_phylum_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Whole_phylum_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Whole_class_PS_dunn.df$diffexpressed_PS[dss2_Whole_class_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Whole_class_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Whole_order_PS_dunn.df$diffexpressed_PS[dss2_Whole_order_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Whole_order_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Whole_family_PS_dunn.df$diffexpressed_PS[dss2_Whole_family_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Whole_family_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Whole_genus_PS_dunn.df$diffexpressed_PS[dss2_Whole_genus_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Whole_genus_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Whole_species_PS_dunn.df$diffexpressed_PS[dss2_Whole_species_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Whole_species_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"
dss2_Whole_asv_PS_dunn.df$diffexpressed_PS[dss2_Whole_asv_PS_dunn.df$`lfc_hsPre-symptomatic` < -0.6 & dss2_Whole_asv_PS_dunn.df$`q_hsPre-symptomatic` < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss2_Whole_phylum_PS_dunn.df$label_PS <- NA
dss2_Whole_phylum_PS_dunn.df$label_PS[dss2_Whole_phylum_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Whole_phylum_PS_dunn.df$taxon[dss2_Whole_phylum_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Whole_class_PS_dunn.df$label_PS <- NA
dss2_Whole_class_PS_dunn.df$label_PS[dss2_Whole_class_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Whole_class_PS_dunn.df$taxon[dss2_Whole_class_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Whole_order_PS_dunn.df$label_PS <- NA
dss2_Whole_order_PS_dunn.df$label_PS[dss2_Whole_order_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Whole_order_PS_dunn.df$taxon[dss2_Whole_order_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Whole_family_PS_dunn.df$label_PS <- NA
dss2_Whole_family_PS_dunn.df$label_PS[dss2_Whole_family_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Whole_family_PS_dunn.df$taxon[dss2_Whole_family_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Whole_genus_PS_dunn.df$label_PS <- NA
dss2_Whole_genus_PS_dunn.df$label_PS[dss2_Whole_genus_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Whole_genus_PS_dunn.df$taxon[dss2_Whole_genus_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Whole_species_PS_dunn.df$label_PS <- NA
dss2_Whole_species_PS_dunn.df$label_PS[dss2_Whole_species_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Whole_species_PS_dunn.df$taxon[dss2_Whole_species_PS_dunn.df$diffexpressed_PS != "NO"]

dss2_Whole_asv_PS_dunn.df$label_PS <- NA
dss2_Whole_asv_PS_dunn.df$label_PS[dss2_Whole_asv_PS_dunn.df$diffexpressed_PS != "NO"] <- 
  dss2_Whole_asv_PS_dunn.df$taxon[dss2_Whole_asv_PS_dunn.df$diffexpressed_PS != "NO"]

## Baseline vs Symptomatic -------------------------------------------

### Epos ----------------------------------------------------------

# add a column of NAs

dss2_Epos_phylum_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Epos_class_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Epos_order_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Epos_family_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Epos_genus_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Epos_species_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Epos_asv_PS_dunn.df$diffexpressed_SY <- "NO"

dss2_Epos_phylum_PS_dunn.df$lfc_hsSymptomatic

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss2_Epos_phylum_PS_dunn.df$diffexpressed_SY[dss2_Epos_phylum_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Epos_phylum_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Epos_class_PS_dunn.df$diffexpressed_SY[dss2_Epos_class_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Epos_class_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Epos_order_PS_dunn.df$diffexpressed_SY[dss2_Epos_order_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Epos_order_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Epos_family_PS_dunn.df$diffexpressed_SY[dss2_Epos_family_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Epos_family_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Epos_genus_PS_dunn.df$diffexpressed_SY[dss2_Epos_genus_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Epos_genus_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Epos_species_PS_dunn.df$diffexpressed_SY[dss2_Epos_species_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Epos_species_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Epos_asv_PS_dunn.df$diffexpressed_SY[dss2_Epos_asv_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Epos_asv_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"

# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss2_Epos_phylum_PS_dunn.df$diffexpressed_SY[dss2_Epos_phylum_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Epos_phylum_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Epos_class_PS_dunn.df$diffexpressed_SY[dss2_Epos_class_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Epos_class_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Epos_order_PS_dunn.df$diffexpressed_SY[dss2_Epos_order_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Epos_order_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Epos_family_PS_dunn.df$diffexpressed_SY[dss2_Epos_family_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Epos_family_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Epos_genus_PS_dunn.df$diffexpressed_SY[dss2_Epos_genus_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Epos_genus_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Epos_species_PS_dunn.df$diffexpressed_SY[dss2_Epos_species_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Epos_species_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Epos_asv_PS_dunn.df$diffexpressed_SY[dss2_Epos_asv_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Epos_asv_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss2_Epos_phylum_PS_dunn.df$label_SY <- NA
dss2_Epos_phylum_PS_dunn.df$label_SY[dss2_Epos_phylum_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Epos_phylum_PS_dunn.df$taxon[dss2_Epos_phylum_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Epos_class_PS_dunn.df$label_SY <- NA
dss2_Epos_class_PS_dunn.df$label_SY[dss2_Epos_class_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Epos_class_PS_dunn.df$taxon[dss2_Epos_class_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Epos_order_PS_dunn.df$label_SY <- NA
dss2_Epos_order_PS_dunn.df$label_SY[dss2_Epos_order_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Epos_order_PS_dunn.df$taxon[dss2_Epos_order_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Epos_family_PS_dunn.df$label_SY <- NA
dss2_Epos_family_PS_dunn.df$label_SY[dss2_Epos_family_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Epos_family_PS_dunn.df$taxon[dss2_Epos_family_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Epos_genus_PS_dunn.df$label_SY <- NA
dss2_Epos_genus_PS_dunn.df$label_SY[dss2_Epos_genus_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Epos_genus_PS_dunn.df$taxon[dss2_Epos_genus_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Epos_species_PS_dunn.df$label_SY <- NA
dss2_Epos_species_PS_dunn.df$label_SY[dss2_Epos_species_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Epos_species_PS_dunn.df$taxon[dss2_Epos_species_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Epos_asv_PS_dunn.df$label_SY <- NA
dss2_Epos_asv_PS_dunn.df$label_SY[dss2_Epos_asv_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Epos_asv_PS_dunn.df$taxon[dss2_Epos_asv_PS_dunn.df$diffexpressed_SY != "NO"]

### Whole ----------------------------------------------------------

# add a column of NAs

dss2_Whole_phylum_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Whole_class_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Whole_order_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Whole_family_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Whole_genus_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Whole_species_PS_dunn.df$diffexpressed_SY <- "NO"
dss2_Whole_asv_PS_dunn.df$diffexpressed_SY <- "NO"

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss2_Whole_phylum_PS_dunn.df$diffexpressed_SY[dss2_Whole_phylum_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Whole_phylum_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Whole_class_PS_dunn.df$diffexpressed_SY[dss2_Whole_class_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Whole_class_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Whole_order_PS_dunn.df$diffexpressed_SY[dss2_Whole_order_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Whole_order_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Whole_family_PS_dunn.df$diffexpressed_SY[dss2_Whole_family_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Whole_family_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Whole_genus_PS_dunn.df$diffexpressed_SY[dss2_Whole_genus_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Whole_genus_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Whole_species_PS_dunn.df$diffexpressed_SY[dss2_Whole_species_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Whole_species_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"
dss2_Whole_asv_PS_dunn.df$diffexpressed_SY[dss2_Whole_asv_PS_dunn.df$lfc_hsSymptomatic > 0.6 & dss2_Whole_asv_PS_dunn.df$q_hsSymptomatic < 0.05] <- "UP"

# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss2_Whole_phylum_PS_dunn.df$diffexpressed_SY[dss2_Whole_phylum_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Whole_phylum_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Whole_class_PS_dunn.df$diffexpressed_SY[dss2_Whole_class_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Whole_class_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Whole_order_PS_dunn.df$diffexpressed_SY[dss2_Whole_order_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Whole_order_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Whole_family_PS_dunn.df$diffexpressed_SY[dss2_Whole_family_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Whole_family_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Whole_genus_PS_dunn.df$diffexpressed_SY[dss2_Whole_genus_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Whole_genus_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Whole_species_PS_dunn.df$diffexpressed_SY[dss2_Whole_species_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Whole_species_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"
dss2_Whole_asv_PS_dunn.df$diffexpressed_SY[dss2_Whole_asv_PS_dunn.df$lfc_hsSymptomatic < -0.6 & dss2_Whole_asv_PS_dunn.df$q_hsSymptomatic < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss2_Whole_phylum_PS_dunn.df$label_SY <- NA
dss2_Whole_phylum_PS_dunn.df$label_SY[dss2_Whole_phylum_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Whole_phylum_PS_dunn.df$taxon[dss2_Whole_phylum_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Whole_class_PS_dunn.df$label_SY <- NA
dss2_Whole_class_PS_dunn.df$label_SY[dss2_Whole_class_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Whole_class_PS_dunn.df$taxon[dss2_Whole_class_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Whole_order_PS_dunn.df$label_SY <- NA
dss2_Whole_order_PS_dunn.df$label_SY[dss2_Whole_order_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Whole_order_PS_dunn.df$taxon[dss2_Whole_order_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Whole_family_PS_dunn.df$label_SY <- NA
dss2_Whole_family_PS_dunn.df$label_SY[dss2_Whole_family_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Whole_family_PS_dunn.df$taxon[dss2_Whole_family_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Whole_genus_PS_dunn.df$label_SY <- NA
dss2_Whole_genus_PS_dunn.df$label_SY[dss2_Whole_genus_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Whole_genus_PS_dunn.df$taxon[dss2_Whole_genus_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Whole_species_PS_dunn.df$label_SY <- NA
dss2_Whole_species_PS_dunn.df$label_SY[dss2_Whole_species_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Whole_species_PS_dunn.df$taxon[dss2_Whole_species_PS_dunn.df$diffexpressed_SY != "NO"]

dss2_Whole_asv_PS_dunn.df$label_SY <- NA
dss2_Whole_asv_PS_dunn.df$label_SY[dss2_Whole_asv_PS_dunn.df$diffexpressed_SY != "NO"] <- 
  dss2_Whole_asv_PS_dunn.df$taxon[dss2_Whole_asv_PS_dunn.df$diffexpressed_SY != "NO"]

## Baseline vs Recovery ----------------------------------------------

### Epos ----------------------------------------------------------

# add a column of NAs

dss2_Epos_phylum_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Epos_class_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Epos_order_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Epos_family_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Epos_genus_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Epos_species_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Epos_asv_PS_dunn.df$diffexpressed_RE <- "NO"

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss2_Epos_phylum_PS_dunn.df$diffexpressed_RE[dss2_Epos_phylum_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Epos_phylum_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Epos_class_PS_dunn.df$diffexpressed_RE[dss2_Epos_class_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Epos_class_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Epos_order_PS_dunn.df$diffexpressed_RE[dss2_Epos_order_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Epos_order_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Epos_family_PS_dunn.df$diffexpressed_RE[dss2_Epos_family_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Epos_family_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Epos_genus_PS_dunn.df$diffexpressed_RE[dss2_Epos_genus_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Epos_genus_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Epos_species_PS_dunn.df$diffexpressed_RE[dss2_Epos_species_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Epos_species_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Epos_asv_PS_dunn.df$diffexpressed_RE[dss2_Epos_asv_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Epos_asv_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"

# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss2_Epos_phylum_PS_dunn.df$diffexpressed_RE[dss2_Epos_phylum_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Epos_phylum_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Epos_class_PS_dunn.df$diffexpressed_RE[dss2_Epos_class_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Epos_class_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Epos_order_PS_dunn.df$diffexpressed_RE[dss2_Epos_order_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Epos_order_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Epos_family_PS_dunn.df$diffexpressed_RE[dss2_Epos_family_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Epos_family_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Epos_genus_PS_dunn.df$diffexpressed_RE[dss2_Epos_genus_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Epos_genus_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Epos_species_PS_dunn.df$diffexpressed_RE[dss2_Epos_species_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Epos_species_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Epos_asv_PS_dunn.df$diffexpressed_RE[dss2_Epos_asv_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Epos_asv_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss2_Epos_phylum_PS_dunn.df$label_RE <- NA
dss2_Epos_phylum_PS_dunn.df$label_RE[dss2_Epos_phylum_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Epos_phylum_PS_dunn.df$taxon[dss2_Epos_phylum_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Epos_class_PS_dunn.df$label_RE <- NA
dss2_Epos_class_PS_dunn.df$label_RE[dss2_Epos_class_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Epos_class_PS_dunn.df$taxon[dss2_Epos_class_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Epos_order_PS_dunn.df$label_RE <- NA
dss2_Epos_order_PS_dunn.df$label_RE[dss2_Epos_order_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Epos_order_PS_dunn.df$taxon[dss2_Epos_order_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Epos_family_PS_dunn.df$label_RE <- NA
dss2_Epos_family_PS_dunn.df$label_RE[dss2_Epos_family_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Epos_family_PS_dunn.df$taxon[dss2_Epos_family_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Epos_genus_PS_dunn.df$label_RE <- NA
dss2_Epos_genus_PS_dunn.df$label_RE[dss2_Epos_genus_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Epos_genus_PS_dunn.df$taxon[dss2_Epos_genus_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Epos_species_PS_dunn.df$label_RE <- NA
dss2_Epos_species_PS_dunn.df$label_RE[dss2_Epos_species_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Epos_species_PS_dunn.df$taxon[dss2_Epos_species_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Epos_asv_PS_dunn.df$label_RE <- NA
dss2_Epos_asv_PS_dunn.df$label_RE[dss2_Epos_asv_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Epos_asv_PS_dunn.df$taxon[dss2_Epos_asv_PS_dunn.df$diffexpressed_RE != "NO"]

### Whole ----------------------------------------------------------

# add a column of NAs

dss2_Whole_phylum_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Whole_class_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Whole_order_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Whole_family_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Whole_genus_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Whole_species_PS_dunn.df$diffexpressed_RE <- "NO"
dss2_Whole_asv_PS_dunn.df$diffexpressed_RE <- "NO"

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss2_Whole_phylum_PS_dunn.df$diffexpressed_RE[dss2_Whole_phylum_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Whole_phylum_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Whole_class_PS_dunn.df$diffexpressed_RE[dss2_Whole_class_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Whole_class_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Whole_order_PS_dunn.df$diffexpressed_RE[dss2_Whole_order_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Whole_order_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Whole_family_PS_dunn.df$diffexpressed_RE[dss2_Whole_family_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Whole_family_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Whole_genus_PS_dunn.df$diffexpressed_RE[dss2_Whole_genus_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Whole_genus_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Whole_species_PS_dunn.df$diffexpressed_RE[dss2_Whole_species_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Whole_species_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"
dss2_Whole_asv_PS_dunn.df$diffexpressed_RE[dss2_Whole_asv_PS_dunn.df$lfc_hsRecovery > 0.6 & dss2_Whole_asv_PS_dunn.df$q_hsRecovery < 0.05] <- "UP"

# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss2_Whole_phylum_PS_dunn.df$diffexpressed_RE[dss2_Whole_phylum_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Whole_phylum_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Whole_class_PS_dunn.df$diffexpressed_RE[dss2_Whole_class_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Whole_class_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Whole_order_PS_dunn.df$diffexpressed_RE[dss2_Whole_order_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Whole_order_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Whole_family_PS_dunn.df$diffexpressed_RE[dss2_Whole_family_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Whole_family_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Whole_genus_PS_dunn.df$diffexpressed_RE[dss2_Whole_genus_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Whole_genus_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Whole_species_PS_dunn.df$diffexpressed_RE[dss2_Whole_species_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Whole_species_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"
dss2_Whole_asv_PS_dunn.df$diffexpressed_RE[dss2_Whole_asv_PS_dunn.df$lfc_hsRecovery < -0.6 & dss2_Whole_asv_PS_dunn.df$q_hsRecovery < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss2_Whole_phylum_PS_dunn.df$label_RE <- NA
dss2_Whole_phylum_PS_dunn.df$label_RE[dss2_Whole_phylum_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Whole_phylum_PS_dunn.df$taxon[dss2_Whole_phylum_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Whole_class_PS_dunn.df$label_RE <- NA
dss2_Whole_class_PS_dunn.df$label_RE[dss2_Whole_class_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Whole_class_PS_dunn.df$taxon[dss2_Whole_class_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Whole_order_PS_dunn.df$label_RE <- NA
dss2_Whole_order_PS_dunn.df$label_RE[dss2_Whole_order_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Whole_order_PS_dunn.df$taxon[dss2_Whole_order_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Whole_family_PS_dunn.df$label_RE <- NA
dss2_Whole_family_PS_dunn.df$label_RE[dss2_Whole_family_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Whole_family_PS_dunn.df$taxon[dss2_Whole_family_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Whole_genus_PS_dunn.df$label_RE <- NA
dss2_Whole_genus_PS_dunn.df$label_RE[dss2_Whole_genus_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Whole_genus_PS_dunn.df$taxon[dss2_Whole_genus_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Whole_species_PS_dunn.df$label_RE <- NA
dss2_Whole_species_PS_dunn.df$label_RE[dss2_Whole_species_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Whole_species_PS_dunn.df$taxon[dss2_Whole_species_PS_dunn.df$diffexpressed_RE != "NO"]

dss2_Whole_asv_PS_dunn.df$label_RE <- NA
dss2_Whole_asv_PS_dunn.df$label_RE[dss2_Whole_asv_PS_dunn.df$diffexpressed_RE != "NO"] <- 
  dss2_Whole_asv_PS_dunn.df$taxon[dss2_Whole_asv_PS_dunn.df$diffexpressed_RE != "NO"]

# Plotting --------------------------------------------

# as per:
# https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

# x-axis is the negative log10 of p-value (-log10)
# y-axis is the log2 of the fold change between the two conditions

# ANCOMBC2 already reports log fold change, but using natural log (ln)

## Phylum --------------------------------------------

### Baseline vs Pre-symptomatic ----------------------------

#### Epos ----------------------------------------------------

dss2_Epos_phylum_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs PS - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_phylum_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs PS - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_phylum_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs SY - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_phylum_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000",  "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs SY - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery --------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_phylum_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs RE - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_phylum_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs RE - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Class --------------------------------------------

### Baseline vs Pre-symptomatic ---------------------

#### Epos ----------------------------------------------------

dss2_Epos_class_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs PS - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_class_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs PS - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic -----------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_class_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs SY - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_class_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000",  "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs SY - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery --------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_class_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs RE - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_class_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs RE - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Order --------------------------------------------

### Baseline vs Pre-symptomatic --------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_order_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs PS - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_order_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs PS - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ---------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_order_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs SY - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_order_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000",  "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs SY - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery --------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_order_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs RE - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_order_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs RE - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Family -------------------------------------------------

### Baseline vs Pre-symptomatic -------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_family_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs PS - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_family_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000",  "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs PS - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ----------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_family_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs SY - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_family_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs SY - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery ------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_family_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs RE - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_family_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs RE - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Genus --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_genus_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "EdU+ BA vs PS - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 26), axis.title=element_text(size=20), 
        axis.text = element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=24))

#### Whole ---------------------------------------------------

dss2_Whole_genus_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "Whole BA vs PS - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 26), axis.title=element_text(size=20), 
        axis.text = element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=24))

### Baseline vs Symptomatic ----------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_genus_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "EdU+ BA vs SY - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 26), axis.title=element_text(size=20), 
        axis.text = element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=24))

#### Whole ---------------------------------------------------

dss2_Whole_genus_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "Whole BA vs SY - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 26), axis.title=element_text(size=20), 
        axis.text = element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=24))

### Baseline vs Recovery --------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_genus_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "EdU+ BA vs RE - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 26), axis.title=element_text(size=20), 
        axis.text = element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=24))

#### Whole ---------------------------------------------------

dss2_Whole_genus_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "Whole BA vs RE - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 26), axis.title=element_text(size=20), 
        axis.text = element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=24))

## Species --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_species_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs PS - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_species_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs PS - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ----------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_species_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs SY - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_species_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs SY - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery --------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_species_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs RE - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_species_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs RE - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## ASV -------------------------------------------------------

# convering tax_table to data frame for Epos
# to explore what the ASVs are that are differentially abundant

Epos_phy_PS_tax_table.df <- data.frame(tax_table(Epos_phy_PS))

# turning row names to a column called ASV
# using tibble command rownames_to_column

Epos_phy_PS_tax_table.df <- rownames_to_column(Epos_phy_PS_tax_table.df, var = "ASV")

# then outputting to Excel
# so i can easily search it

write_xlsx(Epos_phy_PS_tax_table.df, "Epos_phy_PS_tax_table.df.xlsx")

# e82fddd512bc2aba38382d775ad36a64 = Akkermansia_muciniphila
# 8a564e450fc3fc676e7fa23c4703ccc3 = Erysipelatoclostridium Genus
# 18a1fde78282c2abde37741593fe175a = Muribaculaceae Genus
# 98766337b62e23988fae0a9113c2de0c = Muribaculaceae Genus
# 627400a3f8e8cb62bdd9246f06a746cd = Muribaculaceae Genus

### Baseline vs Pre-symptomatic ------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_asv_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs PS - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_asv_PS_dunn.df %>% 
  ggplot(aes(x = `lfc_hsPre-symptomatic`, 
             y = -log10(`q_hsPre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs PS - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic -----------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_asv_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs SY - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_asv_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsSymptomatic, 
             y = -log10(q_hsSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000",  "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs SY - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery -------------------------------------

#### Epos ----------------------------------------------------

dss2_Epos_asv_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Epos - Baseline vs RE - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

#### Whole ---------------------------------------------------

dss2_Whole_asv_PS_dunn.df %>% 
  ggplot(aes(x = lfc_hsRecovery, 
             y = -log10(q_hsRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #2 - Whole - Baseline vs RE - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))
