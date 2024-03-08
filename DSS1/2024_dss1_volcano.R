# ANCOMBC2 volcano plot for DSS #1 data

# Libraries -----------------------------------------

library(tidyverse)
library(ggrepel)

# Data manipulation ---------------------------------

# data frames for ANCOMBC2 Dunnett output

dss1_phylum_out_dunn.df <- dss1_phylum_out[["res_dunn"]]
dss1_class_out_dunn.df <- dss1_class_out[["res_dunn"]]
dss1_order_out_dunn.df <- dss1_order_out[["res_dunn"]]
dss1_family_out_dunn.df <- dss1_family_out[["res_dunn"]]
dss1_genus_out_dunn.df <- dss1_genus_out[["res_dunn"]]
dss1_species_out_dunn.df <- dss1_species_out[["res_dunn"]]
dss1_asv_out_dunn.df <- dss1_asv_out[["res_dunn"]]

dss1_phylum_out_dunn.df
dss1_class_out_dunn.df
dss1_order_out_dunn.df
dss1_family_out_dunn.df
dss1_genus_out_dunn.df
dss1_species_out_dunn.df
dss1_asv_out_dunn.df

# use the below to search for specific taxa
# https://www.statology.org/r-find-value-any-column/

dss1_phy_fix_tax_table.df %>% filter_all(any_vars(. %in% c("RF39", "UCG-005")))

# create column for showing which points are 
# differentially expressed in which direction
# (up or down)

## Baseline vs Pre-symptomatic -------------------------------------

# add a column of NAs

dss1_phylum_out_dunn.df$diffexpressed_PS <- "NO"
dss1_class_out_dunn.df$diffexpressed_PS <- "NO"
dss1_order_out_dunn.df$diffexpressed_PS <- "NO"
dss1_family_out_dunn.df$diffexpressed_PS <- "NO"
dss1_genus_out_dunn.df$diffexpressed_PS <- "NO"
dss1_species_out_dunn.df$diffexpressed_PS <- "NO"
dss1_asv_out_dunn.df$diffexpressed_PS <- "NO"

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss1_phylum_out_dunn.df$diffexpressed_PS[dss1_phylum_out_dunn.df$`lfc_health.statePre-symptomatic` > 0.6 & dss1_phylum_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "UP"
dss1_class_out_dunn.df$diffexpressed_PS[dss1_class_out_dunn.df$`lfc_health.statePre-symptomatic` > 0.6 & dss1_class_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "UP"
dss1_order_out_dunn.df$diffexpressed_PS[dss1_order_out_dunn.df$`lfc_health.statePre-symptomatic` > 0.6 & dss1_order_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "UP"
dss1_family_out_dunn.df$diffexpressed_PS[dss1_family_out_dunn.df$`lfc_health.statePre-symptomatic` > 0.6 & dss1_family_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "UP"
dss1_genus_out_dunn.df$diffexpressed_PS[dss1_genus_out_dunn.df$`lfc_health.statePre-symptomatic` > 0.6 & dss1_genus_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "UP"
dss1_species_out_dunn.df$diffexpressed_PS[dss1_species_out_dunn.df$`lfc_health.statePre-symptomatic` > 0.6 & dss1_species_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "UP"
dss1_asv_out_dunn.df$diffexpressed_PS[dss1_asv_out_dunn.df$`lfc_health.statePre-symptomatic` > 0.6 & dss1_asv_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "UP"


# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss1_phylum_out_dunn.df$diffexpressed_PS[dss1_phylum_out_dunn.df$`lfc_health.statePre-symptomatic` < -0.6 & dss1_phylum_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "DOWN"
dss1_class_out_dunn.df$diffexpressed_PS[dss1_class_out_dunn.df$`lfc_health.statePre-symptomatic` < -0.6 & dss1_class_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "DOWN"
dss1_order_out_dunn.df$diffexpressed_PS[dss1_order_out_dunn.df$`lfc_health.statePre-symptomatic` < -0.6 & dss1_order_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "DOWN"
dss1_family_out_dunn.df$diffexpressed_PS[dss1_family_out_dunn.df$`lfc_health.statePre-symptomatic` < -0.6 & dss1_family_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "DOWN"
dss1_genus_out_dunn.df$diffexpressed_PS[dss1_genus_out_dunn.df$`lfc_health.statePre-symptomatic` < -0.6 & dss1_genus_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "DOWN"
dss1_species_out_dunn.df$diffexpressed_PS[dss1_species_out_dunn.df$`lfc_health.statePre-symptomatic` < -0.6 & dss1_species_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "DOWN"
dss1_asv_out_dunn.df$diffexpressed_PS[dss1_asv_out_dunn.df$`lfc_health.statePre-symptomatic` < -0.6 & dss1_asv_out_dunn.df$`q_health.statePre-symptomatic` < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss1_phylum_out_dunn.df$label_PS <- NA
dss1_phylum_out_dunn.df$label_PS[dss1_phylum_out_dunn.df$diffexpressed_PS != "NO"] <- 
  dss1_phylum_out_dunn.df$taxon[dss1_phylum_out_dunn.df$diffexpressed_PS != "NO"]

dss1_class_out_dunn.df$label_PS <- NA
dss1_class_out_dunn.df$label_PS[dss1_class_out_dunn.df$diffexpressed_PS != "NO"] <- 
  dss1_class_out_dunn.df$taxon[dss1_class_out_dunn.df$diffexpressed_PS != "NO"]

dss1_order_out_dunn.df$label_PS <- NA
dss1_order_out_dunn.df$label_PS[dss1_order_out_dunn.df$diffexpressed_PS != "NO"] <- 
  dss1_order_out_dunn.df$taxon[dss1_order_out_dunn.df$diffexpressed_PS != "NO"]

dss1_family_out_dunn.df$label_PS <- NA
dss1_family_out_dunn.df$label_PS[dss1_family_out_dunn.df$diffexpressed_PS != "NO"] <- 
  dss1_family_out_dunn.df$taxon[dss1_family_out_dunn.df$diffexpressed_PS != "NO"]

dss1_genus_out_dunn.df$label_PS <- NA
dss1_genus_out_dunn.df$label_PS[dss1_genus_out_dunn.df$diffexpressed_PS != "NO"] <- 
  dss1_genus_out_dunn.df$taxon[dss1_genus_out_dunn.df$diffexpressed_PS != "NO"]

dss1_species_out_dunn.df$label_PS <- NA
dss1_species_out_dunn.df$label_PS[dss1_species_out_dunn.df$diffexpressed_PS != "NO"] <- 
  dss1_species_out_dunn.df$taxon[dss1_species_out_dunn.df$diffexpressed_PS != "NO"]

dss1_asv_out_dunn.df$label_PS <- NA
dss1_asv_out_dunn.df$label_PS[dss1_asv_out_dunn.df$diffexpressed_PS != "NO"] <- 
  dss1_asv_out_dunn.df$taxon[dss1_asv_out_dunn.df$diffexpressed_PS != "NO"]

## Baseline vs Symptomatic -------------------------------------------

# add a column of NAs

dss1_phylum_out_dunn.df$diffexpressed_SY <- "NO"
dss1_class_out_dunn.df$diffexpressed_SY <- "NO"
dss1_order_out_dunn.df$diffexpressed_SY <- "NO"
dss1_family_out_dunn.df$diffexpressed_SY <- "NO"
dss1_genus_out_dunn.df$diffexpressed_SY <- "NO"
dss1_species_out_dunn.df$diffexpressed_SY <- "NO"
dss1_asv_out_dunn.df$diffexpressed_SY <- "NO"

dss1_phylum_out_dunn.df$lfc_health.stateSymptomatic

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss1_phylum_out_dunn.df$diffexpressed_SY[dss1_phylum_out_dunn.df$lfc_health.stateSymptomatic > 0.6 & dss1_phylum_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "UP"
dss1_class_out_dunn.df$diffexpressed_SY[dss1_class_out_dunn.df$lfc_health.stateSymptomatic > 0.6 & dss1_class_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "UP"
dss1_order_out_dunn.df$diffexpressed_SY[dss1_order_out_dunn.df$lfc_health.stateSymptomatic > 0.6 & dss1_order_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "UP"
dss1_family_out_dunn.df$diffexpressed_SY[dss1_family_out_dunn.df$lfc_health.stateSymptomatic > 0.6 & dss1_family_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "UP"
dss1_genus_out_dunn.df$diffexpressed_SY[dss1_genus_out_dunn.df$lfc_health.stateSymptomatic > 0.6 & dss1_genus_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "UP"
dss1_species_out_dunn.df$diffexpressed_SY[dss1_species_out_dunn.df$lfc_health.stateSymptomatic > 0.6 & dss1_species_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "UP"
dss1_asv_out_dunn.df$diffexpressed_SY[dss1_asv_out_dunn.df$lfc_health.stateSymptomatic > 0.6 & dss1_asv_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "UP"


# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss1_phylum_out_dunn.df$diffexpressed_SY[dss1_phylum_out_dunn.df$lfc_health.stateSymptomatic < -0.6 & dss1_phylum_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "DOWN"
dss1_class_out_dunn.df$diffexpressed_SY[dss1_class_out_dunn.df$lfc_health.stateSymptomatic < -0.6 & dss1_class_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "DOWN"
dss1_order_out_dunn.df$diffexpressed_SY[dss1_order_out_dunn.df$lfc_health.stateSymptomatic < -0.6 & dss1_order_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "DOWN"
dss1_family_out_dunn.df$diffexpressed_SY[dss1_family_out_dunn.df$lfc_health.stateSymptomatic < -0.6 & dss1_family_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "DOWN"
dss1_genus_out_dunn.df$diffexpressed_SY[dss1_genus_out_dunn.df$lfc_health.stateSymptomatic < -0.6 & dss1_genus_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "DOWN"
dss1_species_out_dunn.df$diffexpressed_SY[dss1_species_out_dunn.df$lfc_health.stateSymptomatic < -0.6 & dss1_species_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "DOWN"
dss1_asv_out_dunn.df$diffexpressed_SY[dss1_asv_out_dunn.df$lfc_health.stateSymptomatic < -0.6 & dss1_asv_out_dunn.df$q_health.stateSymptomatic < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss1_phylum_out_dunn.df$label_SY <- NA
dss1_phylum_out_dunn.df$label_SY[dss1_phylum_out_dunn.df$diffexpressed_SY != "NO"] <- 
  dss1_phylum_out_dunn.df$taxon[dss1_phylum_out_dunn.df$diffexpressed_SY != "NO"]

dss1_class_out_dunn.df$label_SY <- NA
dss1_class_out_dunn.df$label_SY[dss1_class_out_dunn.df$diffexpressed_SY != "NO"] <- 
  dss1_class_out_dunn.df$taxon[dss1_class_out_dunn.df$diffexpressed_SY != "NO"]

dss1_order_out_dunn.df$label_SY <- NA
dss1_order_out_dunn.df$label_SY[dss1_order_out_dunn.df$diffexpressed_SY != "NO"] <- 
  dss1_order_out_dunn.df$taxon[dss1_order_out_dunn.df$diffexpressed_SY != "NO"]

dss1_family_out_dunn.df$label_SY <- NA
dss1_family_out_dunn.df$label_SY[dss1_family_out_dunn.df$diffexpressed_SY != "NO"] <- 
  dss1_family_out_dunn.df$taxon[dss1_family_out_dunn.df$diffexpressed_SY != "NO"]

dss1_genus_out_dunn.df$label_SY <- NA
dss1_genus_out_dunn.df$label_SY[dss1_genus_out_dunn.df$diffexpressed_SY != "NO"] <- 
  dss1_genus_out_dunn.df$taxon[dss1_genus_out_dunn.df$diffexpressed_SY != "NO"]

dss1_species_out_dunn.df$label_SY <- NA
dss1_species_out_dunn.df$label_SY[dss1_species_out_dunn.df$diffexpressed_SY != "NO"] <- 
  dss1_species_out_dunn.df$taxon[dss1_species_out_dunn.df$diffexpressed_SY != "NO"]

dss1_asv_out_dunn.df$label_SY <- NA
dss1_asv_out_dunn.df$label_SY[dss1_asv_out_dunn.df$diffexpressed_SY != "NO"] <- 
  dss1_asv_out_dunn.df$taxon[dss1_asv_out_dunn.df$diffexpressed_SY != "NO"]

## Baseline vs Recovery ----------------------------------------------

# add a column of NAs

dss1_phylum_out_dunn.df$diffexpressed_RE <- "NO"
dss1_class_out_dunn.df$diffexpressed_RE <- "NO"
dss1_order_out_dunn.df$diffexpressed_RE <- "NO"
dss1_family_out_dunn.df$diffexpressed_RE <- "NO"
dss1_genus_out_dunn.df$diffexpressed_RE <- "NO"
dss1_species_out_dunn.df$diffexpressed_RE <- "NO"
dss1_asv_out_dunn.df$diffexpressed_RE <- "NO"

dss1_phylum_out_dunn.df$lfc_health.stateRecovery

# if lfc > 0.6 and q < 0.05, set as "UP" 

dss1_phylum_out_dunn.df$diffexpressed_RE[dss1_phylum_out_dunn.df$lfc_health.stateRecovery > 0.6 & dss1_phylum_out_dunn.df$q_health.stateRecovery < 0.05] <- "UP"
dss1_class_out_dunn.df$diffexpressed_RE[dss1_class_out_dunn.df$lfc_health.stateRecovery > 0.6 & dss1_class_out_dunn.df$q_health.stateRecovery < 0.05] <- "UP"
dss1_order_out_dunn.df$diffexpressed_RE[dss1_order_out_dunn.df$lfc_health.stateRecovery > 0.6 & dss1_order_out_dunn.df$q_health.stateRecovery < 0.05] <- "UP"
dss1_family_out_dunn.df$diffexpressed_RE[dss1_family_out_dunn.df$lfc_health.stateRecovery > 0.6 & dss1_family_out_dunn.df$q_health.stateRecovery < 0.05] <- "UP"
dss1_genus_out_dunn.df$diffexpressed_RE[dss1_genus_out_dunn.df$lfc_health.stateRecovery > 0.6 & dss1_genus_out_dunn.df$q_health.stateRecovery < 0.05] <- "UP"
dss1_species_out_dunn.df$diffexpressed_RE[dss1_species_out_dunn.df$lfc_health.stateRecovery > 0.6 & dss1_species_out_dunn.df$q_health.stateRecovery < 0.05] <- "UP"
dss1_asv_out_dunn.df$diffexpressed_RE[dss1_asv_out_dunn.df$lfc_health.stateRecovery > 0.6 & dss1_asv_out_dunn.df$q_health.stateRecovery < 0.05] <- "UP"


# if lfc < -0.6 and q < 0.05, set as "DOWN"

dss1_phylum_out_dunn.df$diffexpressed_RE[dss1_phylum_out_dunn.df$lfc_health.stateRecovery < -0.6 & dss1_phylum_out_dunn.df$q_health.stateRecovery < 0.05] <- "DOWN"
dss1_class_out_dunn.df$diffexpressed_RE[dss1_class_out_dunn.df$lfc_health.stateRecovery < -0.6 & dss1_class_out_dunn.df$q_health.stateRecovery < 0.05] <- "DOWN"
dss1_order_out_dunn.df$diffexpressed_RE[dss1_order_out_dunn.df$lfc_health.stateRecovery < -0.6 & dss1_order_out_dunn.df$q_health.stateRecovery < 0.05] <- "DOWN"
dss1_family_out_dunn.df$diffexpressed_RE[dss1_family_out_dunn.df$lfc_health.stateRecovery < -0.6 & dss1_family_out_dunn.df$q_health.stateRecovery < 0.05] <- "DOWN"
dss1_genus_out_dunn.df$diffexpressed_RE[dss1_genus_out_dunn.df$lfc_health.stateRecovery < -0.6 & dss1_genus_out_dunn.df$q_health.stateRecovery < 0.05] <- "DOWN"
dss1_species_out_dunn.df$diffexpressed_RE[dss1_species_out_dunn.df$lfc_health.stateRecovery < -0.6 & dss1_species_out_dunn.df$q_health.stateRecovery < 0.05] <- "DOWN"
dss1_asv_out_dunn.df$diffexpressed_RE[dss1_asv_out_dunn.df$lfc_health.stateRecovery < -0.6 & dss1_asv_out_dunn.df$q_health.stateRecovery < 0.05] <- "DOWN"

# create column to contain names of DA taxa

dss1_phylum_out_dunn.df$label_RE <- NA
dss1_phylum_out_dunn.df$label_RE[dss1_phylum_out_dunn.df$diffexpressed_RE != "NO"] <- 
  dss1_phylum_out_dunn.df$taxon[dss1_phylum_out_dunn.df$diffexpressed_RE != "NO"]

dss1_class_out_dunn.df$label_RE <- NA
dss1_class_out_dunn.df$label_RE[dss1_class_out_dunn.df$diffexpressed_RE != "NO"] <- 
  dss1_class_out_dunn.df$taxon[dss1_class_out_dunn.df$diffexpressed_RE != "NO"]

dss1_order_out_dunn.df$label_RE <- NA
dss1_order_out_dunn.df$label_RE[dss1_order_out_dunn.df$diffexpressed_RE != "NO"] <- 
  dss1_order_out_dunn.df$taxon[dss1_order_out_dunn.df$diffexpressed_RE != "NO"]

dss1_family_out_dunn.df$label_RE <- NA
dss1_family_out_dunn.df$label_RE[dss1_family_out_dunn.df$diffexpressed_RE != "NO"] <- 
  dss1_family_out_dunn.df$taxon[dss1_family_out_dunn.df$diffexpressed_RE != "NO"]

dss1_genus_out_dunn.df$label_RE <- NA
dss1_genus_out_dunn.df$label_RE[dss1_genus_out_dunn.df$diffexpressed_RE != "NO"] <- 
  dss1_genus_out_dunn.df$taxon[dss1_genus_out_dunn.df$diffexpressed_RE != "NO"]

dss1_species_out_dunn.df$label_RE <- NA
dss1_species_out_dunn.df$label_RE[dss1_species_out_dunn.df$diffexpressed_RE != "NO"] <- 
  dss1_species_out_dunn.df$taxon[dss1_species_out_dunn.df$diffexpressed_RE != "NO"]

dss1_asv_out_dunn.df$label_RE <- NA
dss1_asv_out_dunn.df$label_RE[dss1_asv_out_dunn.df$diffexpressed_RE != "NO"] <- 
  dss1_asv_out_dunn.df$taxon[dss1_asv_out_dunn.df$diffexpressed_RE != "NO"]

# Plotting --------------------------------------------

# as per:
# https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

# x-axis is the negative log10 of p-value (-log10)
# y-axis is the log2 of the fold change between the two conditions

# ANCOMBC2 already reports log fold change, but using natural log (ln)

## Phylum --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------

dss1_phylum_out_dunn.df %>% 
  ggplot(aes(x = `lfc_health.statePre-symptomatic`, 
             y = -log10(`q_health.statePre-symptomatic`),
                        colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
         title = "DSS #1 - Baseline vs PS - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ------------------------

dss1_phylum_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateSymptomatic, 
             y = -log10(q_health.stateSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs SY - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery ------------------------

dss1_phylum_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateRecovery, 
             y = -log10(q_health.stateRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs RE - Phylum",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Class --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------

dss1_class_out_dunn.df %>% 
  ggplot(aes(x = `lfc_health.statePre-symptomatic`, 
             y = -log10(`q_health.statePre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Pre-symptomatic - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ------------------------

dss1_class_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateSymptomatic, 
             y = -log10(q_health.stateSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Symptomatic - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery ------------------------

dss1_class_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateRecovery, 
             y = -log10(q_health.stateRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Recovery - Class",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Order --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------

dss1_order_out_dunn.df %>% 
  ggplot(aes(x = `lfc_health.statePre-symptomatic`, 
             y = -log10(`q_health.statePre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Pre-symptomatic - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ------------------------

dss1_order_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateSymptomatic, 
             y = -log10(q_health.stateSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Symptomatic - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery ------------------------

dss1_order_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateRecovery, 
             y = -log10(q_health.stateRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Recovery - Order",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Family --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------

dss1_family_out_dunn.df %>% 
  ggplot(aes(x = `lfc_health.statePre-symptomatic`, 
             y = -log10(`q_health.statePre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Pre-symptomatic - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ------------------------

dss1_family_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateSymptomatic, 
             y = -log10(q_health.stateSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Symptomatic - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery ------------------------

dss1_family_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateRecovery, 
             y = -log10(q_health.stateRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Recovery - Family",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## Genus --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------

dss1_genus_out_dunn.df %>% 
  ggplot(aes(x = `lfc_health.statePre-symptomatic`, 
             y = -log10(`q_health.statePre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "BA vs PS - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=26), 
        axis.text.x = element_text(size=20),axis.text.y=element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=26))

### Baseline vs Symptomatic ------------------------

dss1_genus_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateSymptomatic, 
             y = -log10(q_health.stateSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "BA vs SY - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=26), 
        axis.text.x = element_text(size=20),axis.text.y=element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=26))

### Baseline vs Recovery ------------------------

dss1_genus_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateRecovery, 
             y = -log10(q_health.stateRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "BA vs RE - Genus",
       col = "Diff. Exp.") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 30), axis.title=element_text(size=26), 
        axis.text.x = element_text(size=20),axis.text.y=element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=26))

## Species --------------------------------------------

### Baseline vs Pre-symptomatic ------------------------

dss1_species_out_dunn.df %>% 
  ggplot(aes(x = `lfc_health.statePre-symptomatic`, 
             y = -log10(`q_health.statePre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs PS - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ------------------------

dss1_species_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateSymptomatic, 
             y = -log10(q_health.stateSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs SY - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery ------------------------

dss1_species_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateRecovery, 
             y = -log10(q_health.stateRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs RE - Species",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

## ASV --------------------------------------------

# convering tax_table to data frame
# to explore what the ASVs are that are differentially abundant

dss1_phy_fix_tax_table.df <- data.frame(tax_table(dss1_phy_fix))

# turning row names to a column called ASV
# using tibble command rownames_to_column

dss1_phy_fix_tax_table.df <- rownames_to_column(dss1_phy_fix_tax_table.df, var = "ASV")

# then outputting to Excel
# so i can easily search it

write_xlsx(dss1_phy_fix_tax_table.df, "dss1_phy_fix_tax_table.df.xlsx")

### Baseline vs Pre-symptomatic ------------------------

dss1_asv_out_dunn.df %>% 
  ggplot(aes(x = `lfc_health.statePre-symptomatic`, 
             y = -log10(`q_health.statePre-symptomatic`),
             colour = diffexpressed_PS,
             label = label_PS)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Pre-symptomatic - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Symptomatic ------------------------

dss1_asv_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateSymptomatic, 
             y = -log10(q_health.stateSymptomatic),
             colour = diffexpressed_SY,
             label = label_SY)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Symptomatic - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))

### Baseline vs Recovery ------------------------

dss1_asv_out_dunn.df %>% 
  ggplot(aes(x = lfc_health.stateRecovery, 
             y = -log10(q_health.stateRecovery),
             colour = diffexpressed_RE,
             label = label_RE)) +
  geom_text_repel(size = 6) +
  geom_point(size = 2) + 
  scale_color_manual(values=c("#4477AA", "#000000", "#EE6677")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#EE6677") +
  geom_hline(yintercept=-log10(0.05), col="#EE6677") +
  labs(x = "Log fold change", y = "-log10(q-value)", 
       title = "DSS #1 - Baseline vs Recovery - ASV",
       col = "Diff Exp") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(hjust=0.5, size = 24), axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16),axis.text.y=element_text(size=16), 
        legend.text=element_text(size=16), legend.title=element_text(size=20))
