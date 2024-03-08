# relative change in abundances for differentially abundant species
# from ANCOMBC2 results on DSS #2 data

# Libraries ------------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(RColorBrewer)

# Data Manipulation ----------------------------------------------------------------

## Phyloseq object to genus level ---------------------------------------

dss2_phy_fix

dss2_phy_genus <- tax_glom(dss2_phy_fix, taxrank = "Genus")

dss2_phy_genus.df <- psmelt(dss2_phy_genus)

## Merge taxonomy names ---------------------------------------------------------
# into single Classification column
# https://www.marsja.se/how-to-concatenate-two-columns-or-more-in-r-stringr-tidyr/

dss2_phy_genus.df_merged <- dss2_phy_genus.df %>% 
  unite("Classification", Phylum:Genus, remove = FALSE)

## Plot-friendly dataframe (df) ----------------------------------------------------------------

# make new df with only Sample, 
# Classification, and Abundance
# and the relevant metadata

# keeping Genus for easier plotting

dss2_phy_genus.df_classification <- dss2_phy_genus.df_merged %>% 
  select(Sample, Abundance, Classification, Genus, day, hs, cage, sex, sf) %>% 
  filter(sf != "EdU-")

# create new df per sorted fraction

dss2_phy_genus.df_whole <- dss2_phy_genus.df_classification %>% 
  filter(sf == "Whole")

dss2_phy_genus.df_epos <- dss2_phy_genus.df_classification %>% 
  filter(sf == "EdU+")

# Plotting ------------------------------------------------------------------

# Here is the list of the 
# differentially abundant taxa 
# I need to plot:

# Akkermansia
# Erysipelatoclostridium
# Eubacterium_fissicatena_group
# Turicibacter
# Lachnoclostridium
# Roseburia
# Blautia

# A good way to get Genus names easily from the df:
# (in alphabetical order)
sort(unique(dss2_phy_genus.df_classification$Genus))

## Line plots ----------------------------------------------------------------

### Template -------------------------------------------------------------

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "Genus"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(facets = vars(sf)) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Erysipelatoclostridium"), " species"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Plots ---------------------------------------------------------------------

# a good way to see variability
# between the cages

# making day into numeric for proper plotting
dss2_phy_genus.df_classification$day <- as.numeric(dss2_phy_genus.df_classification$day)

# specifying factor levels for health states (hs)

dss2_phy_genus.df_classification$hs <- factor(dss2_phy_genus.df_classification$hs, 
                                              levels = c("Baseline", "Pre-symptomatic",
                                                         "Symptomatic", "Recovery"))

#### Both EdU+ and Whole together ----------------------------------------------------

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "Erysipelatoclostridium"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(facets = vars(sf)) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Erysipelatoclostridium"), " species"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "[Eubacterium]_fissicatena_group"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(facets = vars(sf)) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Eubacterium fissicatena"), " group"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "Akkermansia"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  facet_wrap(facets = vars(sf)) +
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Akkermansia muciniphila")))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "Turicibacter"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  facet_wrap(facets = vars(sf)) +
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Turicibacter"), " species"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "Lachnoclostridium"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  facet_wrap(facets = vars(sf)) +
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Lachnoclostridium"), " species"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "Roseburia"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  facet_wrap(facets = vars(sf)) +
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Roseburia"), " species"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_classification, Genus %in% "Blautia"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  facet_wrap(facets = vars(sf)) +
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Blautia"), " species"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

#### EdU+ alone -----------------------------------------------------------------

# Akkermansia
# Erysipelatoclostridium

dss2_phy_genus.df_epos

dss2_phy_genus.df_epos$hs <- factor(dss2_phy_genus.df_epos$hs,
                                    levels = c("Baseline", "Pre-symptomatic",
                                               "Symptomatic", "Recovery"))

ggplot(subset(dss2_phy_genus.df_epos, Genus %in% "Akkermansia"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Akkermansia"), " species - EdU+"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_epos, Genus %in% "Erysipelatoclostridium"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Erysipelatoclostridium"), " species - EdU+"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

#### Whole alone ----------------------------------------------------------------

dss2_phy_genus.df_whole 

dss2_phy_genus.df_whole$hs <- factor(dss2_phy_genus.df_whole$hs,
                                    levels = c("Baseline", "Pre-symptomatic",
                                               "Symptomatic", "Recovery"))

# Akkermansia
# Erysipelatoclostridium
# Eubacterium_fissicatena_group
# Turicibacter
# Lachnoclostridium
# Roseburia
# Blautia

ggplot(subset(dss2_phy_genus.df_whole, Genus %in% "Akkermansia"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Akkermansia"), " species - Whole"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_whole, Genus %in% "Erysipelatoclostridium"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Erysipelatoclostridium"), " species - Whole"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_whole, Genus %in% "[Eubacterium]_fissicatena_group"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Eubacterium fissicatena"), " group - Whole"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_whole, Genus %in% "Turicibacter"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Turicibacter"), " species - Whole"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_whole, Genus %in% "Lachnoclostridium"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Lachnoclostridium"), " species - Whole"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_whole, Genus %in% "Roseburia"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Roseburia"), " species - Whole"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(dss2_phy_genus.df_whole, Genus %in% "Blautia"), 
       aes(x=day, y=Abundance, group = cage, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 3) + 
  labs(x = "Day", y = "Relative Abundance",
       title =  (expression(paste(italic("Blautia"), " species - Whole"))),
       color = "Health State", shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Heat maps -----------------------------------------------------------------

### Data preparation --------------------------------------------------------------

# I already have a dataframe per sorted fraction -
# dss2_phy_genus.df_epos and
# dss2_phy_genus.df_whole

# getting the median of each bacterium per day, grouping days by using list()
# https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame

# For dss2_phy_genus.df_epos
# grouping by Genus and by day
# https://r-coder.com/aggregate-r/
median_epos <- aggregate(dss2_phy_genus.df_epos$Abundance,
                         list(dss2_phy_genus.df_epos$Genus, 
                              dss2_phy_genus.df_epos$day),median)

# For dss2_phy_genus.df_whole
median_whole <- aggregate(dss2_phy_genus.df_whole$Abundance,
                          list(dss2_phy_genus.df_whole$Genus,
                               dss2_phy_genus.df_whole$day),median) 

# changing column name from Group.1 to Genus
# and from Group.2 to day
# and from x to Abundance

median_epos <- median_epos %>% 
  rename("Genus" = "Group.1") %>%
  rename("day" = "Group.2") %>% 
  rename("Abundance" = "x")

median_whole <- median_whole %>% 
  rename("Genus" = "Group.1") %>%
  rename("day" = "Group.2") %>% 
  rename("Abundance" = "x")

### Prep - data transformation -----------------------------------------------------------

#### Z-score normalization  --------------------------------------------------------------

# first transforming to wide format

median_epos_wide <- median_epos %>% pivot_wider(names_from = Genus, values_from = Abundance)
median_whole_wide <- median_whole %>% pivot_wider(names_from = Genus, values_from = Abundance)

# Using scale(), as per 
# https://stackoverflow.com/questions/62020169/how-to-calculate-z-score-for-each-column-of-dataframe-in-r

# want to only do this on median abundance values
# and not on "day" column (by excluding 1st column using [,-1])
# and converting output to data frame

z_median_epos <- as.data.frame(scale(median_epos_wide[,-1]))
z_median_whole <- as.data.frame(scale(median_whole_wide[,-1]))

# Adding back in the day column to my transformed data
# see https://tibble.tidyverse.org/reference/add_column.html

z_median_epos <- z_median_epos %>% add_column(.before = 1, day = median_epos_wide$day)
z_median_whole <- z_median_whole %>% add_column(.before = 1, day = median_whole_wide$day)

### Heatmap plotting ------------------------------------------------------------------------

#### ggplot heatmap -------------------------------------------------------------------------

# heatmap from littlemissdata

# selecting to only keep differentially abundant (DA) taxa

# The species marked as DA in Epos fraction are:
# Akkermansia & Erysipelatoclostridium

z_median_epos_daa <- z_median_epos %>% 
  select(c("day", "Akkermansia", "Erysipelatoclostridium"))

# to plot EdU abd of taxa that were DA in whole
# but not in EdU fraction

z_median_whole_no_epos_daa <- z_median_epos %>% 
  select(c("day", "Turicibacter","[Eubacterium]_fissicatena_group", 
           "Lachnoclostridium", "Roseburia", "Blautia"))

# The species marked as DA in the Whole community are:
# Erysipelatoclostridium, Akkermansia, Turicibacter, 
# Eubacterium fissicatena group; Lachnoclostridium, Roseburia, Blautia

z_median_whole_daa <- z_median_whole %>% 
  select(c("day", "Erysipelatoclostridium", "Akkermansia", "Turicibacter",
                    "[Eubacterium]_fissicatena_group", "Lachnoclostridium",
                    "Roseburia", "Blautia"))

# DA species found in Whole that are also found in Epos
z_median_whole_ae <- z_median_whole %>% 
  select(c("day", "Akkermansia", "Erysipelatoclostridium"))

# combining DA taxa from epos & whole
# first renaming columns
# so when i combine epos & whole
# they will be distinct

z_median_epos_ae <- rename(z_median_epos_daa, Akkermansia.epos = Akkermansia,
                           Erysipelatoclostridium.epos = Erysipelatoclostridium)

z_median_whole_ae <- rename(z_median_whole_ae, Akkermansia.whole = Akkermansia,
                            Erysipelatoclostridium.whole = Erysipelatoclostridium)

z_median_epos_whole_ae <- inner_join(z_median_epos_ae, z_median_whole_ae, by = "day")

# combining DA taxa from whole community
# not in edu fraction

# first renaming columns
# so when i combine, they will be distinct

z_median_whole_no_epos_daa_renamed <- 
  rename(z_median_whole_no_epos_daa, Turicibacter.epos = Turicibacter,
         `[Eubacterium]_fissicatena_group.epos` = `[Eubacterium]_fissicatena_group`,
         Lachnoclostridium.epos = Lachnoclostridium, Roseburia.epos = Roseburia,
         Blautia.epos = Blautia)

z_median_whole_daa_renamed <- 
  rename(z_median_whole_daa,
         Turicibacter.whole = Turicibacter,
         `[Eubacterium]_fissicatena_group.whole` = `[Eubacterium]_fissicatena_group`,
         Lachnoclostridium.whole = Lachnoclostridium, 
         Roseburia.whole = Roseburia, Blautia.whole = Blautia)

z_median_whole_epos_daa <- inner_join(z_median_whole_no_epos_daa_renamed, 
                                      z_median_whole_daa_renamed, by = "day")

# get rid of Erysipelatoclostridium & Akkermansia

z_median_whole_epos_daa <- select(z_median_whole_epos_daa, 
                                  !c("Erysipelatoclostridium", "Akkermansia"))

## change data from wide format to long format
## see https://environmentalcomputing.net/data-manipulation/reshaping-data/

z_median_epos_daa.long <- gather(z_median_epos_daa, Genus, Abundance, 2:3)
z_median_whole_daa.long <- gather(z_median_whole_daa, Genus, Abundance, 2:8)

z_median_whole_no_epos_daa.long <- gather(z_median_whole_no_epos_daa, Genus, Abundance, 2:6)

z_median_epos_whole_ae.long <- gather(z_median_epos_whole_ae, Genus, Abundance, 2:5)

z_median_whole_epos_daa.long <- gather(z_median_whole_epos_daa, Genus, Abundance, 2:11)

##### Z-transformed data --------------------------------------------------------------------

# see 
# https://vwo.com/blog/heatmap-colors
# and
# https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
# and also ?scale_fill_fermenter()
# for good information on color palettes
# for diverging data (which is what I have)

###### Labels - Epos ----------------------------------------------------------

z_median_epos_daa.long$Genus <- factor(z_median_epos_daa.long$Genus, 
                                  levels = c(
                                    "Erysipelatoclostridium", 
                                    "Akkermansia"))

taxa_labs_epos <- c(
  (expression(paste(bolditalic("Erysipelatoclostridium")))),
  (expression(paste(bolditalic("Akkermansia")))))

###### Plot - Epos ---------------------------------------------------------

# Epos
z_median_epos_daa.long %>% ggplot(aes(as.factor(day), Genus)) +
  geom_tile(aes(fill=Abundance), colour = "white") +
  #scale_fill_gradient() +
  #scale_fill_fermenter(type = "seq", palette = 1, direction = 1) + 
  scale_fill_gradient2(high = "#b2182b",
                       mid = "#f7f7f7",
                       low = "#2166ac")+
  theme_bw() + theme_minimal() +
  labs(title = "", fill = "Z-score",
       x = "Day", y = "") +
  scale_y_discrete(labels = taxa_labs_epos) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

###### Labels - Whole ----------------------------------------------------------

z_median_whole_daa.long$Genus <- factor(z_median_whole_daa.long$Genus, 
                                   levels = c(
                                     # end
                                     "Erysipelatoclostridium",
                                     "Turicibacter",
                                     
                                     # mid
                                     
                                     "Akkermansia",
                                     "[Eubacterium]_fissicatena_group",
                                     
                                     # beg
                                     "Lachnoclostridium",
                                     "Blautia",
                                     "Roseburia"
                                     
                                   ))

taxa_labs_whole <- c(
  # end
  (expression(paste(bolditalic("Erysipelatoclostridium")))),
  (expression(paste(italic("Turicibacter")))),
  (expression(paste(bolditalic("Akkermansia")))),
  (expression(paste(italic("Eubacterium fissicatena")))),
  (expression(paste(italic("Lachnoclostridium")))),
  (expression(paste(italic("Blautia")))),
  (expression(paste(italic("Roseburia"))))
  # beginning
  )
  
###### Plot - Whole ----------------------------------------------------------

# Whole
z_median_whole_daa.long %>% ggplot(aes(as.factor(day), Genus)) +
  geom_tile(aes(fill=Abundance), colour = "white") +
  #scale_fill_gradient(low = "white", high = "steelblue") +
  scale_fill_gradient2(high = "#b2182b",
                       mid = "#f7f7f7",
                       low = "#2166ac")+
  theme_bw() + theme_minimal() +
  scale_y_discrete(labels = taxa_labs_whole) +
  labs(title = "", fill = "Z-score",
       x = "Day", y = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

###### Labels - Epos from Whole ----------------------------------------------------------

z_median_whole_no_epos_daa.long$Genus <- factor(z_median_whole_no_epos_daa.long$Genus, 
                                        levels = c(
                                          # end
                                          "Turicibacter",
                                          
                                          # middle
                                          
                                          "[Eubacterium]_fissicatena_group",
                                          
                                          # beginning
                                          "Blautia",
                                          "Lachnoclostridium",
                                          "Roseburia"
                                          
                                        ))

taxa_labs_whole <- c(
  # end
  (expression(paste(italic("Turicibacter")))),
  (expression(paste(italic("Eubacterium fissicatena")))),
  (expression(paste(italic("Blautia")))),
  (expression(paste(italic("Lachnoclostridium")))),
  (expression(paste(italic("Roseburia"))))
  # beginning
)

###### Plot - Epos from Whole ----------------------------------------------------------

# Epos abundances for taxa DA in Whole sorted fraction (sf)

z_median_whole_no_epos_daa.long %>% ggplot(aes(as.factor(day), Genus)) +
  geom_tile(aes(fill=Abundance), colour = "white") +
  scale_fill_gradient2(high = "#b2182b",
                       mid = "#f7f7f7",
                       low = "#2166ac")+
  theme_bw() + theme_minimal() +
  #scale_y_discrete(labels = taxa_labs_whole) +
  labs(title = "", fill = "Z-score",
       x = "Day", y = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

###### Labels - Akker & Erysi, EdU & Whole ----------------------------------------------------------

z_median_epos_whole_ae.long$Genus <- factor(z_median_epos_whole_ae.long$Genus, 
                                                levels = c(
                                                  # end
                                                  "Erysipelatoclostridium.whole",
                                                  
                                                  # middle
                                                  
                                                  "Erysipelatoclostridium.epos",
                                                  
                                                  # beginning
                                                  "Akkermansia.whole",
                                                  "Akkermansia.epos"
                                                ))

taxa_labs_epos_whole_ae <- c(
  # end
  (expression(paste(italic("Erysipelatoclostridium"), " - Whole"))),
  (expression(paste(italic("Erysipelatoclostridium"), " - EdU+"))),
  (expression(paste(italic("Akkermansia"), " - Whole"))),
  (expression(paste(italic("Akkermansia"), " - EdU+")))
  # beg
)

###### Plot - Akker & Erysi, EdU & Whole ----------------------------------------------------------

z_median_epos_whole_ae.long %>% ggplot(aes(as.factor(day), Genus)) +
  geom_tile(aes(fill=Abundance), colour = "white") +
  scale_fill_gradient2(high = "#b2182b",
                       mid = "#f7f7f7",
                       low = "#2166ac")+
  theme_bw() + theme_minimal() +
  scale_y_discrete(labels = taxa_labs_epos_whole_ae) +
  labs(title = "", fill = "Z-score",
       x = "Day", y = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

###### Labels - Epos from Whole + Whole ---------------------------------------------------------

z_median_whole_epos_daa.long$Genus <- factor(z_median_whole_epos_daa.long$Genus, 
                                                levels = c(
                                                  # end
                                                  "Turicibacter.whole",
                                                  "Turicibacter.epos",
                                                  
                                                  # middle
                                                  "[Eubacterium]_fissicatena_group.whole",
                                                  "[Eubacterium]_fissicatena_group.epos",
                                                  
                                                  # beginning
                                                  "Blautia.whole",
                                                  "Blautia.epos",
                                                  "Lachnoclostridium.whole",
                                                  "Lachnoclostridium.epos",
                                                  "Roseburia.whole",
                                                  "Roseburia.epos"
                                                  
                                                ))

taxa_labs_whole_epos <- c(
  # end
  (expression(paste(italic("Turicibacter"), " - Whole"))),
  (expression(paste(italic("Turicibacter"), " - EdU+"))),
  (expression(paste(italic("Eubacterium fissicatena"), " - Whole"))),
  (expression(paste(italic("Eubacterium fissicatena"), " - EdU+"))),
  (expression(paste(italic("Blautia"), " - Whole"))),
  (expression(paste(italic("Blautia"), " - EdU+"))),
  (expression(paste(italic("Lachnoclostridium"), " - Whole"))),
  (expression(paste(italic("Lachnoclostridium"), " - EdU+"))),
  (expression(paste(italic("Roseburia"), " - Whole"))),
  (expression(paste(italic("Roseburia"), " - EdU+")))
  # beginning
)

###### Plot - Epos from Whole + Whole ---------------------------------------------------------

z_median_whole_epos_daa.long %>% ggplot(aes(as.factor(day), Genus)) +
  geom_tile(aes(fill=Abundance), colour = "white") +
  scale_fill_gradient2(high = "#b2182b",
                       mid = "#f7f7f7",
                       low = "#2166ac")+
  theme_bw() + theme_minimal() +
  scale_y_discrete(labels = taxa_labs_whole_epos) +
  labs(title = "", fill = "Z-score",
       x = "Day", y = "") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))
