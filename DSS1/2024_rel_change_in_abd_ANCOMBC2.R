# relative change in abundances 
# for differentially abundant taxa
# from ANCOMBC2 results

# Libraries ------------------------------------------------------------

library(tidyverse)
library(phyloseq)

# Data Manipulation ----------------------------------------------------------------

## Phyloseq object to genus level ---------------------------------------

dss1_phy_fix

dss1_phy_genus <- tax_glom(dss1_phy_fix, taxrank = "Genus")

dss1_phy_genus.df <- psmelt(dss1_phy_genus)

## Merge taxonomy names ---------------------------------------------------------
# into single Classification column
# https://www.marsja.se/how-to-concatenate-two-columns-or-more-in-r-stringr-tidyr/

dss1_phy_genus.df_merged <- dss1_phy_genus.df %>% 
  unite("Classification", Phylum:Genus, remove = FALSE)

## Pivot_wider ----------------------------------------------------------------

# make new dataframe (df) with only 
# Sample, Classification, and Abundance
# and the relevant metadata
# keeping Genus for easier plotting

dss1_phy_genus.df_classification <- dss1_phy_genus.df_merged %>% 
  select(Sample, Abundance, Classification, Genus, day, health.state, cage, sex)

# specifying order of health states
dss1_phy_genus.df_classification$health.state <- 
  factor(dss1_phy_genus.df_classification$health.state,
         levels = c("Baseline", "Pre-symptomatic",
                    "Symptomatic", "Recovery"))

# then make df with Sample as first column,
# Classification as second column,
# and Abundance as values

dss1_phy_genus.df_classification_wider <- dss1_phy_genus.df_classification %>% 
  pivot_wider(names_from = "Classification", values_from = "Abundance")

# Plotting ------------------------------------------------------------------

# Here is the list of taxa I need to plot

# *Erysipelatoclostridium
# *Lactobacillus
# *[Eubacterium]_fissicatena_group
# *Peptococcaceae Family
# *Clostridia_vadinBB60_group
# *Clostridia_UCG-014
# *Parasutterella
# *Muribaculaceae
# *Oscillibacter
# *Lachnospiraceae_UCG-008
# *RF39 (Bacilli)
# *Faecalibaculum
# *Lachnospiraceae_FCS020_group
# *Incertae_Sedis
# *Ruminococcaceae Family
# *UCG-005 (Oscillospiraceae)
# *[Eubacterium]_coprostanoligenes_group
# *Anaeroplasma

# A good way to get Genus names easily from the df:
# (in alphabetical order)
sort(unique(dss1_phy_genus.df_classification$Genus))

## Heat maps -----------------------------------------------------------------

### Data preparation --------------------------------------------------------------

# I already have a dataframe at genus level

dss1_phy_genus.df_classification

# getting the median of each bacterium per day, grouping days by using list()
# https://stackoverflow.com/questions/21982987/mean-per-group-in-a-data-frame

# grouping by Genus and by day
# https://r-coder.com/aggregate-r/
median_abd <- aggregate(dss1_phy_genus.df_classification$Abundance,
                         list(dss1_phy_genus.df_classification$Genus, 
                              dss1_phy_genus.df_classification$day),median)

# changing column name from Group.1 to Genus
# and from Group.2 to day
# and from x to Abundance
median_abd <- median_abd %>% 
  rename("Genus" = "Group.1") %>%
  rename("day" = "Group.2") %>% 
  rename("Abundance" = "x")

### Prep - data transformation -----------------------------------------------------------

#### Z-score normalization  --------------------------------------------------------------

# first transforming to wide format

median_abd_wide <- median_abd %>% pivot_wider(names_from = Genus, values_from = Abundance)

# Using scale(), as per 
# https://stackoverflow.com/questions/62020169/how-to-calculate-z-score-for-each-column-of-dataframe-in-r

# want to only do this on median abundance values
# and not on "day" column (by excluding 1st column using [,-1])
# and converting output to data frame

z_median_abd <- as.data.frame(scale(median_abd_wide[,-1]))

# Adding back in the day column to my transformed data
# see https://tibble.tidyverse.org/reference/add_column.html

z_median_abd <- z_median_abd %>% add_column(.before = 1, day = median_abd_wide$day)

### Heatmap plotting ------------------------------------------------------------------------

#### ggplot heatmap -------------------------------------------------------------------------

# heatmap from littlemissdata

# selecting to only keep differentially abundant taxa

# The species marked as DA are:
# Erysipelatoclostridium
# Lactobacillus
# [Eubacterium]_fissicatena_group
# Peptococcaceae Family
# Clostridia_vadinBB60_group
# Clostridia_UCG-014
# Parasutterella
# Muribaculaceae
# Oscillibacter
# Lachnospiraceae_UCG-008
# RF39 (Bacilli)
# Faecalibaculum
# Lachnospiraceae_FCS020_group
# Incertae_Sedis
# Ruminococcaceae Family
# UCG-005 (Oscillospiraceae)
# [Eubacterium]_coprostanoligenes_group
# Anaeroplasma

z_median_abd_daa <- z_median_abd %>% 
  select(c("day", "Erysipelatoclostridium", "Lactobacillus", 
           "[Eubacterium]_fissicatena_group", "Peptococcaceae Family",
           "Clostridia_vadinBB60_group", "Clostridia_UCG-014", 
           "Parasutterella", "Muribaculaceae", "Oscillibacter", 
           "Lachnospiraceae_UCG-008", "RF39", "Faecalibaculum",
            "Lachnospiraceae_FCS020_group", "Incertae_Sedis",
           "Ruminococcaceae Family", "UCG-005", 
           "[Eubacterium]_coprostanoligenes_group",
           "Anaeroplasma"))

## change data from wide format back to long format
## see https://environmentalcomputing.net/data-manipulation/reshaping-data/

z_median_abd_daa.long <- gather(z_median_abd_daa, Genus, Abundance, 2:19)

##### Z-transformed data --------------------------------------------------------------------

# see 
# https://vwo.com/blog/heatmap-colors
# and
# https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
# and also ?scale_fill_fermenter()
# for good information on color palettes
# for diverging data (which is what I have)

###### Labels ----------------------------------------------------------

z_median_abd_daa.long$Genus <- factor(z_median_abd_daa.long$Genus, 
                                       levels = c(
                                         # bottom of graph
                                         "Faecalibaculum", 
                                         "Lachnospiraceae_UCG-008",
                                         "Lactobacillus", 
                                         "Parasutterella",
                                         "[Eubacterium]_fissicatena_group", 
                                         "Erysipelatoclostridium",
                                         "Lachnospiraceae_FCS020_group",
                                         "Clostridia_UCG-014", 
                                         "Clostridia_vadinBB60_group", 
                                         "Oscillibacter", 
                                         "Muribaculaceae", 
                                         "Peptococcaceae Family",
                                         "Anaeroplasma", 
                                         "Incertae_Sedis",
                                         "[Eubacterium]_coprostanoligenes_group",
                                         "Ruminococcaceae Family", 
                                         "UCG-005", 
                                         "RF39"
                                         # top of graph
                                         ))

taxa_labs_abd <- c(
  # bottom of graph
  (expression(paste(italic("Faecalibaculum")))),
  (expression(paste(italic("Lachnospiraceae"), " UCG-008"))),
  (expression(paste(italic("Lactobacillus")))),
  (expression(paste(italic("Parasutterella")))),
  (expression(paste(italic("Eubacterium fissicatena"), " group"))),
  (expression(paste(italic("Erysipelatoclostridium")))),
  (expression(paste(italic("Lachnospiraceae"), " FCS020_group"))),
  (expression(paste(italic("Clostridia"), " UCG-014"))),
  (expression(paste(italic("Clostridia vadin"), " BB60_group"))),
  (expression(paste(italic("Oscillibacter")))),
  (expression(paste(italic("Muribaculaceae")))),
  (expression(paste(italic("Peptococcaceae")))),
  (expression(paste(italic("Anaeroplasma")))),
  (expression(paste(italic("Incertae sedis")))),
  (expression(paste(italic("Eubacterium coprostanoligenes"), " group"))),
  (expression(paste(italic("Ruminococcaceae")))),
  (expression(paste(italic("Oscillospiraceae"), " UCG-005"))),
  (expression(paste(italic("Bacilli"), " RF39")))
  # top of graph
  )

###### Plot ---------------------------------------------------------

z_median_abd_daa.long %>% ggplot(aes(as.factor(day), Genus)) +
  geom_tile(aes(fill=Abundance), colour = "white") +
  scale_fill_gradient2(high = "#b2182b",
                       mid = "#f7f7f7",
                       low = "#2166ac")+
  theme_bw() + theme_minimal() +
  labs(title = "", fill = "Z-score",
       x = "Day", y = "") +
  scale_y_discrete(labels = taxa_labs_abd) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))
