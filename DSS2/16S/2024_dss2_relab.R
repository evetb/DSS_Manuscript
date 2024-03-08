# relative abundances & their associated bar plots for DSS #2 16S data
# using microshades (https://karstenslab.github.io/microshades/) for taxonomy plots
# following vignette: https://karstenslab.github.io/microshades/articles/microshades-GP.html

# DSS 2

# Libraries --------------------------------------------------------------------------------------------------
library(speedyseq)
library(tidyverse)
library(qiime2R)
library(microshades)
library(cowplot) # for plot_grid()

# QIIME2R -----------------------------------------------------------------------------------
dss2_phy <- qza_to_phyloseq(features = "nnn-filtered-table.qza",
                            tree = "rooted-tree.qza",
                            taxonomy = "taxonomy.qza",
                            metadata = "DSS2_metadata2.tsv")

dss2_meta <- read.table("DSS2_metadata2.tsv", header = TRUE)

# Genus abundances per cage --------------------------------------------------------

## Subsetting by cage ------------------------------------------------------------------------------

F1 <- subset_samples(dss2_phy, cage == 'F1')
F2 <- subset_samples(dss2_phy, cage == 'F2')
M1 <- subset_samples(dss2_phy, cage == 'M1')
M2 <- subset_samples(dss2_phy, cage == 'M2')

## Subsetting by cage & sf -----------------------------------------------------------

F1_df_gen_epos <- filter(F1_df_gen, sf == "EdU+")
F1_df_gen_whole <- filter(F1_df_gen, sf == "Whole")

F2_df_gen_epos <- filter(F2_df_gen, sf == "EdU+")
F2_df_gen_whole <- filter(F2_df_gen, sf == "Whole")

M1_df_gen_epos <- filter(M1_df_gen, sf == "EdU+")
M1_df_gen_whole <- filter(M1_df_gen, sf == "Whole")

M2_df_gen_epos <- filter(M2_df_gen, sf == "EdU+")
M2_df_gen_whole <- filter(M2_df_gen, sf == "Whole")

# getting rid of the one day that got contaminated
M2_df_gen_whole <- filter(M2_df_gen_whole, Sample != "M2W8")

## relative abundance overall & per cage at genus level as dataframes (df) ------------------------------------------

dss2_df_gen <- tax_glom(dss2_phy, "Genus", NArm = FALSE) %>% 
  transform_sample_counts(., function(x) x / sum(x)) %>% 
  psmelt()

F1_df_gen <- tax_glom(F1, "Genus", NArm = FALSE) %>% 
  transform_sample_counts(., function(x) x / sum(x)) %>% 
  psmelt()

F2_df_gen <- tax_glom(F2, "Genus", NArm = FALSE) %>% 
  transform_sample_counts(., function(x) x / sum(x)) %>% 
  psmelt()

M1_df_gen <- tax_glom(M1, "Genus", NArm = FALSE) %>% 
  transform_sample_counts(., function(x) x / sum(x)) %>% 
  psmelt()

M2_df_gen <- tax_glom(M2, "Genus", NArm = FALSE) %>% 
  transform_sample_counts(., function(x) x / sum(x)) %>% 
  psmelt()

## Mean genera abundance per cage ---------------------------------------------------------------------

### F1 --------------------

F1_df_gen_mean <- setNames(aggregate(F1_df_gen$Abundance,list(F1_df_gen$Genus),mean), c("Genus","Mean.Abd"))
F1_df_gen_sdv <- setNames(aggregate(F1_df_gen$Abundance,list(F1_df_gen$Genus),sd), c("Genus","Sdv.Abd"))

F1_df_gen_mean_sdv <- inner_join(F1_df_gen_mean,F1_df_gen_sdv)
F1_df_gen_mean_sdv$Cage <- "F1"

### F2 --------------------
F2_df_gen_mean <- setNames(aggregate(F2_df_gen$Abundance,list(F2_df_gen$Genus),mean), c("Genus","Mean.Abd"))
F2_df_gen_sdv <- setNames(aggregate(F2_df_gen$Abundance,list(F2_df_gen$Genus),sd), c("Genus","Sdv.Abd"))

F2_df_gen_mean_sdv <- inner_join(F2_df_gen_mean,F2_df_gen_sdv)
F2_df_gen_mean_sdv$Cage <- "F2"

### M1 --------------------
M1_df_gen_mean <- setNames(aggregate(M1_df_gen$Abundance,list(M1_df_gen$Genus),mean), c("Genus","Mean.Abd"))
M1_df_gen_sdv <- setNames(aggregate(M1_df_gen$Abundance,list(M1_df_gen$Genus),sd), c("Genus","Sdv.Abd"))

M1_df_gen_mean_sdv <- inner_join(M1_df_gen_mean,M1_df_gen_sdv)
M1_df_gen_mean_sdv$Cage <- "M1"

### M2 --------------------
M2_df_gen_mean <- setNames(aggregate(M2_df_gen$Abundance,list(M2_df_gen$Genus),mean), c("Genus","Mean.Abd"))
M2_df_gen_sdv <- setNames(aggregate(M2_df_gen$Abundance,list(M2_df_gen$Genus),sd), c("Genus","Sdv.Abd"))

M2_df_gen_mean_sdv <- inner_join(M2_df_gen_mean,M2_df_gen_sdv)
M2_df_gen_mean_sdv$Cage <- "M2"

## Mean genera abundance per cage per sf ---------------------------------------------------------------------

### F1 ----------------------------------------------

#### EdU+ ------------------------------------------

F1_df_gen_epos_mean <- setNames(aggregate(F1_df_gen_epos$Abundance,list(F1_df_gen_epos$Genus),mean), 
                                c("Genus","Mean.Abd"))
F1_df_gen_epos_sdv <- setNames(aggregate(F1_df_gen_epos$Abundance,list(F1_df_gen_epos$Genus),sd), 
                               c("Genus","Sdv.Abd"))

F1_df_gen_epos_mean_sdv <- inner_join(F1_df_gen_epos_mean,F1_df_gen_epos_sdv)
F1_df_gen_epos_mean_sdv$Cage <- "F1"

#### Whole ------------------------------------------

F1_df_gen_whole_mean <- setNames(aggregate(F1_df_gen_whole$Abundance,list(F1_df_gen_whole$Genus),mean), 
                                 c("Genus","Mean.Abd"))
F1_df_gen_whole_sdv <- setNames(aggregate(F1_df_gen_whole$Abundance,list(F1_df_gen_whole$Genus),sd), 
                                c("Genus","Sdv.Abd"))

F1_df_gen_whole_mean_sdv <- inner_join(F1_df_gen_whole_mean,F1_df_gen_whole_sdv)
F1_df_gen_whole_mean_sdv$Cage <- "F1"

### F2 ----------------------------------------------

#### EdU+ ------------------------------------------

F2_df_gen_epos_mean <- setNames(aggregate(F2_df_gen_epos$Abundance,list(F2_df_gen_epos$Genus),mean), 
                                c("Genus","Mean.Abd"))
F2_df_gen_epos_sdv <- setNames(aggregate(F2_df_gen_epos$Abundance,list(F2_df_gen_epos$Genus),sd), 
                               c("Genus","Sdv.Abd"))

F2_df_gen_epos_mean_sdv <- inner_join(F2_df_gen_epos_mean,F2_df_gen_epos_sdv)
F2_df_gen_epos_mean_sdv$Cage <- "F2"

#### Whole ------------------------------------------

F2_df_gen_whole_mean <- setNames(aggregate(F2_df_gen_whole$Abundance,list(F2_df_gen_whole$Genus),mean), 
                                 c("Genus","Mean.Abd"))
F2_df_gen_whole_sdv <- setNames(aggregate(F2_df_gen_whole$Abundance,list(F2_df_gen_whole$Genus),sd), 
                                c("Genus","Sdv.Abd"))

F2_df_gen_whole_mean_sdv <- inner_join(F2_df_gen_whole_mean,F2_df_gen_whole_sdv)
F2_df_gen_whole_mean_sdv$Cage <- "F2"

### M1 ----------------------------------------------

#### EdU+ ------------------------------------------

M1_df_gen_epos_mean <- setNames(aggregate(M1_df_gen_epos$Abundance,list(M1_df_gen_epos$Genus),mean), 
                                c("Genus","Mean.Abd"))
M1_df_gen_epos_sdv <- setNames(aggregate(M1_df_gen_epos$Abundance,list(M1_df_gen_epos$Genus),sd), 
                               c("Genus","Sdv.Abd"))

M1_df_gen_epos_mean_sdv <- inner_join(M1_df_gen_epos_mean,M1_df_gen_epos_sdv)
M1_df_gen_epos_mean_sdv$Cage <- "M1"

#### Whole ------------------------------------------

M1_df_gen_whole_mean <- setNames(aggregate(M1_df_gen_whole$Abundance,list(M1_df_gen_whole$Genus),mean), 
                                 c("Genus","Mean.Abd"))
M1_df_gen_whole_sdv <- setNames(aggregate(M1_df_gen_whole$Abundance,list(M1_df_gen_whole$Genus),sd), 
                                c("Genus","Sdv.Abd"))

M1_df_gen_whole_mean_sdv <- inner_join(M1_df_gen_whole_mean,M1_df_gen_whole_sdv)
M1_df_gen_whole_mean_sdv$Cage <- "M1"

### M2 ----------------------------------------------

#### EdU+ ------------------------------------------

M2_df_gen_epos_mean <- setNames(aggregate(M2_df_gen_epos$Abundance,list(M2_df_gen_epos$Genus),mean), 
                                c("Genus","Mean.Abd"))
M2_df_gen_epos_sdv <- setNames(aggregate(M2_df_gen_epos$Abundance,list(M2_df_gen_epos$Genus),sd), 
                               c("Genus","Sdv.Abd"))

M2_df_gen_epos_mean_sdv <- inner_join(M2_df_gen_epos_mean,M2_df_gen_epos_sdv)
M2_df_gen_epos_mean_sdv$Cage <- "M2"

#### Whole ------------------------------------------

M2_df_gen_whole_mean <- setNames(aggregate(M2_df_gen_whole$Abundance,list(M2_df_gen_whole$Genus),mean), 
                                 c("Genus","Mean.Abd"))
M2_df_gen_whole_sdv <- setNames(aggregate(M2_df_gen_whole$Abundance,list(M2_df_gen_whole$Genus),sd), 
                                c("Genus","Sdv.Abd"))

M2_df_gen_whole_mean_sdv <- inner_join(M2_df_gen_whole_mean,M2_df_gen_whole_sdv)
M2_df_gen_whole_mean_sdv$Cage <- "M2"

## Top Genera per cage & sf --------------------------------------------------------

# arranging the df in descending order
# (highest values first)
# then choosing only the first 10 rows using head()
# and saving the result into a new df

# https://dplyr.tidyverse.org/reference/arrange.html

## Epos

F1_top10_genera_epos <- F1_df_gen_epos_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

F2_top10_genera_epos <- F2_df_gen_epos_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

M1_top10_genera_epos <- M1_df_gen_epos_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

M2_top10_genera_epos <- M2_df_gen_epos_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

## Whole

F1_top10_genera_whole <- F1_df_gen_whole_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

F2_top10_genera_whole <- F2_df_gen_whole_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

M1_top10_genera_whole <- M1_df_gen_whole_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

M2_top10_genera_whole <- M2_df_gen_whole_mean_sdv %>% 
  arrange(desc(Mean.Abd)) %>% 
  head(10)

# finding the commonality between all of these
# to get top 5 genera shared between all of them

# https://sparkbyexamples.com/r-programming/r-join-multiple-data-frames/?expand_article=1
# https://stackoverflow.com/questions/68702694/how-can-i-keep-data-from-one-dataframe-based-on-its-appearance-in-a-second-dataf

all_cages_top10_epos <- reduce(list(F1_top10_genera_epos, 
                                    F2_top10_genera_epos,
                                    M1_top10_genera_epos,
                                    M2_top10_genera_epos),
                               inner_join, by = "Genus")

all_cages_top10_whole <- reduce(list(F1_top10_genera_whole, 
                                     F2_top10_genera_whole,
                                     M1_top10_genera_whole,
                                     M2_top10_genera_whole),
                                inner_join, by = "Genus")

# now arranging it all by the Mean.Abd columns

all_cages_top10_epos <- all_cages_top10_epos %>% 
  arrange(desc(Mean.Abd.x), desc(Mean.Abd.y), desc(Mean.Abd.x.x), desc(Mean.Abd.y.y))

all_cages_top10_whole <- all_cages_top10_whole %>% 
  arrange(desc(Mean.Abd.x), desc(Mean.Abd.y), desc(Mean.Abd.x.x), desc(Mean.Abd.y.y))

# Microshades --------------------------------------------------------------------------------------
## Data manipulation -------------------------------------------------------------------------------------------

# Agglomerate and normalize the phyloseq object, and melt to a data frame
dss2_df <- prep_mdf(dss2_phy)

# Get rid of sample M2W8, which is contaminated
dss2_df <- filter(dss2_df, Sample != "M2W8")

# Changing phyla names, as per:
# https://datatofish.com/replace-values-dataframe-r/

dss2_df["Phylum"][dss2_df["Phylum"] == "Firmicutes"] <- "Bacillota"
dss2_df["Phylum"][dss2_df["Phylum"] == "Actinobacteriota"] <- "Actinomycetota"
dss2_df["Phylum"][dss2_df["Phylum"] == "Proteobacteria"] <- "Pseudomonadota"

# Generate a color object for the specified data
color_objs_dss2 <- 
  create_color_dfs(dss2_df,
                   selected_groups =
                     c("Verrucomicrobiota", "Pseudomonadota",
                       "Actinomycetota", "Bacteroidota",
                       "Bacillota") , cvd = TRUE)

# Extract
mdf_dss2 <- color_objs_dss2$mdf
cdf_dss2 <- color_objs_dss2$cdf

# turning day into a factor & specifying levels
mdf_dss2$day <- factor(mdf_dss2$day, 
                       levels=c("1", "5", "7", "8", "11",
                                "12", "13", "14", "15", 
                                "18", "22", "26", "29"))

# turning hs into a factor & specifying levels
mdf_dss2$hs <- factor(mdf_dss2$hs,
                                levels = c("Baseline", "Pre-symptomatic", 
                                           "Symptomatic", "Recovery"))

## Data manipulation - per sf (all cages) ---------------------------------------------------

# making separate data frames per cage & per sf

dss2_df_epos <- dss2_df %>% filter(sf == "EdU+")
dss2_df_whole <- dss2_df %>% filter(sf == "Whole")

# Generate a color object for the specified data
color_objs_epos <- 
  create_color_dfs(dss2_df_epos,
                   selected_groups =
                     c("Pseudomonadota", "Verrucomicrobiota",
                       "Actinomycetota", "Bacteroidota",
                       "Bacillota") , cvd = TRUE)

color_objs_whole <- 
  create_color_dfs(dss2_df_whole,
                   selected_groups =
                     c("Pseudomonadota", "Verrucomicrobiota",
                       "Actinomycetota", "Bacteroidota",
                       "Bacillota") , cvd = TRUE)

# Extract
mdf_epos <- color_objs_epos$mdf
cdf_epos <- color_objs_epos$cdf

mdf_whole <- color_objs_whole$mdf
cdf_whole <- color_objs_whole$cdf

# turning day into a factor & specifying levels
mdf_epos$day <- factor(mdf_epos$day, 
                          levels=c("1", "5", "7", "8", "11",
                                   "12", "13", "14", "15", 
                                   "18", "22", "26", "29"))

mdf_whole$day <- factor(mdf_whole$day, 
                           levels=c("1", "5", "7", "8", "11",
                                    "12", "13", "14", "15", 
                                    "18", "22", "26", "29"))

# turning hs into a factor & specifying levels
mdf_epos$hs <- factor(mdf_epos$hs,
                         levels = c("Baseline", "Pre-symptomatic", 
                                    "Symptomatic", "Recovery"),
                         labels = c("BA", "PS", 
                                    "SY", "RE"))

mdf_whole$hs <- factor(mdf_whole$hs,
                          levels = c("Baseline", "Pre-symptomatic", 
                                     "Symptomatic", "Recovery"),
                          labels = c("BA", "PS", 
                                     "SY", "RE"))

## Plot ---------------------------

### Per sf, all cages --------------------------------------

# custom legends
epos_legend <- custom_legend(mdf_epos, cdf_epos, 
                           legend_key_size = 0.6,
                           legend_text_size = 18)
whole_legend <- custom_legend(mdf_whole, cdf_whole, 
                           legend_key_size = 0.6,
                           legend_text_size = 18)

#### Epos -----------------------------------

epos_plot <- plot_microshades(mdf_epos, cdf_epos, x = "day") + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  labs(x = "Day", y = "Relative Abundance", title = "DSS #2 - EdU+") +
  facet_grid(rows = vars(cage), 
             cols = vars(hs), 
             scales = "free", space = "free") +
  theme_bw() +
  theme(text = element_text(size = 20), 
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        panel.spacing.y = unit(1.5, "lines")) 

plot_grid(epos_plot, epos_legend, rel_widths = c(.75, .3))

#### Whole ----------------------------------

whole_plot <- plot_microshades(mdf_whole, cdf_whole, x = "day") + 
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  labs(x = "Day", y = "Relative Abundance", title = "DSS #2 - Whole") +
  facet_grid(rows = vars(cage), 
             cols = vars(hs), 
             scales = "free", space = "free") +
  theme_bw() +
  theme(text = element_text(size = 20), 
        plot.title = element_text(size = 24, hjust = 0.5),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        panel.spacing.y = unit(1.5, "lines")) 

plot_grid(whole_plot, whole_legend, rel_widths = c(.75, .3))
