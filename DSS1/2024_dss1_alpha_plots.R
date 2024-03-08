# alpha diversity plots for DSS #1

# Libraries --------------------------------------------------------------

library(tidyverse)

# Data manipulation --------------------------------------------------------

# .tsv files are from QIIME2 output
# for values of various alpha diversity metrics

OTUs <- read_tsv("dss1_obvs_otus_day.tsv")
shannon <- read_tsv("dss1_shannon_day.tsv")
simpson <- read_tsv("dss1_simpson_day.tsv")
simpsonE <- read_tsv("dss1_simpsonE_day.tsv")

# compile together, as per: 
# https://sparkbyexamples.com/r-programming/r-join-multiple-data-frames/

list_alpha <- list(OTUs, shannon, simpson, simpsonE)
all_alpha <- list_alpha %>% 
  reduce(inner_join,
         by = c("id", "sex", "cage", 
                `health-state`, "day", "extraction-date",
                "collection-date"))

# getting rid of first row
# (descriptor for QIIME2)

all_alpha <- all_alpha[-1,]

# fixing health state
# first making day numeric

all_alpha$day <- as.numeric(all_alpha$day)

all_alpha$`health-state` <- 
  ifelse(all_alpha$day >= 1 & all_alpha$day <=8, "Baseline",
         ifelse(all_alpha$day >= 9 & all_alpha$day <= 11, "Pre-symptomatic",
                ifelse(all_alpha$day >= 12 & all_alpha$day <= 15, "Symptomatic", "Recovery")))

# making health state a factor & specifying order of levels
# and giving labels for the levels for graphing purposes

all_alpha$`health-state` <- factor(all_alpha$`health-state`,
                                   levels = c("Baseline", "Pre-symptomatic",
                                              "Symptomatic", "Recovery"),
                                   labels = c("BA", "PS", "SY", "RE"))

# Plots -----------------------------------------------------------------------

# making all alpha diversity metrics into numeric

all_alpha$observed_features <- as.numeric(all_alpha$observed_features)
all_alpha$shannon_entropy <- as.numeric(all_alpha$shannon_entropy)
all_alpha$simpson <- as.numeric(all_alpha$simpson)
all_alpha$simpson_e <- as.numeric(all_alpha$simpson_e)

# Boxplots ---------------------------------------------------------------------------

# per health state

## F1 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1")), 
                     mapping = aes(x=`health-state`, y=observed_features, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "DSS #1 - Cage F1 OTUs",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1")), 
                     mapping = aes(x=`health-state`, y=shannon_entropy, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "DSS #1 - Cage F1 Shannon",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1")), 
                     mapping = aes(x=`health-state`, y=simpson_e, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "DSS #1 - Cage F1 Simpson E",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1")), 
                     mapping = aes(x=`health-state`, y=simpson, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "DSS #1 - Cage F1 Simpson",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## F2 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2")), 
                     mapping = aes(x=`health-state`, y=observed_features, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "DSS #1 - Cage F2 OTUs",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2")), 
                     mapping = aes(x=`health-state`, y=shannon_entropy, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "DSS #1 - Cage F2 Shannon",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2")), 
                     mapping = aes(x=`health-state`, y=simpson_e, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "DSS #1 - Cage F2 Simpson E",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2")), 
                     mapping = aes(x=`health-state`, y=simpson, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "DSS #1 - Cage F2 Simpson",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## M1 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1")), 
                     mapping = aes(x=`health-state`, y=observed_features, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "DSS #1 - Cage M1 OTUs",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1")), 
                     mapping = aes(x=`health-state`, y=shannon_entropy, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "DSS #1 - Cage M1 Shannon",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1")), 
                     mapping = aes(x=`health-state`, y=simpson_e, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "DSS #1 - Cage M1 Simpson E",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1")), 
                     mapping = aes(x=`health-state`, y=simpson, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "DSS #1 - Cage M1 Simpson",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## M2 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2")), 
                     mapping = aes(x=`health-state`, y=observed_features, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "DSS #1 - Cage M2 OTUs",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2")), 
                     mapping = aes(x=`health-state`, y=shannon_entropy, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "DSS #1 - Cage M2 Shannon",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2")), 
                     mapping = aes(x=`health-state`, y=simpson_e, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "DSS #1 - Cage M2 Simpson E",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2")), 
                     mapping = aes(x=`health-state`, y=simpson, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "DSS #1 - Cage M2 Simpson",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## Females --------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female")), 
                     mapping = aes(x=`health-state`, y=observed_features, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "DSS #1 - Female OTUs",
       shape = "Cage",
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female")), 
                     mapping = aes(x=`health-state`, y=shannon_entropy, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "DSS #1 - Female Shannon",
       shape = "Cage",
       fill = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female")), 
                     mapping = aes(x=`health-state`, y=simpson_e, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "DSS #1 - Female Simpson E",
       shape = "Cage",
       fill = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female")), 
                     mapping = aes(x=`health-state`, y=simpson, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "DSS #1 - Female Simpson",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Males ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male")), 
                     mapping = aes(x=`health-state`, y=observed_features, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "DSS #1 - Male OTUs",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male")), 
                     mapping = aes(x=`health-state`, y=shannon_entropy, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "DSS #1 - Male Shannon",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male")), 
                     mapping = aes(x=`health-state`, y=simpson_e, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "DSS #1 - Male Simpson E",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male")), 
                     mapping = aes(x=`health-state`, y=simpson, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "DSS #1 - Male Simpson",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## All cages ------------------------------------------------------------------------------

all_alpha %>% ggplot(., 
                     mapping = aes(x=`health-state`, y=observed_features, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Observed Features", 
       #title = "DSS #1 - All OTUs",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

all_alpha %>% ggplot(., 
                     mapping = aes(x=`health-state`, y=shannon_entropy, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Shannon Index", 
       #title = "DSS #1 - All Shannon",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

# also Friedman version - with medians per cage per hs

all.cages_shannon %>% ggplot(., 
                     mapping = aes(x=health.state, y=Med.shannon, fill = health.state)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(x = "Health State", 
       y = "Shannon", 
       title = "DSS #1 - Median Shannon",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(., 
                     mapping = aes(x=`health-state`, y=simpson_e, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Simpson Evenness", 
       #title = "DSS #1 - All Simpson E",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

all_alpha %>% ggplot(., 
                     mapping = aes(x=`health-state`, y=simpson, fill = `health-state`)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "DSS #1 - All Simpson",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))
