# making line graphs and box plots of DSS #2 
# flow data - median %EdU+ 

# Libraries ----------------------------

library(tidyverse)
library(stringr)
library(lmerTest)
library(lspline)
library(segmented)

# Data manipulation ----------------------

## import ------------

edu <- read_csv("prop_epos.csv")

## turning data from wide to long format ----------------------------

edu.long <- pivot_longer(edu, !day, names_to = "cage", values_to = "prop_edu")

## changing cage names (removing numbers past period)------------------

edu.long$cage <- str_remove(edu.long$cage,regex("\\.[0-9]",dotall = T))

## median & SD & MAD ----------------------------------

# first making data frame for each cage, 
# since I want to get the median %EdU+  
# of each replicate in each cage per day
# https://dplyr.tidyverse.org/reference/filter.html

F1_edu <- filter(edu.long, cage == "F1")
F2_edu <- filter(edu.long, cage == "F2")
M1_edu <- filter(edu.long, cage == "M1")
M2_edu <- filter(edu.long, cage == "M2")

# getting the median edu 
# & sd per replicate per day w/in each cage

# grouping days by using list()
# changing names of columns using setNames(), as per
# https://statisticsglobe.com/set-column-names-within-aggregate-function-in-r

F1_median_edu <- setNames(aggregate(F1_edu$prop_edu,list(F1_edu$day),median), c("Day","Med.edu"))
F1_mad_edu <- setNames(aggregate(F1_edu$prop_edu,list(F1_edu$day),mad), c("Day","MAD.edu"))

F1_sd_edu <- setNames(aggregate(F1_edu$prop_edu,list(F1_edu$day),sd), c("Day","SD.edu"))

F2_median_edu <- setNames(aggregate(F2_edu$prop_edu,list(F2_edu$day),median), c("Day", "Med.edu"))
F2_mad_edu <- setNames(aggregate(F2_edu$prop_edu,list(F2_edu$day),mad), c("Day", "MAD.edu"))

F2_sd_edu <- setNames(aggregate(F2_edu$prop_edu,list(F2_edu$day),sd), c("Day","SD.edu"))

M1_median_edu <- setNames(aggregate(M1_edu$prop_edu,list(M1_edu$day),median), c("Day","Med.edu"))
M1_mad_edu <- setNames(aggregate(M1_edu$prop_edu,list(M1_edu$day),mad), c("Day","MAD.edu"))

M1_sd_edu <- setNames(aggregate(M1_edu$prop_edu,list(M1_edu$day),sd), c("Day","SD.edu"))

M2_median_edu <- setNames(aggregate(M2_edu$prop_edu,list(M2_edu$day),median), c("Day","Med.edu"))
M2_mad_edu <- setNames(aggregate(M2_edu$prop_edu,list(M2_edu$day),mad), c("Day","MAD.edu"))

M2_sd_edu <- setNames(aggregate(M2_edu$prop_edu,list(M2_edu$day),sd), c("Day","SD.edu"))

# creating data frames with the values of median & sd
# using inner_join()

F1_med_mad_edu <- inner_join(F1_median_edu,F1_mad_edu)
F1_med_mad_edu$Cage <- "F1"

F2_med_mad_edu <- inner_join(F2_median_edu,F2_mad_edu)
F2_med_mad_edu$Cage <- "F2"

M1_med_mad_edu <- inner_join(M1_median_edu,M1_mad_edu)
M1_med_mad_edu$Cage <- "M1"

M2_med_mad_edu <- inner_join(M2_median_edu,M2_mad_edu)
M2_med_mad_edu$Cage <- "M2"

# Combining all above medians & SDs together into one data frame

all.cages_edu_med_mad <- rbind(F1_med_mad_edu,F2_med_mad_edu,M1_med_mad_edu,M2_med_mad_edu)

# median & mad per health state 
# for all cages

all.cages_edu_med_hs <- setNames(aggregate(all.cages_edu_med_mad$Med.edu,
                                               list(all.cages_edu_med_mad$hs),median), 
                                     c("hs","Med.edu"))

all.cages_edu_mad_hs <- setNames(aggregate(all.cages_edu_med_mad$Med.edu,
                                               list(all.cages_edu_med_mad$hs),mad), 
                                     c("hs","MAD.edu"))


all.cages_edu_med_mad_hs <- inner_join(all.cages_edu_med_hs, all.cages_edu_mad_hs)

## Variability ----------------------------------------------

# quantifying variability in median 

all.cages_edu_med_mad$Day <- factor(all.cages_edu_med_mad$Day,
                                    levels = c("1", "5", "7", "8",
                                               "11", "12", "13", "14",
                                               "15", "18", "26", "29",
                                               "33"))

### within cages ------------------

all.cages_edu_med_mad %>% ggplot(aes(x = Day, y = MAD.edu, 
                                     group = Cage, colour = Cage)) +
  geom_point(size = 4) + 
  geom_line(linewidth = 2) +
  #geom_errorbar(aes(ymin=Med.edu-MAD.edu, ymax=Med.edu+MAD.edu),linewidth = 2) +
  labs(x = "Day", y = "MAD", title = "MAD for Intra-Median %EdU+") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))

### between cages ------------------------------

# MAD for the median each day for all cages
# run the below to recreate Days_MAD if turning
# Day from factor back into numeric

Days_MAD <- setNames(aggregate(all.cages_edu_med_mad$Med.edu,
                             list(all.cages_edu_med_mad$Day), mad), c("Day","InterMAD.edu"))

Days_MAD$Day <- factor(Days_MAD$Day,
                                    levels = c("1", "5", "7", "8",
                                               "11", "12", "13", "14",
                                               "15", "18", "26", "29",
                                               "33"))

# adding column so I can group them after using factor

Days_MAD$Cage <- "all"

# need to add health states

Days_MAD$hs <- cut(Days_MAD$Day, c(-Inf, 7, 11, 15, Inf), 
                      c("Baseline", "Pre-symptomatic",
                        "Symptomatic", "Recovery"))

# ordering health states

Days_MAD$hs <- factor(Days_MAD$hs, 
                      levels = c("Baseline", "Pre-symptomatic",
                                 "Symptomatic", "Recovery"), 
                      labels = c("BA", "PS", "SY", "RE"))

#### plotting - line plot -----------------------
Days_MAD %>% ggplot(aes(x = Day, y = InterMAD.edu, group = Cage)) +
  geom_point(size = 4) + 
  geom_line(linewidth = 2) +
  #geom_errorbar(aes(ymin=Med.edu-MAD.edu, ymax=Med.edu+MAD.edu),linewidth = 2) +
  labs(x = "Day", y = "MAD", title = "Inter-cage MAD of %EdU+ cells") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 34, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))

#### plotting - box plot per hs -------------------------
Days_MAD %>% ggplot(aes(x = hs, y = InterMAD.edu, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "MAD", title = "Inter-cage MAD of %EdU+ cells") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  #facet_grid(. ~ hs, scales = "free", space = "free") + 
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 34, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.position = "none")

all.cages_edu %>% ggplot(aes(x = hs, y = Med.edu, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Health State", 
       y = "Median %EdU+ Cells", fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  #facet_grid(. ~ hs, scales = "free", space = "free") + 
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

# Plots ------------------------------------

## Frequency histogram -------------------------

# https://www.sfu.ca/~mjbrydon/tutorials/BAinR/visualize.html

# for values of Med.edu
hist(all.cages_edu_med_mad$Med.edu) 

## Line plots ---------------

all.cages_edu %>% ggplot(aes(x = Day, y = Med.edu, colour = Cage)) +
  geom_point(size = 4) + 
  geom_line(linewidth = 2) +
  geom_errorbar(aes(ymin=Med.edu-SD.edu, ymax=Med.edu+SD.edu), 
                linewidth = 2) +
  labs(x = "Day", y = "Median %EdU+", title = "") +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))

# day as a factor

all.cages_edu_med_mad$Day <- factor(all.cages_edu_med_mad$Day,
                                    levels = c("1", "5", "7", "8",
                                               "11", "12", "13", "14",
                                               "15", "18", "26", "29",
                                               "33"))

all.cages_edu_med_mad %>% ggplot(aes(x = Day, y = Med.edu, 
                                     group = Cage, colour = Cage)) +
  geom_point(size = 4) + 
  geom_line(linewidth = 2) +
  geom_errorbar(aes(ymin=Med.edu-MAD.edu, ymax=Med.edu+MAD.edu), 
                linewidth = 2) +
  labs(x = "Day", y = "Median %EdU+", title = "") + #DSS #2 - Median %EdU+
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))

## Box plots --------------------

### All data -----------------------

# need to add health states in first
# from https://stackoverflow.com/questions/50988447/add-column-with-values-depending-on-another-column-to-a-dataframe

all.cages_edu$hs <- ifelse(all.cages_edu$Day >= 1 & all.cages_edu$Day <=7, "Baseline",
                           ifelse(all.cages_edu$Day >= 8 & all.cages_edu$Day <= 11, "Pre-symptomatic",
                           ifelse(all.cages_edu$Day >= 12 & all.cages_edu$Day <= 15, "Symptomatic", "Recovery")))

# using cut() might be less tedious, e.g.
# df$c4 <- cut(df$c2, c(-Inf,4,9,Inf), c("low", "medium", "high"))

# converting hs into a factor and specifying levels
all.cages_edu$hs <- factor(all.cages_edu$hs, levels = c("Baseline", "Pre-symptomatic",
                                                        "Symptomatic", "Recovery"),
                           labels = c("BA", "PS", "SY", "RE"))

# Median EdU

all.cages_edu %>% ggplot(aes(x = hs, y = Med.edu, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Health State", 
       y = "Median %EdU+ Cells", fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  #facet_grid(. ~ hs, scales = "free", space = "free") + 
  theme_minimal() + theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.position = "none")

### Median per cage per hs -----------------------

# Median EdU

all.cages_edu_hs %>% ggplot(aes(x = hs, y = Med.edu.hs, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Health State", 
       y = "Median %EdU+ Cells", fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  #facet_grid(. ~ hs, scales = "free", space = "free") + 
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
    axis.text = element_text(size = 16),
    axis.title=element_text(size=24), legend.position = "none")

### Mean per cage per hs ---------------------------------

all.cages_mean_edu_hs %>% ggplot(aes(x = hs, y = Mean.edu.hs, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Health State", 
       y = "Median %EdU+ Cells", fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  #facet_grid(. ~ hs, scales = "free", space = "free") + 
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")