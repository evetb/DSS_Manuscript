# statistics for edu data from DSS #2

# wilcoxon sign-ranked test with post-hoc
# adjustment for multiple comparisons

# Libraries ----------------------------

library(tidyverse)
library(rstatix) # for pipe-friendly friedman test

# Data manipulation ----------------------

## Single-computed Medians --------------------------------

all.cages_edu

# turning cage into a factor
all.cages_edu$Cage <- factor(all.cages_edu$Cage, 
                                levels = c("F1", "F2", "M1", "M2"))

## Double-computed Medians -----------------------------------

# using medians calculated in 2024_dss2_edu_graphs.R script

# need to compute the median EdU values
# for each time period, per cage

# adding hs to the data frames

# Baseline: 1, 5, 7
# Pre-symp: 8, 11
# Symp: 12, 13, 14, 15
# Rec: 18, 26, 29, 33

F1_median_edu$hs <- cut(F1_median_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                         c("Baseline", "Pre-symptomatic",
                           "Symptomatic", "Recovery"))

F2_median_edu$hs <- cut(F2_median_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))

M1_median_edu$hs <- cut(M1_median_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))

M2_median_edu$hs <- cut(M2_median_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))


# calculating median & sd per health state (hs)

F1_med_edu_hs <- setNames(aggregate(F1_median_edu$Med.edu,list(F1_median_edu$hs),median), c("hs","Med.edu.hs"))
F1_sd_edu_hs <- setNames(aggregate(F1_median_edu$Med.edu,list(F1_median_edu$hs),sd), c("hs","SD.edu.hs"))

F2_med_edu_hs <- setNames(aggregate(F2_median_edu$Med.edu,list(F2_median_edu$hs),median), c("hs","Med.edu.hs"))
F2_sd_edu_hs <- setNames(aggregate(F2_median_edu$Med.edu,list(F2_median_edu$hs),sd), c("hs","SD.edu.hs"))

M1_med_edu_hs <- setNames(aggregate(M1_median_edu$Med.edu,list(M1_median_edu$hs),median), c("hs","Med.edu.hs"))
M1_sd_edu_hs <- setNames(aggregate(M1_median_edu$Med.edu,list(M1_median_edu$hs),sd), c("hs","SD.edu.hs"))

M2_med_edu_hs <- setNames(aggregate(M2_median_edu$Med.edu,list(M2_median_edu$hs),median), c("hs","Med.edu.hs"))
M2_sd_edu_hs <- setNames(aggregate(M2_median_edu$Med.edu,list(M2_median_edu$hs),sd), c("hs","SD.edu.hs"))


# joining the median and sd for each cage together
# and recreating the Cage column for future needs

F1_med_sd_edu_hs <- inner_join(F1_med_edu_hs, F1_sd_edu_hs)
F1_med_sd_edu_hs$Cage <- "F1"

F2_med_sd_edu_hs <- inner_join(F2_med_edu_hs, F2_sd_edu_hs)
F2_med_sd_edu_hs$Cage <- "F2"

M1_med_sd_edu_hs <- inner_join(M1_med_edu_hs, M1_sd_edu_hs)
M1_med_sd_edu_hs$Cage <- "M1"

M2_med_sd_edu_hs <- inner_join(M2_med_edu_hs, M2_sd_edu_hs)
M2_med_sd_edu_hs$Cage <- "M2"

# joining info of median edu stats from all cages together
all.cages_edu_hs <- rbind(F1_med_sd_edu_hs, F2_med_sd_edu_hs, M1_med_sd_edu_hs, M2_med_sd_edu_hs)

# turning cage into a factor
all.cages_edu_hs$Cage <- factor(all.cages_edu_hs$Cage, 
                                levels = c("F1", "F2", "M1", "M2"))

## Medians of Double-computed Medians -----------------------------------

all_med_edu_hs <- setNames(aggregate(all.cages_edu_hs$Med.edu.hs,list(all.cages_edu_hs$hs),median), c("hs","Med2.edu.hs"))
all_sd_edu_hs <- setNames(aggregate(F1_median_edu$Med.edu,list(F1_median_edu$hs),sd), c("hs","SD2.edu.hs"))

all_med_sd_edu_hs <- inner_join(all_med_edu_hs, all_sd_edu_hs)

## Double-computed means -------------------------------

# getting the mean edu & sd per replicate per day w/in each cage
# grouping days by using list()
# changing names of columns using setNames(), as per
# https://statisticsglobe.com/set-column-names-within-aggregate-function-in-r

F1_mean_edu <- setNames(aggregate(F1_edu$prop_edu,list(F1_edu$day),mean), c("Day","Mean.edu"))
F1_mean_sd_edu <- setNames(aggregate(F1_edu$prop_edu,list(F1_edu$day),sd), c("Day","SD.edu"))

F2_mean_edu <- setNames(aggregate(F2_edu$prop_edu,list(F2_edu$day),mean), c("Day", "Mean.edu"))
F2_mean_sd_edu <- setNames(aggregate(F2_edu$prop_edu,list(F2_edu$day),sd), c("Day","SD.edu"))

M1_mean_edu <- setNames(aggregate(M1_edu$prop_edu,list(M1_edu$day),mean), c("Day","Mean.edu"))
M1_mean_sd_edu <- setNames(aggregate(M1_edu$prop_edu,list(M1_edu$day),sd), c("Day","SD.edu"))

M2_mean_edu <- setNames(aggregate(M2_edu$prop_edu,list(M2_edu$day),mean), c("Day","Mean.edu"))
M2_mean_sd_edu <- setNames(aggregate(M2_edu$prop_edu,list(M2_edu$day),sd), c("Day","SD.edu"))

# creating data frames with the values of mean & sd
# using inner_join()

F1_mean_std_edu <- inner_join(F1_mean_edu,F1_mean_sd_edu)
F1_mean_std_edu$Cage <- "F1"

F2_mean_std_edu <- inner_join(F2_mean_edu,F2_mean_sd_edu)
F2_mean_std_edu$Cage <- "F2"

M1_mean_std_edu <- inner_join(M1_mean_edu,M1_mean_sd_edu)
M1_mean_std_edu$Cage <- "M1"

M2_mean_std_edu <- inner_join(M2_mean_edu,M2_mean_sd_edu)
M2_mean_std_edu$Cage <- "M2"

# Combining all above means & SDs together into one data frame

all.cages_mean_edu <- rbind(F1_mean_std_edu,F2_mean_std_edu,M1_mean_std_edu,M2_mean_std_edu)

# now calculating means per health state per cage

F1_mean_edu$hs <- cut(F1_mean_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))

F2_mean_edu$hs <- cut(F2_mean_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))

M1_mean_edu$hs <- cut(M1_mean_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))

M2_mean_edu$hs <- cut(M2_mean_edu$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))


# calculating mean & sd per health state (hs)

F1_mean_edu_hs <- setNames(aggregate(F1_mean_edu$Mean.edu,list(F1_mean_edu$hs),mean), c("hs","Mean.edu.hs"))
F1_mean_sd_edu_hs <- setNames(aggregate(F1_mean_edu$Mean.edu,list(F1_mean_edu$hs),sd), c("hs","SD.edu.hs"))

F2_mean_edu_hs <- setNames(aggregate(F2_mean_edu$Mean.edu,list(F2_mean_edu$hs),mean), c("hs","Mean.edu.hs"))
F2_mean_sd_edu_hs <- setNames(aggregate(F2_mean_edu$Mean.edu,list(F2_mean_edu$hs),sd), c("hs","SD.edu.hs"))

M1_mean_edu_hs <- setNames(aggregate(M1_mean_edu$Mean.edu,list(M1_mean_edu$hs),mean), c("hs","Mean.edu.hs"))
M1_mean_sd_edu_hs <- setNames(aggregate(M1_mean_edu$Mean.edu,list(M1_mean_edu$hs),sd), c("hs","SD.edu.hs"))

M2_mean_edu_hs <- setNames(aggregate(M2_mean_edu$Mean.edu,list(M2_mean_edu$hs),mean), c("hs","Mean.edu.hs"))
M2_mean_sd_edu_hs <- setNames(aggregate(M2_mean_edu$Mean.edu,list(M2_mean_edu$hs),sd), c("hs","SD.edu.hs"))

# joining mean & sds together

F1_mean_sd_edu_hs <- inner_join(F1_mean_edu_hs, F1_mean_sd_edu_hs)
F1_mean_sd_edu_hs$Cage <- "F1"

F2_mean_sd_edu_hs <- inner_join(F2_mean_edu_hs, F2_mean_sd_edu_hs)
F2_mean_sd_edu_hs$Cage <- "F2"

M1_mean_sd_edu_hs <- inner_join(M1_mean_edu_hs, M1_mean_sd_edu_hs)
M1_mean_sd_edu_hs$Cage <- "M1"

M2_mean_sd_edu_hs <- inner_join(M2_mean_edu_hs, M2_mean_sd_edu_hs)
M2_mean_sd_edu_hs$Cage <- "M2"

# joining info of edu/fitc stats from all cages together
all.cages_mean_edu_hs <- rbind(F1_mean_sd_edu_hs, F2_mean_sd_edu_hs, M1_mean_sd_edu_hs, M2_mean_sd_edu_hs)

# turning cage into a factor

all.cages_mean_edu_hs$Cage <- factor(all.cages_mean_edu_hs$Cage, 
                                     levels = c("F1", "F2",
                                                "M1", "M2"))

# Stats ---------------------------------

# info on friedman & follow-up wilcoxon test:
# https://www.datanovia.com/en/lessons/friedman-test-in-r/

## Double-computed medians -------------------

### EdU -------------

# friedman test
all.cages_edu_hs %>% friedman_test(Med.edu.hs ~ hs | Cage)

# follow-up wilcoxon sign-rank test
all.cages_edu_hs %>% wilcox_test(Med.edu.hs ~ hs, 
                                 ref.group = "Baseline", 
                                 paired = TRUE, 
                                 p.adjust.method = "bonferroni")

## Means ----------------------------

### EdU ------------------------------

# friedman test
all.cages_mean_edu_hs %>% friedman_test(Mean.edu.hs ~ hs | Cage)

# follow-up wilcoxon sign-rank test
all.cages_mean_edu_hs %>% wilcox_test(Mean.edu.hs ~ hs, 
                                  ref.group = "Baseline", 
                                  paired = TRUE, 
                                  p.adjust.method = "bonferroni")