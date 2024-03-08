# Statistical tests on alpha diversity data
# For DSS #1

# Libraries --------------------------------------------------------------

library(tidyverse)
library(rstatix) # for pipe-friendly friedman test

# Data manipulation --------------------------------------------------------

## Medians -----------------------------------------------------------------

# compute median [alpha diversity metric] per health state per cage
# as well as median absolute deviation (mad)

### Observed features ("OTUs") --------------------

#### F1 -------------
F1_med_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$observed_features,
                                    list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                    median), c("health.state","Med.OTU"))

F1_mad_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$observed_features,
                                    list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                    mad), c("health.state","MAD.OTU"))

#### F2 -------------
F2_med_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$observed_features,
                                 list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                 median), c("health.state","Med.OTU"))

F2_mad_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$observed_features,
                                    list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                    mad), c("health.state","MAD.OTU"))

#### M1 -------------
M1_med_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$observed_features,
                                 list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                 median), c("health.state","Med.OTU"))

M1_mad_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$observed_features,
                                    list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                    mad), c("health.state","MAD.OTU"))

#### M2 -------------
M2_med_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$observed_features,
                                 list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                 median), c("health.state","Med.OTU"))

M2_mad_otu <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$observed_features,
                                    list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                    mad), c("health.state","MAD.OTU"))

### Shannon --------------------

all_alpha$shannon_entropy

#### F1 -------------
F1_med_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$shannon_entropy,
                                 list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                 median), c("health.state","Med.shannon"))

F1_mad_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$shannon_entropy,
                                    list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                    mad), c("health.state","MAD.shannon"))

#### F2 -------------
F2_med_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$shannon_entropy,
                                 list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                 median), c("health.state","Med.shannon"))

F2_mad_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$shannon_entropy,
                                    list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                    mad), c("health.state","MAD.shannon"))

#### M1 -------------
M1_med_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$shannon_entropy,
                                 list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                 median), c("health.state","Med.shannon"))

M1_mad_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$shannon_entropy,
                                    list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                    mad), c("health.state","MAD.shannon"))

#### M2 -------------
M2_med_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$shannon_entropy,
                                 list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                 median), c("health.state","Med.shannon"))

M2_mad_shannon <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$shannon_entropy,
                                    list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                    mad), c("health.state","MAD.shannon"))

### Simpson E --------------------

all_alpha$simpson_e

#### F1 -------------
F1_med_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$simpson_e,
                                 list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                 median), c("health.state","Med.simp_e"))

F1_mad_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$simpson_e,
                                    list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                    mad), c("health.state","MAD.simp_e"))

#### F2 -------------
F2_med_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$simpson_e,
                                 list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                 median), c("health.state","Med.simp_e"))

F2_mad_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$simpson_e,
                                    list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                    mad), c("health.state","MAD.simp_e"))

#### M1 -------------
M1_med_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$simpson_e,
                                 list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                 median), c("health.state","Med.simp_e"))

M1_mad_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$simpson_e,
                                    list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                    mad), c("health.state","MAD.simp_e"))

#### M2 -------------
M2_med_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$simpson_e,
                                 list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                 median), c("health.state","Med.simp_e"))

M2_mad_simp_e <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$simpson_e,
                                    list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                    mad), c("health.state","MAD.simp_e"))

### Simpson --------------------

#### F1 -------------
F1_med_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$simpson,
                                 list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                 median), c("health.state","Med.simp"))

F1_mad_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1"))$simpson,
                                    list(subset(all_alpha, subset = (cage == "F1"))$`health-state`),
                                    mad), c("health.state","MAD.simp"))

#### F2 -------------
F2_med_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$simpson,
                                 list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                 median), c("health.state","Med.simp"))

F2_mad_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2"))$simpson,
                                    list(subset(all_alpha, subset = (cage == "F2"))$`health-state`),
                                    mad), c("health.state","MAD.simp"))

#### M1 -------------
M1_med_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$simpson,
                                 list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                 median), c("health.state","Med.simp"))

M1_mad_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1"))$simpson,
                                    list(subset(all_alpha, subset = (cage == "M1"))$`health-state`),
                                    mad), c("health.state","MAD.simp"))

#### M2 -------------
M2_med_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$simpson,
                                 list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                 median), c("health.state","Med.simp"))

M2_mad_simp <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2"))$simpson,
                                    list(subset(all_alpha, subset = (cage == "M2"))$`health-state`),
                                    mad), c("health.state","MAD.simp"))

### Joining ---------------------------------

# joining the median and mad for each cage & alpha metric together
# and recreating the Cage column for future needs

#### OTU ----------------

F1_med_mad_otu <- inner_join(F1_med_otu, F1_mad_otu)
F1_med_mad_otu$cage <- "F1"

F2_med_mad_otu <- inner_join(F2_med_otu, F2_mad_otu)
F2_med_mad_otu$cage <- "F2"

M1_med_mad_otu <- inner_join(M1_med_otu, M1_mad_otu)
M1_med_mad_otu$cage <- "M1"

M2_med_mad_otu <- inner_join(M2_med_otu, M2_mad_otu)
M2_med_mad_otu$cage <- "M2"

#### Shannon ----------------

F1_med_mad_shannon <- inner_join(F1_med_shannon, F1_mad_shannon)
F1_med_mad_shannon$cage <- "F1"

F2_med_mad_shannon <- inner_join(F2_med_shannon, F2_mad_shannon)
F2_med_mad_shannon$cage <- "F2"

M1_med_mad_shannon <- inner_join(M1_med_shannon, M1_mad_shannon)
M1_med_mad_shannon$cage <- "M1"

M2_med_mad_shannon <- inner_join(M2_med_shannon, M2_mad_shannon)
M2_med_mad_shannon$cage <- "M2"

#### Simpson E ----------------

F1_med_mad_simp_e <- inner_join(F1_med_simp_e, F1_mad_simp_e)
F1_med_mad_simp_e$cage <- "F1"

F2_med_mad_simp_e <- inner_join(F2_med_simp_e, F2_mad_simp_e)
F2_med_mad_simp_e$cage <- "F2"

M1_med_mad_simp_e <- inner_join(M1_med_simp_e, M1_mad_simp_e)
M1_med_mad_simp_e$cage <- "M1"

M2_med_mad_simp_e <- inner_join(M2_med_simp_e, M2_mad_simp_e)
M2_med_mad_simp_e$cage <- "M2"

#### Simpson ----------------

F1_med_mad_simp <- inner_join(F1_med_simp, F1_mad_simp)
F1_med_mad_simp$cage <- "F1"

F2_med_mad_simp <- inner_join(F2_med_simp, F2_mad_simp)
F2_med_mad_simp$cage <- "F2"

M1_med_mad_simp <- inner_join(M1_med_simp, M1_mad_simp)
M1_med_mad_simp$cage <- "M1"

M2_med_mad_simp <- inner_join(M2_med_simp, M2_mad_simp)
M2_med_mad_simp$cage <- "M2"

#### Joining from all cages together (per alpha metric) --------------------------
all.cages_otu <- rbind(F1_med_mad_otu, F2_med_mad_otu, M1_med_mad_otu, M2_med_mad_otu)
all.cages_shannon <- rbind(F1_med_mad_shannon, F2_med_mad_shannon, M1_med_mad_shannon, M2_med_mad_shannon)
all.cages_simp_e <- rbind(F1_med_mad_simp_e, F2_med_mad_simp_e, M1_med_mad_simp_e, M2_med_mad_simp_e)
all.cages_simp <- rbind(F1_med_mad_simp, F2_med_mad_simp, M1_med_mad_simp, M2_med_mad_simp)

# turning cage into a factor
all.cages_otu$cage <- factor(all.cages_otu$cage, 
                                levels = c("F1", "F2", "M1", "M2"))

all.cages_shannon$cage <- factor(all.cages_shannon$cage, 
                             levels = c("F1", "F2", "M1", "M2"))

all.cages_simp_e$cage <- factor(all.cages_simp_e$cage, 
                             levels = c("F1", "F2", "M1", "M2"))

all.cages_simp$cage <- factor(all.cages_simp$cage, 
                             levels = c("F1", "F2", "M1", "M2"))

# Stats ---------------------------------

# info on friedman & follow-up wilcoxon test:
# https://www.datanovia.com/en/lessons/friedman-test-in-r/

## Friedman test ---------------------------------------------
all.cages_otu %>% friedman_test(Med.OTU ~ health.state | cage)
# p = 0.272

all.cages_shannon %>% friedman_test(Med.shannon ~ health.state | cage)
# p = 0.0440

all.cages_simp_e %>% friedman_test(Med.simp_e ~ health.state | cage)
# p = 0.127

all.cages_simp %>% friedman_test(Med.simp ~ health.state | cage)
# p = 0.145

## Wilcoxon sign-rank test ----------------------------------------
all.cages_otu %>% wilcox_test(Med.OTU ~ health.state, 
                                 ref.group = "BA", 
                                 paired = TRUE, 
                                 p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_shannon %>% wilcox_test(Med.shannon ~ health.state, 
                              ref.group = "BA", 
                              paired = TRUE, 
                              p.adjust.method = "bonferroni")
# p.adj = 0.375 for all

all.cages_simp_e %>% wilcox_test(Med.simp_e ~ health.state, 
                              ref.group = "BA", 
                              paired = TRUE, 
                              p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_simp %>% wilcox_test(Med.simp ~ health.state, 
                              ref.group = "BA", 
                              paired = TRUE, 
                              p.adjust.method = "bonferroni")
# p.adj > 0.05 for all
