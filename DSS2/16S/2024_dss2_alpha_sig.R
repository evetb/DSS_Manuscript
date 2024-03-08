# stats on alpha significance
# per cage & sorted fraction (sf)
# between health states (hs)

# Libraries ----------------------------------------------------

library(tidyverse)
library(rstatix) # for pipe-friendly Friedman test

# Data Manipulation --------------------------------------------

## Medians -----------------------------------------------------------------

# compute median [alpha diversity metric] 
# per health state per cage per sf

### OTUs --------------------

subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))

#### F1 -------------
F1_med_otu_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$observed_features,
                                 list(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$hs),
                                 median), c("health.state","Med.OTU"))

F1_med_otu_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$observed_features,
                                      list(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$hs),
                                      median), c("health.state","Med.OTU"))


#### F2 -------------

F2_med_otu_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$observed_features,
                                      list(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$hs),
                                      median), c("health.state","Med.OTU"))

F2_med_otu_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$observed_features,
                                       list(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$hs),
                                       median), c("health.state","Med.OTU"))

#### M1 -------------

M1_med_otu_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$observed_features,
                                      list(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$hs),
                                      median), c("health.state","Med.OTU"))

M1_med_otu_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$observed_features,
                                       list(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$hs),
                                       median), c("health.state","Med.OTU"))

#### M2 -------------

M2_med_otu_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$observed_features,
                                      list(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$hs),
                                      median), c("health.state","Med.OTU"))

M2_med_otu_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$observed_features,
                                       list(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$hs),
                                       median), c("health.state","Med.OTU"))


### Shannon --------------------

#### F1 -------------
F1_med_shannon_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$shannon_entropy,
                                     list(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$hs),
                                     median), c("health.state","Med.shannon"))

F1_med_shannon_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$shannon_entropy,
                                          list(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$hs),
                                          median), c("health.state","Med.shannon"))

#### F2 -------------

F2_med_shannon_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$shannon_entropy,
                                          list(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$hs),
                                          median), c("health.state","Med.shannon"))

F2_med_shannon_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$shannon_entropy,
                                           list(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$hs),
                                           median), c("health.state","Med.shannon"))

#### M1 -------------

M1_med_shannon_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$shannon_entropy,
                                          list(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$hs),
                                          median), c("health.state","Med.shannon"))

M1_med_shannon_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$shannon_entropy,
                                           list(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$hs),
                                           median), c("health.state","Med.shannon"))

#### M2 -------------

M2_med_shannon_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$shannon_entropy,
                                          list(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$hs),
                                          median), c("health.state","Med.shannon"))

M2_med_shannon_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$shannon_entropy,
                                           list(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$hs),
                                           median), c("health.state","Med.shannon"))

### Simpson E --------------------

#### F1 -------------
F1_med_simp_e_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$simpson_e,
                                    list(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$hs),
                                    median), c("health.state","Med.simp_e"))


F1_med_simp_e_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$simpson_e,
                                         list(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$hs),
                                         median), c("health.state","Med.simp_e"))

#### F2 -------------

F2_med_simp_e_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$simpson_e,
                                         list(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$hs),
                                         median), c("health.state","Med.simp_e"))


F2_med_simp_e_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$simpson_e,
                                          list(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$hs),
                                          median), c("health.state","Med.simp_e"))


#### M1 -------------

M1_med_simp_e_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$simpson_e,
                                         list(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$hs),
                                         median), c("health.state","Med.simp_e"))


M1_med_simp_e_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$simpson_e,
                                          list(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$hs),
                                          median), c("health.state","Med.simp_e"))


#### M2 -------------
M2_med_simp_e_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$simpson_e,
                                         list(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$hs),
                                         median), c("health.state","Med.simp_e"))


M2_med_simp_e_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$simpson_e,
                                          list(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$hs),
                                          median), c("health.state","Med.simp_e"))

### Simpson --------------------

#### F1 -------------
F1_med_simp_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$simpson,
                                  list(subset(all_alpha, subset = (cage == "F1" & sf == "EdU+"))$hs),
                                  median), c("health.state","Med.simp"))

F1_med_simp_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$simpson,
                                       list(subset(all_alpha, subset = (cage == "F1" & sf == "Whole"))$hs),
                                       median), c("health.state","Med.simp"))

#### F2 -------------

F2_med_simp_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$simpson,
                                       list(subset(all_alpha, subset = (cage == "F2" & sf == "EdU+"))$hs),
                                       median), c("health.state","Med.simp"))

F2_med_simp_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$simpson,
                                        list(subset(all_alpha, subset = (cage == "F2" & sf == "Whole"))$hs),
                                        median), c("health.state","Med.simp"))

#### M1 -------------

M1_med_simp_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$simpson,
                                       list(subset(all_alpha, subset = (cage == "M1" & sf == "EdU+"))$hs),
                                       median), c("health.state","Med.simp"))

M1_med_simp_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$simpson,
                                        list(subset(all_alpha, subset = (cage == "M1" & sf == "Whole"))$hs),
                                        median), c("health.state","Med.simp"))

#### M2 ------------------

M2_med_simp_epos <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$simpson,
                                       list(subset(all_alpha, subset = (cage == "M2" & sf == "EdU+"))$hs),
                                       median), c("health.state","Med.simp"))

M2_med_simp_whole <- setNames(aggregate(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$simpson,
                                        list(subset(all_alpha, subset = (cage == "M2" & sf == "Whole"))$hs),
                                        median), c("health.state","Med.simp"))

### Joining ---------------------------------

# recreating the Cage column for future needs

#### OTU ----------------

F1_med_otu_epos$cage <- "F1"
F2_med_otu_epos$cage <- "F2"
M1_med_otu_epos$cage <- "M1"
M2_med_otu_epos$cage <- "M2"

F1_med_otu_whole$cage <- "F1"
F2_med_otu_whole$cage <- "F2"
M1_med_otu_whole$cage <- "M1"
M2_med_otu_whole$cage <- "M2"

#### Shannon ----------------

F1_med_shannon_epos$cage <- "F1"
F2_med_shannon_epos$cage <- "F2"
M1_med_shannon_epos$cage <- "M1"
M2_med_shannon_epos$cage <- "M2"

F1_med_shannon_whole$cage <- "F1"
F2_med_shannon_whole$cage <- "F2"
M1_med_shannon_whole$cage <- "M1"
M2_med_shannon_whole$cage <- "M2"

#### Simpson E ----------------

F1_med_simp_e_epos$cage <- "F1"
F2_med_simp_e_epos$cage <- "F2"
M1_med_simp_e_epos$cage <- "M1"
M2_med_simp_e_epos$cage <- "M2"

F1_med_simp_e_whole$cage <- "F1"
F2_med_simp_e_whole$cage <- "F2"
M1_med_simp_e_whole$cage <- "M1"
M2_med_simp_e_whole$cage <- "M2"

#### Simpson ----------------

F1_med_simp_epos$cage <- "F1"
F2_med_simp_epos$cage <- "F2"
M1_med_simp_epos$cage <- "M1"
M2_med_simp_epos$cage <- "M2"

F1_med_simp_whole$cage <- "F1"
F2_med_simp_whole$cage <- "F2"
M1_med_simp_whole$cage <- "M1"
M2_med_simp_whole$cage <- "M2"

#### Joining from all cages together (per alpha metric) --------------------------
all.cages_otu_epos <- rbind(F1_med_otu_epos, F2_med_otu_epos, M1_med_otu_epos, M2_med_otu_epos)
all.cages_shannon_epos <- rbind(F1_med_shannon_epos, F2_med_shannon_epos, M1_med_shannon_epos, M2_med_shannon_epos)
all.cages_simp_e_epos <- rbind(F1_med_simp_e_epos, F2_med_simp_e_epos, M1_med_simp_e_epos, M2_med_simp_e_epos)
all.cages_simp_epos <- rbind(F1_med_simp_epos, F2_med_simp_epos, M1_med_simp_epos, M2_med_simp_epos)

all.cages_otu_whole <- rbind(F1_med_otu_whole, F2_med_otu_whole, M1_med_otu_whole, M2_med_otu_whole)
all.cages_shannon_whole <- rbind(F1_med_shannon_whole, F2_med_shannon_whole, M1_med_shannon_whole, M2_med_shannon_whole)
all.cages_simp_e_whole <- rbind(F1_med_simp_e_whole, F2_med_simp_e_whole, M1_med_simp_e_whole, M2_med_simp_e_whole)
all.cages_simp_whole <- rbind(F1_med_simp_whole, F2_med_simp_whole, M1_med_simp_whole, M2_med_simp_whole)

# turning cage into a factor
all.cages_otu_epos$cage <- factor(all.cages_otu_epos$cage, 
                             levels = c("F1", "F2", "M1", "M2"))
all.cages_shannon_epos$cage <- factor(all.cages_shannon_epos$cage, 
                                 levels = c("F1", "F2", "M1", "M2"))
all.cages_simp_e_epos$cage <- factor(all.cages_simp_e_epos$cage, 
                                levels = c("F1", "F2", "M1", "M2"))
all.cages_simp_epos$cage <- factor(all.cages_simp_epos$cage, 
                              levels = c("F1", "F2", "M1", "M2"))

all.cages_otu_whole$cage <- factor(all.cages_otu_whole$cage, 
                                  levels = c("F1", "F2", "M1", "M2"))
all.cages_shannon_whole$cage <- factor(all.cages_shannon_whole$cage, 
                                      levels = c("F1", "F2", "M1", "M2"))
all.cages_simp_e_whole$cage <- factor(all.cages_simp_e_whole$cage, 
                                     levels = c("F1", "F2", "M1", "M2"))
all.cages_simp_whole$cage <- factor(all.cages_simp_whole$cage, 
                                   levels = c("F1", "F2", "M1", "M2"))

# Stats ---------------------------------

# info on friedman & follow-up wilcoxon test:
# https://www.datanovia.com/en/lessons/friedman-test-in-r/

## Friedman test ---------------------------------------------
all.cages_otu_epos %>% friedman_test(Med.OTU ~ health.state | cage)
# p = 0.127

all.cages_otu_whole %>% friedman_test(Med.OTU ~ health.state | cage)
# p = 0.308

all.cages_shannon_epos %>% friedman_test(Med.shannon ~ health.state | cage)
# p = 0.0576 

all.cages_shannon_whole %>% friedman_test(Med.shannon ~ health.state | cage)
# p = 0.187

all.cages_simp_e_epos %>% friedman_test(Med.simp_e ~ health.state | cage)
# p = 0.112

all.cages_simp_e_whole %>% friedman_test(Med.simp_e ~ health.state | cage)
# p = 0.165

all.cages_simp_epos %>% friedman_test(Med.simp ~ health.state | cage)
# p = 0.0576

all.cages_simp_whole %>% friedman_test(Med.simp ~ health.state | cage)
# p = 0.127

## Wilcoxon sign-rank test ----------------------------------------
all.cages_otu_epos %>% wilcox_test(Med.OTU ~ health.state, 
                              ref.group = "BA", 
                              paired = TRUE, 
                              p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_shannon_epos %>% wilcox_test(Med.shannon ~ health.state, 
                                  ref.group = "BA", 
                                  paired = TRUE, 
                                  p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_simp_e_epos %>% wilcox_test(Med.simp_e ~ health.state, 
                                 ref.group = "BA", 
                                 paired = TRUE, 
                                 p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_simp_epos %>% wilcox_test(Med.simp ~ health.state, 
                               ref.group = "BA", 
                               paired = TRUE, 
                               p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_otu_whole %>% wilcox_test(Med.OTU ~ health.state, 
                                   ref.group = "BA", 
                                   paired = TRUE, 
                                   p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_shannon_whole %>% wilcox_test(Med.shannon ~ health.state, 
                                       ref.group = "BA", 
                                       paired = TRUE, 
                                       p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_simp_e_whole %>% wilcox_test(Med.simp_e ~ health.state, 
                                      ref.group = "BA", 
                                      paired = TRUE, 
                                      p.adjust.method = "bonferroni")
# p.adj > 0.05 for all

all.cages_simp_whole %>% wilcox_test(Med.simp ~ health.state, 
                                    ref.group = "BA", 
                                    paired = TRUE, 
                                    p.adjust.method = "bonferroni")
# p.adj > 0.05 for all