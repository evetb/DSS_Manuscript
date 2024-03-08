# Graphs of DSS #2 mouse health parameters 
# (weight, LCN2)

# Libraries ---------------------------------------------------------------

library(tidyverse)

# LCN2 --------------------------------------------------------------------

lcn2 <- read_csv("DSS2_LCN2.csv")

## LCN2 data manipulation -------------------------------------------------------------------

# make day into a factor & specifying the order
lcn2$day <- factor(lcn2$day, 
                   levels = c("8", "9", "10", "11", "12", "13", "14", "15", "16",
                              "17", "18", "26", "29", "33"))  

# turning pg/mL to ng/g
# there are 1000 pg in one ng
# i diluted 100 mg of feces in 1 mL
# (pg/mL) * (ng/1000 pg) * (mL/0.1g) = ng/100g 
# divide everything by 100 to get ng LCN-2 per g feces

lcn2$lcn2_ngg <- (lcn2$lcn2)/100

## LCN2 plots --------------------------------------------------------------------------------

### Log10 scale on ng/g ---------------------------------------------

lcn2 %>% ggplot(aes(x=day, y=lcn2_ngg, color = cage, group = cage)) +
  geom_point(size = 4) + geom_line(linewidth = 2) +
  labs(#title = "DSS #2 LCN-2", 
       x = "Day", y = "LCN-2 (ng/g)",
       color = "Cage") +
  theme_minimal() + theme_bw() +
  scale_y_log10() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

#### coloured by health states (hs) ---------------------------------------------------

lcn2$day <- as.numeric(levels(lcn2$day))[lcn2$day]

lcn2$hs <- cut(lcn2$day, c(-Inf, 7, 11, 15, Inf), 
               c("Baseline", "Pre-symptomatic",
                 "Symptomatic", "Recovery"))

lcn2$day <- factor(lcn2$day, 
                   levels=c("8", "9", "10", "11", "12", "13", "14", 
                            "15", "16", "17", "18",
                            "26", "29", "33"))

lcn2$sex <- with(lcn2, ifelse(cage == "F1", "Female",
                              ifelse(cage == "F2", "Female",
                                     ifelse(cage == "M1", "Male", "Male"))))

lcn2 %>% ggplot(aes(x=as.numeric(levels(day)[day]), y=lcn2_ngg, group = cage)) +
  geom_point((aes(color = hs, shape = sex)), size = 4) + 
  geom_line((aes(color = hs)), linewidth = 2) +
  scale_color_manual(values=c("#CCBB44", "#EE6677", "#228833")) +
  labs(title = "DSS #2 LCN-2", x = "Day", y = "LCN-2 (ng/g)",
       color = "Health State", shape = "Sex") +
  theme_minimal() + theme_bw() +
  scale_y_log10() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))

# Weights ------------------------------------------------------------------

weights <- read_csv("DSS2_weights.csv")

weights$day <- factor(weights$day, 
                      levels = c("1", "5", "7", "8", 
                                 "9", "10", "11", "12", 
                                 "13", "14", "15", "16", 
                                 "17", "18", "26", "29", "33"))

## Median & SD weights ----------------------------------------------------------

# first making data frame for each cage, 
# since I want to get the median weight of each mouse 
# in each cage per day
# https://dplyr.tidyverse.org/reference/filter.html

F1_tbl <- filter(weights, cage == "F1")
F2_tbl <- filter(weights, cage == "F2")
M1_tbl <- filter(weights, cage == "M1")
M2_tbl <- filter(weights, cage == "M2")

# getting the median weight & sd per mouse per day, grouping days by using list()
# changing names of columns using setNames(), as per
# https://statisticsglobe.com/set-column-names-within-aggregate-function-in-r

F1_median_weights <- setNames(aggregate(F1_tbl$weight,list(F1_tbl$day),median), c("Day","Med.Wt"))
F1_sd_weights <- setNames(aggregate(F1_tbl$weight,list(F1_tbl$day),sd), c("Day","SD.Wt"))

F2_median_weights <- setNames(aggregate(F2_tbl$weight,list(F2_tbl$day),median), c("Day", "Med.Wt"))
F2_sd_weights <- setNames(aggregate(F2_tbl$weight,list(F2_tbl$day),sd), c("Day","SD.Wt"))

M1_median_weights <- setNames(aggregate(M1_tbl$weight,list(M1_tbl$day),median), c("Day","Med.Wt"))
M1_sd_weights <- setNames(aggregate(M1_tbl$weight,list(M1_tbl$day),sd), c("Day","SD.Wt"))

M2_median_weights <- setNames(aggregate(M2_tbl$weight,list(M2_tbl$day),median), c("Day","Med.Wt"))
M2_sd_weights <- setNames(aggregate(M2_tbl$weight,list(M2_tbl$day),sd), c("Day","SD.Wt"))

# creating data frames with the values of median & sd
# using inner_join() 

F1_med_std_weights <- inner_join(F1_median_weights,F1_sd_weights)
F1_med_std_weights$Cage <- "F1"

F2_med_std_weights <- inner_join(F2_median_weights,F2_sd_weights)
F2_med_std_weights$Cage <- "F2"

M1_med_std_weights <- inner_join(M1_median_weights,M1_sd_weights)
M1_med_std_weights$Cage <- "M1"

M2_med_std_weights <- inner_join(M2_median_weights,M2_sd_weights)
M2_med_std_weights$Cage <- "M2"

# Combining all above medians & SDs together into one data frame

all.cages <- rbind(F1_med_std_weights,F2_med_std_weights,M1_med_std_weights,M2_med_std_weights)

## Relative change in weight -------------------------------------------------------------

# Getting the relative change in weight per day, per mouse
weights.rel <- weights %>% group_by(day) %>% group_by(mouse) %>%
  mutate(rel.change.wt = weight - first(weight))

# First extracting individual cages from the data frame
F1_weights.rel <- weights.rel %>% filter(cage == "F1")
F2_weights.rel <- weights.rel %>% filter(cage == "F2")
M1_weights.rel <- weights.rel %>% filter(cage == "M1")
M2_weights.rel <- weights.rel %>% filter(cage == "M2")

# Then calculating median & SD 
# of the relative change in weight 
# Per day, per mouse

F1_weights.rel.med <- setNames(aggregate(F1_weights.rel$rel.change.wt,list(F1_weights.rel$day), median),
                               c("day","rel.med.wt"))
F1_weights.rel.sd <- setNames(aggregate(F1_weights.rel$rel.change.wt,list(F1_weights.rel$day),sd), 
                              c("day","rel.sd.wt"))

F2_weights.rel.med <- setNames(aggregate(F2_weights.rel$rel.change.wt,list(F2_weights.rel$day), median),
                               c("day","rel.med.wt"))
F2_weights.rel.sd <- setNames(aggregate(F2_weights.rel$rel.change.wt,list(F2_weights.rel$day),sd), 
                              c("day","rel.sd.wt"))

M1_weights.rel.med <- setNames(aggregate(M1_weights.rel$rel.change.wt,list(M1_weights.rel$day), median),
                               c("day","rel.med.wt"))
M1_weights.rel.sd <- setNames(aggregate(M1_weights.rel$rel.change.wt,list(M1_weights.rel$day),sd), 
                              c("day","rel.sd.wt"))

M2_weights.rel.med <- setNames(aggregate(M2_weights.rel$rel.change.wt,list(M2_weights.rel$day), median),
                               c("day","rel.med.wt"))
M2_weights.rel.sd <- setNames(aggregate(M2_weights.rel$rel.change.wt,list(M2_weights.rel$day),sd), 
                              c("day","rel.sd.wt"))

# Then joining sd data 
# to the _weights.rel data frames

F1_weights.rel.med.sd <- inner_join(F1_weights.rel.med,F1_weights.rel.sd)
F1_weights.rel.med.sd$cage <- "F1"

F2_weights.rel.med.sd <- inner_join(F2_weights.rel.med,F2_weights.rel.sd)
F2_weights.rel.med.sd$cage <- "F2"

M1_weights.rel.med.sd <- inner_join(M1_weights.rel.med,M1_weights.rel.sd)
M1_weights.rel.med.sd$cage <- "M1"

M2_weights.rel.med.sd <- inner_join(M2_weights.rel.med,M2_weights.rel.sd)
M2_weights.rel.med.sd$cage <- "M2"

# Then putting all cages 
# into one data frame for plotting

all.cages.rel <- rbind(F1_weights.rel.med.sd,F2_weights.rel.med.sd,
                       M1_weights.rel.med.sd,M2_weights.rel.med.sd)

## Weight plots -----------------------------------------------------------------------------

### Relative Weights ----------------------------------------------------------------------

all.cages.rel.plot <- all.cages.rel %>% 
  ggplot(aes(x=day, y=rel.med.wt, group=cage, color = cage)) +
  geom_point(size = 4) + geom_line(linewidth = 2) +
  labs(#title = "DSS #2 Weights", 
       x = "Day", y = "Rel. Change in Weight",
       color = "Cage") +
  geom_errorbar(aes(ymin=rel.med.wt-rel.sd.wt, ymax=rel.med.wt+rel.sd.wt), linewidth = 1) +
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all.cages.rel.plot