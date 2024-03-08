# making line plots of gRodon output 
# for whole community per day

# Loading libraries -----------------------------------------------------------------------------------

library(tidyverse) 
library(rstatix)
library(lme4) # for lmer
library(lmerTest) # for stats from lmer
library(olsrr) # for normality checks

# Setting working directory
setwd("")

# Plotting n_le = 1000 with new health states (hs) -----------------------------

## Importing data ----------------------------------

F2_1k_new_hs <- read_csv("grodon_vals_1k_F2_new_hs.csv")
M2_1k_new_hs <- read_csv("grodon_vals_1k_M2_new_hs.csv")

# making hs into a factor

F2_1k_new_hs$hs <- factor(F2_1k_new_hs$hs, 
                          levels = c("Baseline", "Pre-symptomatic", 
                                     "Symptomatic", "Recovery"))

M2_1k_new_hs$hs <- factor(M2_1k_new_hs$hs, 
                          levels = c("Baseline", "Pre-symptomatic", 
                                     "Symptomatic", "Recovery"))

## Plotting -----------------------------------------

### Line plots -------------------------------------

F2_1k_new_hs %>% ggplot(aes(x = day, y = min_DT, group = cage, colour = hs)) +
  geom_point(size = 4) + geom_line(linewidth = 2) +
  labs(title = "Cage F2", x = "Day", 
       y = "Min. Doubling Time (Hrs)", colour = "Health State") +
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), linewidth = 1) +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  #facet_grid(. ~ hs, scales = "free", space = "free") + 
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))


M2_1k_new_hs %>% ggplot(aes(x = day, y = min_DT, group = cage, colour = hs)) +
  geom_point(size = 4) + geom_line(linewidth = 2) +
  labs(title = "Cage M2", x = "Day", 
       y = "Min. Doubling Time (Hrs)", colour = "Health State") +
  geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI), linewidth = 1) +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  #facet_grid(. ~ hs, scales = "free", space = "free") + 
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))

### Regression plots -----------------------------

F2_1k_new_hs %>% ggplot(aes(x = day, y = min_DT, group = cage, colour = hs)) +
  geom_smooth(method=lm, color= "black", se=TRUE) +
  geom_point(size = 4) +
  labs(title = "Cage F2", x = "Day", 
       y = "Min. Doubling Time (Hrs)", colour = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))


M2_1k_new_hs %>% ggplot(aes(x = day, y = min_DT, group = cage, colour = hs)) +
  geom_smooth(method=lm, color= "black", se=TRUE) +
  geom_point(size = 4) +
  labs(title = "Cage M2", x = "Day", 
       y = "Min. Doubling Time (Hrs)", colour = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))


#### Combining cages ------------------------------------

both_1k_new_hs <- full_join(F2_1k_new_hs, M2_1k_new_hs)


both_1k_new_hs$hs <- factor(both_1k_new_hs$hs, 
                            levels = c("Baseline", "Pre-symptomatic", 
                                       "Symptomatic", "Recovery"))
## making a column called sample
## which will be numeric
## one value per cage

both_1k_new_hs$sample <-  with(both_1k_new_hs,
                               ifelse(cage == "F1", 1,
                                      ifelse(cage == "F2", 2,
                                             ifelse(cage == "M1", 3, 4))))

both.lm <- lm(formula = min_DT~day, data = both_1k_new_hs)
summary(both.lm)

#### checking for normality ------------------------------------

# plot several residual graphs 
# (including residuals vs fitted, QQ plot)

plot(both.lm)

# plot histogram for normal distribution
# lme output doesn't play well w/olssr
ols_plot_resid_hist(both.lm)

# box plot
ols_plot_resid_box(both.lm)

# running tests for normality

ols_test_normality(both.lm)

# seems normal!

#### lmer on both cages ----------------------------------

# https://janajarecki.com/blog/repeated-measures-regression-in-r/
# https://diposit.ub.edu/dspace/bitstream/2445/188685/1/mixedmodels.pdf

# both.lmer <- lmer(formula = min_DT~day + (1|cage), 
#                 data = both_1k_new_hs)

both.lmer <- lmer(formula = min_DT~hs + (1|cage) + (1|day),
                  data = both_1k_new_hs)

summary(both.lmer)

# correlation using day 
# (since hs is not numeric)

correlation <- cor(both_1k_new_hs$day,
                   both_1k_new_hs$min_DT,
                   method = 'pearson',
                   use = "complete.obs")

# r = -0.6095734
# r^2 = 0.3715797

##### Plotting ---------------------------------------------------------------

both_1k_new_hs %>% ggplot(aes(x = day, y = min_DT, colour = hs)) +
  geom_smooth(method=lm, color= "black", se=TRUE) +
  geom_point(aes(shape = cage), size = 4) +
  labs(title = "Both Cages", x = "Day", 
       y = "Min. Doubling Time (Hrs)", colour = "Health State",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=24)
  )

# gRodon vs LCN2 -----------------------------------------

## Data import/export ---------------------------------------

# export both_1k_new_hs as made above
# to add in LCN2 quantities per day

write_csv(both_1k_new_hs, "both_1k_new_hs.csv")

grodon_lcn2 <- read_csv("grodon_lcn2.csv")

# specifying factor order

grodon_lcn2$hs <- factor(grodon_lcn2$hs, 
                         levels = c("Baseline", "Pre-symptomatic",
                                    "Symptomatic", "Recovery"))

## Plots -------------------------------------------

### Scatterplot ---------------------------------------

# standardized from log2 values

grodon_lcn2 %>% ggplot(aes(x = scale(log2(lcn2_ngg)), y = scale(log2(min_DT)))) +
  geom_point((aes(colour = hs)), size = 4, shape = 16) +
  scale_colour_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  geom_smooth(method=lm, color= "black", se=TRUE) +
  labs(x = "Std.log2(LCN2) (ng/g)", 
       y = "Std.log2(Min.DT) (Hrs)", 
       colour = "Health State") +
  #scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 30, hjust = 0.5), axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.text=element_text(size=16),
        legend.title=element_text(size=24))

# correlation
corr_grodon_lcn2 <- cor(grodon_lcn2$lcn2_ngg, 
                        grodon_lcn2$min_DT, 
                                  method = 'pearson',
                                  use = "complete.obs")
# r = -0.09411664
# r^2 = 0.008857943

## lmer --------------------------------------------

# on standardized log2 values

lmer_grodon_lcn2 <- grodon_lcn2 %>% 
  lmer(formula = scale(log2(min_DT))~scale(log2(lcn2_ngg)) + 
         (1|cage) + (1|day))

summary(lmer_grodon_lcn2)