# plotting alpha diversity for DSS #2

# Libraries --------------------------------------------------------------------
library(tidyverse)

# Data manipulation ------------------------------------------------------------

Epos_OTUs <- read_tsv("nnn-Epos-observed-otus.tsv")
Epos_shannon <- read_tsv("nnn-Epos-shannon-vector.tsv")
Epos_simp <- read_tsv("nnn-Epos-simpson-index.tsv")
Epos_simpE <- read_tsv("nnn-Epos-simpsonE-index.tsv")

W_OTUs <- read_tsv("nnn-W-observed-otus.tsv") 
W_shannon <- read_tsv("nnn-W-shannon-vector.tsv")
W_simp <- read_tsv("nnn-W-simpson-index.tsv")
W_simpE <- read_tsv("nnn-W-simpsonE-index.tsv")

# compile together, as per: 
# https://sparkbyexamples.com/r-programming/r-join-multiple-data-frames/

# need to put in data regarding day, health state (hs), sorted fraction (sf)

list_Epos_alpha <- list(Epos_OTUs, Epos_shannon, Epos_simp, Epos_simpE)
Epos_all_alpha <- list_Epos_alpha %>% reduce(inner_join, by = "...1")

list_W_alpha <- list(W_OTUs, W_shannon, W_simp, W_simpE)
W_all_alpha <- list_W_alpha %>% reduce(inner_join, by = "...1")

# compile together Epos and W - from:
# https://dplyr.tidyverse.org/reference/bind.html

all_alpha <- bind_rows(Epos_all_alpha, W_all_alpha)

# fixing first column

all_alpha <- all_alpha %>% rename("sample" = "...1")

# adding sorted fraction - from:
# https://stackoverflow.com/questions/39903376/if-column-contains-string-then-enter-value-for-that-row

all_alpha$sf <- ifelse(grepl("Epos", all_alpha$sample), "EdU+", "Whole")

# adding cage
# https://stackoverflow.com/questions/39903376/if-column-contains-string-then-enter-value-for-that-row
# https://www.marsja.se/r-add-column-to-dataframe-based-on-other-columns-conditions-dplyr/

all_alpha <- all_alpha %>% 
  mutate(cage = case_when(str_detect(sample, "F1") ~ "F1",
                          str_detect(sample, "F2") ~ "F2",
                          str_detect(sample, "M1") ~ "M1", TRUE ~ "M2"))

# adding sex

all_alpha$sex <- ifelse(grepl("F", all_alpha$sample), "Female", "Male")

# adding day

all_alpha <- all_alpha %>% 
  mutate(day = case_when(str_detect(sample, "Epos10") ~ 18,
                         str_detect(sample, "W10") ~ 18,
                         str_detect(sample, "Epos11") ~ 26,
                         str_detect(sample, "W11") ~ 26,
                         str_detect(sample, "Epos12") ~ 29,
                         str_detect(sample, "W12") ~ 29,
                         str_detect(sample, "Epos13") ~ 33,
                         str_detect(sample, "W13") ~ 33,
                         str_detect(sample, "Epos2") ~ 5,
                         str_detect(sample, "W2") ~ 5,
                         str_detect(sample, "3") ~ 7,
                         str_detect(sample, "4") ~ 8,
                         str_detect(sample, "5") ~ 11,
                         str_detect(sample, "6") ~ 12,
                         str_detect(sample, "7") ~ 13,
                         str_detect(sample, "8") ~ 14,
                         str_detect(sample, "9") ~ 15,
                         TRUE ~ 1))

# adding health state

all_alpha$hs <- ifelse(all_alpha$day >= 1 & all_alpha$day <=7, "Baseline",
                           ifelse(all_alpha$day >= 8 & all_alpha$day <= 11, "Pre-symptomatic",
                                  ifelse(all_alpha$day >= 12 & all_alpha$day <= 15, "Symptomatic", "Recovery")))


# Plots -----------------------------------------------------------------------

# Line plots ------------------------------------------------------------------

## F1 -------------------------------------------------------------------------

# subsetting notation from:
# https://stackoverflow.com/questions/5794414/using-multiple-criteria-in-subset-function-and-logical-operators

### Observed OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
       mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 OTUs: EdU+",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
    axis.text = element_text(size = 20),
    axis.title=element_text(size=24), legend.text=element_text(size=20),
    legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 OTUs: Whole",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 Shannon: EdU+",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 Shannon: Whole",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 Simpson E: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 Simpson E: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 Simpson: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F1 Simpson: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## F2 -------------------------------------------------------------------------

# subsetting notation from:
# https://stackoverflow.com/questions/5794414/using-multiple-criteria-in-subset-function-and-logical-operators

### Observed OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 OTUs: EdU+",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 OTUs: Whole",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 Shannon: EdU+",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 Shannon: Whole",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 Simpson E: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 Simpson E: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 Simpson: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage F2 Simpson: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## M1 -------------------------------------------------------------------------

# subsetting notation from:
# https://stackoverflow.com/questions/5794414/using-multiple-criteria-in-subset-function-and-logical-operators

### Observed OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 OTUs: EdU+",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 OTUs: Whole",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 Shannon: EdU+",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 Shannon: Whole",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 Simpson E: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 Simpson E: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 Simpson: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M1 Simpson: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## M2 -------------------------------------------------------------------------

# subsetting notation from:
# https://stackoverflow.com/questions/5794414/using-multiple-criteria-in-subset-function-and-logical-operators

### Observed OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 OTUs: EdU+",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 OTUs: Whole",
       color = "Health State",
       x = "Day",
       y = "Observed Features") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 Shannon: EdU+",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 Shannon: Whole",
       color = "Health State",
       x = "Day",
       y = "Shannon") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 Simpson E: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 Simpson E: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson E") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 Simpson: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(title = "Cage M2 Simpson: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Females -------------------------------------------------------------------------

## OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female OTUs: EdU+",
       color = "Health State",
       x = "Day",
       y = "Observed Features",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female OTUs: Whole",
       color = "Health State",
       x = "Day",
       y = "Observed Features",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female Shannon: EdU+",
       color = "Health State",
       x = "Day",
       y = "Shannon",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female Shannon: Whole",
       color = "Health State",
       x = "Day",
       y = "Shannon",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female Simpson E: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson E",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female Simpson E: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson E",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female Simpson: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "Female Simpson: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Males ---------------------------------------------------------------------------

## OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male OTUs: EdU+",
       color = "Health State",
       x = "Day",
       y = "Observed Features",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male OTUs: Whole",
       color = "Health State",
       x = "Day",
       y = "Observed Features",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male Shannon: EdU+",
       color = "Health State",
       x = "Day",
       y = "Shannon",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male Shannon: Whole",
       color = "Health State",
       x = "Day",
       y = "Shannon",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male Simpson E: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson E",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male Simpson E: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson E",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male Simpson: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(1, 2)) +
  labs(title = "Male Simpson: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## All cages -------------------------------------------------------------------------

## OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All OTUs: EdU+",
       color = "Health State",
       x = "Day",
       y = "Observed Features",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=day, y=observed_features, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All OTUs: Whole",
       color = "Health State",
       x = "Day",
       y = "Observed Features",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All Shannon: EdU+",
       color = "Health State",
       x = "Day",
       y = "Shannon",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=day, y=shannon_entropy, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All Shannon: Whole",
       color = "Health State",
       x = "Day",
       y = "Shannon",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All Simpson E: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson E",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=day, y=simpson_e, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All Simpson E: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson E",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All Simpson: EdU+",
       color = "Health State",
       x = "Day",
       y = "Simpson",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=day, y=simpson, group = cage, color = hs)) + 
  geom_line(linewidth = 1.5) +
  geom_point(aes(shape = cage), size = 4) +
  scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(title = "All Simpson: Whole",
       color = "Health State",
       x = "Day",
       y = "Simpson",
       shape = "Cage") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# Boxplots ---------------------------------------------------------------------------

# per health state

# make hs a factor & organize levels & labels

all_alpha$hs <- factor(all_alpha$hs, levels = c("Baseline", "Pre-symptomatic",
                                                "Symptomatic", "Recovery"),
                       labels = c("BA", "PS", "SY", "RE"))

## F1 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage F1 OTUs: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage F1 OTUs: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage F1 Shannon: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage F1 Shannon: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage F1 Simpson E: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage F1 Simpson E: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage F1 Simpson: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage F1 Simpson: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## F2 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage F2 OTUs: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage F2 OTUs: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage F2 Shannon: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage F2 Shannon: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage F2 Simpson E: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage F2 Simpson E: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage F2 Simpson: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "F2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage F2 Simpson: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## M1 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage M1 OTUs: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage M1 OTUs: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage M1 Shannon: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage M1 Shannon: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage M1 Simpson E: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage M1 Simpson E: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage M1 Simpson: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M1" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage M1 Simpson: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## M2 ----------------------------------------------------------------------------------

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage M2 OTUs: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Cage M2 OTUs: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage M2 Shannon: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Cage M2 Shannon: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage M2 Simpson E: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Cage M2 Simpson E: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage M2 Simpson: EdU+",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

all_alpha %>% ggplot(data = subset(all_alpha, subset = (cage == "M2" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Cage M2 Simpson: Whole",
       color = "Health States") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24), legend.position = "none")

## Females -------------------------------------------------------------------------

## OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Female OTUs: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Female OTUs: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Female Shannon: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Female Shannon: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Female Simpson E: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Female Simpson E: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Female Simpson: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Female" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Female Simpson: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Males ---------------------------------------------------------------------------

## OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Male OTUs: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Observed Features", 
       title = "Male OTUs: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Male Shannon: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Shannon", 
       title = "Male Shannon: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Male Simpson E: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Simpson E", 
       title = "Male Simpson E: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Male Simpson: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sex == "Male" & sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(1, 2)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "Male Simpson: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

## All cages -------------------------------------------------------------------------

## OTUs

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
  scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Observed Features", 
       #title = "All OTUs: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=hs, y=observed_features, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
   scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Observed Features", 
       #title = "All OTUs: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

## Shannon

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
   scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Shannon Index", 
       #title = "All Shannon: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=hs, y=shannon_entropy, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
   scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Shannon Index", 
       #title = "All Shannon: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

## Simpson E

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
   scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Simpson Evenness", 
       #title = "All Simpson E: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson_e, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
   scale_shape_manual(values = c(16, 1, 17, 2)) +
  labs(x = "Day", 
       y = "Simpson Evenness", 
       #title = "All Simpson E: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(#plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 20),
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=20))

## Simpson

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "EdU+")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
   scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "All Simpson: EdU+",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))

all_alpha %>% ggplot(data = subset(all_alpha, subset = (sf == "Whole")), 
                     mapping = aes(x=hs, y=simpson, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = cage), width = 0.2, size = 3) +
   scale_shape_manual(values = c(16, 17, 1, 2)) +
  labs(x = "Day", 
       y = "Simpson", 
       title = "All Simpson: Whole",
       fill = "Health States",
       shape = "Cage") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() + theme_bw() +
  theme(plot.title = element_text(size = 36, hjust = 0.5), 
        axis.text = element_text(size = 16),
        axis.title=element_text(size=24))
