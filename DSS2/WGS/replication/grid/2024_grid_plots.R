# to make plots of replication rates
# per genome for GRiD output

# first need to import all the data.
# each replication rate for each genome for each day
# is in its own .txt file

# Libraries ---------------

library(tidyverse)
library(data.table)
library(gridExtra)

# data manipulation -------------

## F2 ----------------------------

# reading in one file (as a test)
F2D1_S48_bin.1 <- read_delim("F2D1_S48-F2-bin.1-alignedReads.GRiD.tsv")

# now i want to read them ALL in
# importing all .tsv files in all subdirectories 
# adding each one as a new row to a tibble

# using rbindlist from data.table, as seen here:
# https://stackoverflow.com/questions/39377370/bind-rows-of-different-data-types

# Set the working directory to the 
# parent directory of all the subdirectories
# which contain the GRiD values

setwd("")

# Create an empty tibble to hold all the data
all_F2_grid <- tibble()

# Use the list.files function to find all the .tsv files in all subdirectories
F2_grid_file_paths <- list.files(recursive = TRUE, pattern = "\\.tsv$", full.names = TRUE)

# Loop over the file paths and read each .tsv file into a tibble using read_tsv from dplyr
for (files in F2_grid_file_paths) {
  file_data <- read_tsv(files)
  all_F2_grid <- rbindlist((list(all_F2_grid, file_data)))
}

# View the final tibble with all the data
all_F2_grid

# now i want to add a "day" column and a "bin" column 

# first separating using "-"
# and getting rid of the alignedReads.GRiD part
# to have Sample, Cage, Bin, and Tail

all_F2_grid_sep <- all_F2_grid %>% 
  separate(Sample, into = c("Sample", "Cage", "Bin", "Tail"), sep = "-")

# removing Tail, which is the alignedReads.GRiD part
all_F2_grid_sep <- all_F2_grid_sep %>% select(-Tail)

# get Day column based on Sample name
all_F2_grid_sep <- all_F2_grid_sep %>% 
  mutate(Day = 
           case_when(Sample == "F2D1_S48" ~ 1,
                     Sample == "F2D5_S53" ~ 5,
                     Sample == "F2D7_S1" ~ 7,
                     Sample == "F2D8_S1" ~ 8,
                     Sample == "F2D11_S3" ~ 11,
                     Sample == "F2D12_S4" ~ 12,
                     Sample == "F2D13_S52" ~ 13,
                     Sample == "F2D14_S42" ~ 14,
                     Sample == "F2D15_S47" ~ 15,
                     Sample == "F2D18_S55" ~ 18))

# adding in health states, using cut()

all_F2_grid_sep$hs <- cut(all_F2_grid_sep$Day, c(-Inf, 7, 11, 15, Inf), 
                        c("Baseline", "Pre-symptomatic",
                          "Symptomatic", "Recovery"))

# exporting as csv to easily combine WGS & GRiD analyses

write.csv(all_F2_grid_sep, "all_F2_grid_sep.csv")

# finding all cases of GRiD >= 3
# which does not make biological sense
# https://sparkbyexamples.com/r-programming/r-select-rows-by-condition/

grid_over3_F2 <- all_F2_grid_sep[all_F2_grid_sep$GRiD >= 3,]
unique(grid_over3_F2$Bin)

# bin.138, bin.46, bin.112, bin.34, bin.147, bin.75, bin.47 

# making filtered, wide table 
# to have publication-friendly table
# of bin GRiD values per day
# for those bins with GRiD values <=3

all_F2_grid_sep.wide <- all_F2_grid_sep %>% 
  subset(GRiD <= 3) %>% # only keep GRiD values <= 3
  select(Bin, Day, GRiD) %>% # only selecting relevant columns
  pivot_wider(names_from = Day, values_from = GRiD) %>% # changing from long to wide format
  na.omit # removing rows w/NAs (which will be all the ones with GRiD values > 3)
  
write.csv(all_F2_grid_sep.wide, "all_F2_grid_sep.wide.csv")

## M2 --------------------------------

# Set the working directory to the 
# parent directory of all the subdirectories
# which contain the GRiD values

setwd("")

# Create an empty tibble to hold all the data
all_M2_grid <- tibble()

# Use the list.files function to find all the .tsv files in all subdirectories
M2_grid_file_paths <- list.files(recursive = TRUE, pattern = "\\.txt$", full.names = TRUE)

# Loop over the file paths and read each .tsv file into a tibble using read_tsv from dplyr
for (files in M2_grid_file_paths) {
  file_data <- read_tsv(files)
  all_M2_grid <- rbindlist((list(all_M2_grid, file_data)))
}

# View the final tibble with all the data
all_M2_grid

# now i want to add a "day" column and a "bin" column 

# first separating using "-"
# and getting rid of the alignedReads.GRiD part
# to have Sample, Cage, Bin, and Tail

all_M2_grid_sep <- all_M2_grid %>% 
  separate(Sample, into = c("Sample", "Cage", "Bin", "Tail"), sep = "-")

# removing Tail, which is the alignedReads.GRiD part
all_M2_grid_sep <- all_M2_grid_sep %>% select(-Tail)

# get Day column based on Sample name
all_M2_grid_sep <- all_M2_grid_sep %>% 
  mutate(Day = 
           case_when(Sample == "M2D1_S43" ~ 1,
                     Sample == "M2D5_S44" ~ 5,
                     Sample == "M2D7_S54" ~ 7,
                     Sample == "M2D8_S2" ~ 8,
                     Sample == "M2D11_S5" ~ 11,
                     Sample == "M2D12_S6" ~ 12,
                     Sample == "M2D13_S49" ~ 13,
                     Sample == "M2D14_S46" ~ 14,
                     Sample == "M2D15_S45" ~ 15,
                     Sample == "M2D18_S50" ~ 18))

# adding in health states, using cut()

all_M2_grid_sep$hs <- cut(all_M2_grid_sep$Day, c(-Inf, 7, 11, 15, Inf), 
                          c("Baseline", "Pre-symptomatic",
                            "Symptomatic", "Recovery"))

# exporting as csv to easily combine WGS & GRiD analyses

write.csv(all_M2_grid_sep, "all_M2_grid_sep.csv")

# finding all cases of GRiD >= 3
# which does not make biological sense
# https://sparkbyexamples.com/r-programming/r-select-rows-by-condition/

grid_over3_M2 <- all_M2_grid_sep[all_M2_grid_sep$GRiD >= 3,]
unique(grid_over3_M2$Bin)

# "bin.114"        "bin.41"         "maxbin.114_sub" "bin.96"         "bin.97"

# making filtered, wide table 
# to have publication-friendly table
# of bin GRiD values per day
# for those bins with GRiD values <=3

all_M2_grid_sep.wide <- all_M2_grid_sep %>% 
  subset(GRiD <= 3) %>% # only keep GRiD values <= 3
  select(Bin, Day, GRiD) %>% # only selecting relevant columns
  pivot_wider(names_from = Day, values_from = GRiD) %>% # changing from long to wide format
  na.omit # removing rows w/NAs (which will be all the ones with GRiD values > 3)

write.csv(all_M2_grid_sep.wide, "all_M2_grid_sep.wide.csv")

# plotting -------------------

# for each bin, plot its GRiD values over time
# colouring by health state

# next want to associate taxonomy with bin
# this was done manually, by looking at
# an excel sheet of GTDB-Tk taxonomy names 
# & their associated bins

## Box plots ----------------------------------------

# only for taxa highlighted in the paper

# Akkermansia muciniphila
# Duncaniella sp910589485
# Bifidobacterium globosum
# Muribaculaceae CAG-485 sp002362485

### F2 --------------------------------------

ggplot(subset(all_F2_grid_sep, Bin %in% "maxbin.006_sub"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Akkermansia muciniphila"), " - F2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.122"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Duncaniella"), " sp910589485 - F2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "maxbin.016_sub"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Bifidobacterium globosum"), " - F2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.94"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Muribaculaceae"), " CAG-485 sp002362485 - F2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### M2 --------------------------------------
ggplot(subset(all_M2_grid_sep, Bin %in% "bin.3"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Akkermansia muciniphila"), " - M2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.126_sub"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Duncaniella"), " sp910589485 - M2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "maxbin.010"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Bifidobacterium globosum"), " - M2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.135"), 
       aes(x=Day, y=GRiD, fill = hs)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, size = 2, shape = 16) +
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Muribaculaceae"), " CAG-485 sp002362485 - M2"))),
       fill = "Health State") +
  scale_fill_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

## Line plots -------------------------------------------

### M2 ---------------------

# 43 bins

# "bin.10"         "bin.101"        "bin.110"        "bin.114"        "bin.125"

# missing many days
ggplot(subset(all_M2_grid_sep, Bin %in% "bin.10"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceace"), " MGBC164599 sp910575625 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.101"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acutalibacter"), " sp910587995 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.110"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " MD308 sp010206225 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.114"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " 1XD8-76 sp910573755 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.125"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Eubacterium"), " sp910584245 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "bin.126_sub"    "bin.129"        "bin.13"         "bin.135"        "bin.137"       

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.126_sub"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Duncaniella"), " sp910589485 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.129"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acutalibacter"), " sp009936035 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.13"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Choladocola"), " sp009774135 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.135"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Muribaculaceae"), " CAG-485 sp002362485 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.137"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Paramuribaculum"), " sp910579675 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "bin.140"        "bin.152"        "bin.155"        "bin.166"        "bin.20"

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.140"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " 1XD8-76 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.152"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Kineothrix"), " sp000403275 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.155"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " RGIG7193 sp910586125 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.166"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acetatifactor"), " sp910584235 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.20"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " UBA3282 sp910579735 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "bin.21"         "bin.25"         "bin.27"         "bin.3"          "bin.30"

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.21"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " VSOB01 sp910587635 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.25"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " UBA3282 sp009774585 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.27"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Ruminococcaceae"), " CAG-115 sp910587465 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.3"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Akkermansia muciniphila"), " - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.30"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acetatifactor"), " sp910579755 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "bin.35_sub"     "bin.40"         "bin.41"         "bin.46"         "bin.47"  

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.35_sub"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Adlercreutzia muris"), " - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.40"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Ruminiclostridium"), " sp910585505 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.41"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("bin.41"), " - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.46"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Caccovivens"), " sp - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.47"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " 14-2 sp009774455 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "bin.56_sub"     "bin.60"         "bin.71"         "bin.73"         "bin.76"   

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.56_sub"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Eggerthellaceae"), " D16-63 sp910588095 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.60"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " UBA3282 sp910577735 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.71"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " UBA3282 sp009774655 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.73"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Caccovicinus"), " sp910575565 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.76"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " COE1 sp009774375 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "bin.79"         "bin.80"         "bin.81"         "bin.82"         "bin.89"     

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.79"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " CAG-95 sp910587295 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.80"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Butyribacter"), " sp009774235 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.81"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " An181 sp910585545 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.82"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Anaerotignum"), " sp910576545 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.89"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acutalibacter muris"), " - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "bin.9"          "bin.94"         "bin.95"         "bin.96"         "bin.97"  

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.9"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " 14-2 sp910585625 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.94"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Bacilli"), " UBA5026 sp910586425- M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.95"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " CAG-95 sp910579425 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.96"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Choladocola"), " sp910575445 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "bin.97"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " COE1 sp009774345 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# "maxbin.003"     "maxbin.010"     "maxbin.114_sub"

ggplot(subset(all_M2_grid_sep, Bin %in% "maxbin.003"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Muribaculaceae"), " CAG-873 sp910577315 - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_M2_grid_sep, Bin %in% "maxbin.010"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Bifidobacterium globosum"), " - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# values too high (10-35)
ggplot(subset(all_M2_grid_sep, Bin %in% "maxbin.114_sub"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("maxbin.114_sub"), " - M2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

### F2 ---------------------

# Set the working directory to the parent directory of all the subdirectories
setwd("E:/DESKTOP-UFR1DAD/Documents/McGill/DSS_Manuscript/DSS2/WGS/grid/grid_out/F2")

# 32 bins

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.1"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Scatovivens"), " sp. - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.109"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " CAG-95 sp910587295 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.112"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Pelethomonas"), " sp910587645 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.114"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acutalibacter"), " sp009917525 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.118"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Sporofaciens"), " sp910585725 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.122"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Duncaniella"), " sp910589485 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.127"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " VSOB01 sp910587635 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

# very high values
ggplot(subset(all_F2_grid_sep, Bin %in% "bin.138"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Eubacterium G"), " sp - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.14"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " UBA3282 sp910575775 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.145"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " 1XD8-76 sp. - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.147"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Avoscillospira"), " sp - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.29"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Christensenellales"), " UBA3700 RACS-045 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.32"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " COE1 sp009774345 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.34"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Christensenellales"), " UBA3700 MGBC161649 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.35"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Christensenellales"), " UBA3700 MGBC161649 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.38"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acetatifactor"), " sp003612485 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.4"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Ruminococcaceae"), " CAG-115 sp910587465 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.45"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Eggerthellaceae"), " D16-63 sp910588095- F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.46"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acetatifactor"), " sp910585615 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.47"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Enterenecus"), " sp910587255 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.5"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acutalibacter"), " sp009936035 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.52"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Adlercreutzia muris"), " - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.55"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Duncaniella"), " sp910576785 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.7"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Acutalibacter"), " sp910587995 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.71"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Gallimonas"), " sp910585595 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.75"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lawsonibacter"), " sp910577805 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.79"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " 14-2 sp910585625 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.8"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lawsonibacter"), " sp910588635 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.9"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Lachnospiraceae"), " UBA3282 sp009774585 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "bin.94"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Muribaculaceae"), " CAG-485 sp002362485 - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "maxbin.006_sub"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Akkermansia muciniphila"), " - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))

ggplot(subset(all_F2_grid_sep, Bin %in% "maxbin.016_sub"), 
       aes(x=Day, y=GRiD, group = Bin, color = hs)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) + 
  labs(x = "Day", y = "GRiD value", 
       title =  (expression(paste(italic("Bifidobacterium globosum"), " - F2"))),
       color = "Health State") +
  scale_color_manual(values = c("#66CCEE", "#CCBB44", "#EE6677", "#228833")) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(size = 24, hjust = 0.5), axis.text = element_text(size = 20),
        axis.title=element_text(size=24), legend.text=element_text(size=20),
        legend.title=element_text(size=24))