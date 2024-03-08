#!/bin/bash

# running GRiD on the lab virtual machine
# looping through each SAM file
# this does all SAMS from one day at a time
# due to lack of storage space 

# READS, DAY and BIN were changed for each
# run as needed 
# all bins for one day, then all days for one cage
# and repeated for the other cage (cage F2 & M2)

source activate grid

READS=M2
DAY=M2D7_S54

SAMS=/home/ebeauchemin/$READS-day-bin-bams/$DAY

# ensure the output directory (-o) is empty
# thus you need a new directory for each run

# for loop

cd /home/ebeauchemin/grid_bins/$READS

for i in *.fa; do BASE="${i##*/}"; BIN=${BASE%.fa}; grid single -r $SAMS/$DAY-$BIN-bams -e sam -g /home/ebeauchemin/grid_bins/$READS/$BIN.fa -o /home/ebeauchemin/grid_out/$READS/$DAY/$DAY-$BIN -n 8; done

# one at a time

BIN=F2-bin.1 # change this each time

grid single -r $SAMS/$DAY-$BIN-bams -e sam -g /home/ebeauchemin/grid_bins/$READS/$BIN.fa -o /home/ebeauchemin/grid_out/$READS/$DAY/$DAY-$BIN -n 8