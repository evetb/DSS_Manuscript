#!/bin/bash
#SBATCH --time=24:0:0
#SBATCH --tasks=8
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# making bowtie indices for each bin
# for cage F2, as needed to run bowtie2 
# to  map the reads of each day
# to an index made of a single bin
# for every bin w/in one cage
# as needed to run GRiD.
# READS was changed to "M2"
# to run the same script for
# cage M2

READS=F2
cwd=/home/etb/projects/def-corinnem/etb/dss2_WGS

module load StdEnv/2020 gcc/9.3.0 bowtie2/2.4.4 samtools/1.13

# you must create the index folder that bowtie will output to
mkdir -p $cwd/bowtie2/bowtie_align/$READS/index_bins

cd $cwd/grid_bins/$READS

# base and bin to get appropriate bin names
# which are used for the output files

for i in *.fa
do
        base="${i##*/}"
        bin=${base%.fa}
        bowtie2-build $i $cwd/bowtie2/bowtie_align/$READS/index_bins/index_$bin
done