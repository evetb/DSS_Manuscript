#!/bin/bash
#SBATCH --time=2:0:0
#SBATCH --tasks=1
#SBATCH --mem=2G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# getting coverage information for 
# the bins that were run with GTDB-Tk
# for taxonomic assignment
# this is part of the process of 
# getting relative abundances of these bins

module load StdEnv/2020 gcc/9.3.0 samtools/1.13

cwd=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2
READS=M2
BAMS=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/GTDB-Tk_bams

cd $cwd

SAMPLES=( $(tail -n +2 samples_$READS.txt ) )

echo sample$'\t'contig$'\t'startpos$'\t'endpos$'\t'numreads$'\t'covbases$'\t'coverage$'\t'meandepth$'\t'meanbaseq$'\t'meanmapq > $BAMS/bin-coverage_$READS.tsv

for SAMPLE in "${SAMPLES[@]}"
do
  samtools coverage -H $BAMS/${SAMPLE}-alignedReads.bam | awk -v sample=${SAMPLE} '{print sample"\t"$0}' >> $BAMS/bin-coverage_$READS.tsv
done
