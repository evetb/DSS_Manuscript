#!/bin/bash
#SBATCH --time 0:30:0
#SBATCH --tasks=1
#SBATCH --mem=1G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# script to get bam coverage information 
# from bams of .ffn alignments (from prokka output)
# as part of the process of getting depth of 
# coverage information per gene 
# (needed for accurate community-wide 
# minimal doubling time with gRodon)

# using a for loop to go thorugh all bam files
# READS was changed to "M2"
# to run the same script for
# cage M2

cwd=/home/etb/projects/def-corinnem/etb/dss2_WGS/bamcov
READS=F2
BAMS=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/ffn_bams

cd $cwd

for bam in $BAMS/*.bam
do 
	base="${bam##*/}"
	sample=${base%-alignedReads.bam}
	./bamcov -H $bam >> $sample.tsv
	sed -i '1 i\#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq' $sample.tsv
done
