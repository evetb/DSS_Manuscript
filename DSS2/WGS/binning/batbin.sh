#!/bin/bash
#SBATCH --account=def-corinnem
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running metabat on the co-assembly
# from cage F2; changed READS to "M2"
# to run the same script for cage M2

cwd=/project/def-corinnem/etb/dss2_WGS
READS=F2

cd $cwd/metabat2/$READS

module load StdEnv/2020 gcc/9.3.0 metabat/2.14

jgi_summarize_bam_contig_depths --outputDepth $READS-depth.txt $cwd/bowtie2/$READS-alignedReads.bam

metabat -i $cwd/megahit/megahit_$READS/final.contigs.fa -a $READS-depth.txt -o $cwd/metabat2/$READS/bins/$READS-bin -m 1500 -t 8
