#!/bin/bash
#SBATCH --time=2:0:0
#SBATCH --tasks=8
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# for storing sam files as indexed bams

module load StdEnv/2020 gcc/9.3.0 samtools/1.13

samtools sort F2-alignedReads.sam -O bam -o F2-alignedReads.bam
samtools index F2-alignedReads.bam

samtools sort M2-alignedReads.sam -O bam -o M2-alignedReads.bam
samtools index M2-alignedReads.bam

# convert bam file back to sam file
for file in ./*.bam; do echo $file; samtools view -h $file > ${file/.bam/.sam}; done

# READS was changed to "F2"
# to run the same script for
# cage F2

READS=M2
cd /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/$READS-assembly_bams

for file in ./*.bam; do echo $file; samtools view -h $file > ${file/.bam/.sam}; done
