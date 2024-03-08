#!/bin/bash
#SBATCH --time=24:0:0
#SBATCH --mem=10G
#SBATCH --array=1-10
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# script for mapping the reads of each day
# to an index made of a single bin
# for every bin w/in one cage
# (with per-bin indices having already being made)
# this is needed for running GRiD
# READS was changed to "M2"
# to run the same script for
# cage M2

READS=F2
cwd=/home/etb/projects/def-corinnem/etb/dss2_WGS
FASTA_DIR=$cwd/bowtie2/bowtie_remove_mouse/$READS

cd $cwd/bowtie2
SAMPLE=( $(tail -n +2 samples_$READS.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') )

module load StdEnv/2020 gcc/9.3.0 bowtie2/2.4.4 samtools/1.13

#mkdir -p bowtie_align/$READS/logs

# you must create the index folder that bowtie will output to
mkdir -p $cwd/bowtie2/bowtie_align/$READS/index_bins

cd $cwd/grid_bins/$READS

# base and bin to get appropriate bin names
# which are used for the relevant bam files 
# and other files (input & output)

for i in *.fa
do
        base="${i##*/}"
        bin=${base%.fa}
        mkdir $cwd/bowtie2/bowtie_align/$READS/${SAMPLE}-$bin-bams
        R1s=$FASTA_DIR/${SAMPLE}.mouse_removed.1.fastq.gz
        R2s=$FASTA_DIR/${SAMPLE}.mouse_removed.2.fastq.gz
        bowtie2 -x $cwd/bowtie2/bowtie_align/$READS/index_bins/index_$bin -1 $R1s -2 $R2s -S $cwd/bowtie2/bowtie_align/$READS/${SAMPLE}-$bin-bams/${SAMPLE}-$bin-alignedReads.sam --reorder --no-unal --threads ${SLURM_TASKS_PER_NODE} 1> /dev/null 2> $cwd/bowtie2/bowtie_align/$READS/logs/${SAMPLE}-$bin-align.log
        samtools sort $cwd/bowtie2/bowtie_align/$READS/${SAMPLE}-$bin-bams/${SAMPLE}-$bin-alignedReads.sam -O bam -o $cwd/bowtie2/bowtie_align/$READS/${SAMPLE}-$bin-bams/${SAMPLE}-$bin-alignedReads.bam
        samtools index $cwd/bowtie2/bowtie_align/$READS/${SAMPLE}-$bin-bams/${SAMPLE}-$bin-alignedReads.bam
        rm $cwd/bowtie2/bowtie_align/$READS/${SAMPLE}-$bin-bams/${SAMPLE}-$bin-alignedReads.sam
done
