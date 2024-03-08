#!/bin/bash
#SBATCH --time=10:0:0
#SBATCH --tasks=16
#SBATCH --mem=10G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# script to align each day's reads 
# to each day's assemblies

# READS was changed to "M2"
# to run the same script for
# cage M2

# looping protocol from https://gist.github.com/meren/5c8616959f3eea1f0632b50e2f02fb1e
# naming protocol from https://www.biostars.org/p/467896/

READS=F2
cwd=/home/etb/projects/def-corinnem/etb/dss2_WGS
FASTA_DIR=$cwd/bowtie2/bowtie_remove_mouse/$READS
MEGA=$cwd/megahit/megahit_$READS-daily

module load StdEnv/2020 gcc/9.3.0 bowtie2/2.4.4 samtools/1.13

#mkdir -p bowtie_align/$READS/logs

# you must create the index folder that bowtie will output to
mkdir -p $cwd/bowtie2/bowtie_align/$READS/index_daily

cd $FASTA_DIR

# base and sample are used here to get appropriate names
# for each sample, to be used for file input/output

for i in *.1.fastq.gz
do
        base="${i##*/}"
        sample=${base%.mouse_removed.1.fastq.gz}
        mkdir $cwd/bowtie2/bowtie_align/$READS/$sample-bams
        bowtie2-build $MEGA/$sample-megahit_out/$sample.contigs.fa $cwd/bowtie2/bowtie_align/$READS/index_daily/index_$sample
        R1s=$FASTA_DIR/$sample.mouse_removed.1.fastq.gz
        R2s=$FASTA_DIR/$sample.mouse_removed.2.fastq.gz
        bowtie2 -x $cwd/bowtie2/bowtie_align/$READS/index_daily/index_$sample -1 $R1s -2 $R2s -S $cwd/bowtie2/bowtie_align/$READS/$sample-bams/$sample-alignedReads.sam --reorder --no-unal --threads ${SLURM_TASKS_PER_NODE} 1> /dev/null 2> $cwd/bowtie2/bowtie_align/$READS/logs/$sample-align.log
        samtools sort $cwd/bowtie2/bowtie_align/$READS/$sample-bams/$sample-alignedReads.sam -O bam -o $cwd/bowtie2/bowtie_align/$READS/$sample-bams/$sample-alignedReads.bam
        samtools index $cwd/bowtie2/bowtie_align/$READS/$sample-bams/$sample-alignedReads.bam
        rm $cwd/bowtie2/bowtie_align/$READS/$sample-bams/$sample-alignedReads.sam
done
