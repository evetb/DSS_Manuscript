#!/bin/bash
#SBATCH --time=10:0:0
#SBATCH --tasks=16
#SBATCH --mem=15G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# aligning reads to M2 .ffn files, as part of 
# the process of getting depth of coverage 
# information per gene (needed for accurate 
# community-wide minimal doubling times with gRodon)

# READS was changed to "F2"
# to run the same script for
# cage F2

# looping protocol from https://gist.github.com/meren/5c8616959f3eea1f0632b50e2f02fb1e
# naming protocol from https://www.biostars.org/p/467896/

READS=M2
FASTA_DIR=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/$READS
FFN=/home/etb/projects/def-corinnem/etb/dss2_WGS/prokka/$READS-ffns

cd /project/def-corinnem/etb/dss2_WGS/bowtie2

module load StdEnv/2020 gcc/9.3.0 bowtie2/2.4.4 samtools/1.13

mkdir -p bowtie_align/$READS/logs_ffns

# you must create the index folder that bowtie will output to
mkdir -p bowtie_align/$READS/index_ffns

# the "index_$READS" gives the prefix "index_$READS" to each output file
bowtie2-build $FFN/all_$READS-ffns.fa bowtie_align/$READS/index_ffns/index_$READS

cd $FASTA_DIR

# base and sample are used here to get appropriate names
# for each sample, to be used for file input/output

for i in *.1.fastq.gz
do
        base="${i##*/}"
        sample=${base%.mouse_removed.1.fastq.gz}
        R1s=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/$READS/$sample.mouse_removed.1.fastq.gz
        R2s=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/$READS/$sample.mouse_removed.2.fastq.gz
        bowtie2 -x /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/index_ffns/index_$READS -1 $R1s -2 $R2s -S /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/ffn_bams/$sample-alignedReads.sam --reorder --no-unal --threads ${SLURM_TASKS_PER_NODE} 1> /dev/null 2> /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/logs_ffns/$sample-align.log
        samtools sort /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/ffn_bams/$sample-alignedReads.sam -O bam -o /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/ffn_bams/$sample-alignedReads.bam
        samtools index /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/ffn_bams/$sample-alignedReads.bam
        rm /home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_align/$READS/ffn_bams/$sample-alignedReads.sam
done
