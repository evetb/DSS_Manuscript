#!/bin/bash
#SBATCH --time=2:0:0
#SBATCH --tasks=16
#SBATCH --mem=6G
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# script to align reads of each cage 
# to its appropriate co-assembly

# READS was changed to "M2"
# to run the same script for
# cage M2

cwd=/project/def-corinnem/etb/dss2_WGS/megahit/megahit_F2
FASTA_DIR=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/F2
READS=F2
cd $FASTA_DIR

R1s=$(ls $FASTA_DIR/*.1.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')
R2s=$(ls $FASTA_DIR/*.2.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')

cd /project/def-corinnem/etb/dss2_WGS/bowtie2

module load StdEnv/2020 bowtie2/2.4.4

#mkdir -p bowtie_align/$READS/logs

# you must create the index folder that bowtie will output to
#mkdir -p bowtie_align/$READS/index

# the "index_$READS" gives the prefix "index_$READS" to each output file
#bowtie2-build $cwd/final.contigs.fa bowtie_align/$READS/index/index_$READS

bowtie2 -x bowtie_align/$READS/index/index_$READS -1 $R1s -2 $R2s -S $READS-alignedReads.sam --reorder --no-unal --threads ${SLURM_TASKS_PER_NODE} 1> /dev/null 2> bowtie_align/$READS/logs/$READS-align.log
