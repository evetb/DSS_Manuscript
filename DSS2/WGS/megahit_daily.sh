#!/bin/bash
#SBATCH --account=def-corinnem
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running megahit to assemble reads
# for each day separately, in each cage,
# being cage F2 and cage M2.
# READS was changed to "M2"
# to run the same script for cage M2

READS=F2
FASTA_DIR=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/$READS
ASSEM_DIR=/home/etb/projects/def-corinnem/etb/dss2_WGS/megahit

cd $ASSEM_DIR
mkdir megahit_$READS-daily

module load StdEnv/2020 megahit/1.2.9

cd $FASTA_DIR

for i in *.1.fastq.gz
do
        base="${i##*/}"
        sample=${base%.mouse_removed.1.fastq.gz}
        R1s=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/$READS/$sample.mouse_removed.1.fastq.gz
        R2s=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/$READS/$sample.mouse_removed.2.fastq.gz
        megahit -1 $R1s -2 $R2s -o $ASSEM_DIR/megahit_$READS-daily/$sample-megahit_out --out-prefix $sample -t 16
done
