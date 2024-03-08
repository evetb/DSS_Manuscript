#!/bin/bash
#SBATCH --account=def-corinnem
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running maxbin2 on the bins
# created per day for cage F2;
# READS was changed to "M2"
# to run the same script for cage M2

READS=F2
cwd=/home/etb/projects/def-corinnem/etb/dss2_WGS
FASTA_DIR=$cwd/bowtie2/bowtie_remove_mouse/$READS
MEGA=$cwd/megahit/megahit_$READS-daily

module load StdEnv/2020 gcc/9.3.0 maxbin/2.2.7

cd $FASTA_DIR

# base & sample are to get the appropriate names
# for each read, (which includes cage
# and day) based on their names
# from the Bowtie2-decontaminated reads 

for i in *.1.fastq.gz
do
        base="${i##*/}"
        sample=${base%.mouse_removed.1.fastq.gz}
        cd $cwd/maxbin2/$READS
        mkdir $sample-bins
        run_MaxBin.pl -contig $MEGA/$sample-megahit_out/$sample.contigs.fa -abund $cwd/metabat2/$READS/$sample-out/$sample-depth.txt -out $cwd/maxbin2/$READS/$sample-bins/$sample-maxbin -min_contig_length 1500 -thread 16
done

# to run maxbin2, i first need metabat2 depth output
