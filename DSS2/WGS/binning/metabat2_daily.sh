#!/bin/bash
#SBATCH --account=def-corinnem
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running metabat on the per-day
# assemblies for cage F2;
# READS was changed to "M2"
# to run the same script 
# for cage M2

READS=F2
FASTA_DIR=/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/$READS
MEGA=/home/etb/projects/def-corinnem/etb/dss2_WGS/megahit/megahit_$READS-daily
cwd=/project/def-corinnem/etb/dss2_WGS

#cd $cwd/metabat2/$READS

module load StdEnv/2020 gcc/9.3.0 metabat/2.14

cd $FASTA_DIR

# base & sample are to get the appropriate names
# for each read, (which includes cage
# and day) based on their names
# from the Bowtie2-decontaminated reads 

for i in *.1.fastq.gz
do
        base="${i##*/}"
        sample=${base%.mouse_removed.1.fastq.gz}
        mkdir $cwd/metabat2/$READS/$sample-out
        cd $cwd/metabat2/$READS/$sample-out
        jgi_summarize_bam_contig_depths --outputDepth $sample-depth.txt $cwd/bowtie2/bowtie_align/$READS/$sample-bams/$sample-alignedReads.bam
        metabat -i $MEGA/$sample-megahit_out/$sample.contigs.fa -a $sample-depth.txt -o $cwd/metabat2/$READS/$sample-out/bins/$sample-bin -m 1500 -t 8
done
