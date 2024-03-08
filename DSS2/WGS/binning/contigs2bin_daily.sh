#!/bin/bash

#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# this is needed for dastool input
# running dastool's Fasta_to_Contig2Bin script
# on the bins created from maxbin2 and metabat
# for cage F2; READS was changed to "M2"
# to run this same script for cage M2

module load singularity/3.8

cwd=/scratch/etb
READS=F2
FASTA_DIR=$cwd/bowtie2/bowtie_remove_mouse/$READS

cd $FASTA_DIR

# base & sample are to get the appropriate names
# for each read, (which includes cage
# and day) based on their names
# from the Bowtie2-decontaminated reads 

for i in *.1.fastq.gz
do
        base="${i##*/}"
        sample=${base%.mouse_removed.1.fastq.gz}
        cd $cwd/dastool
        #maxbin2
        singularity exec -B /home -B /project -B /scratch das_tool-1.1.4.img Fasta_to_Contig2Bin.sh -i $cwd/maxbin2/$READS/$sample-bins -e fasta > $cwd/dastool/dastool_daily/$READS/contigs/maxbin_$sample.contigs2bin.tsv
        #metabat
        singularity exec -B /home -B /project -B /scratch das_tool-1.1.4.img Fasta_to_Contig2Bin.sh -i $cwd/metabat2/$READS/$sample-out/bins -e fa > $cwd/dastool/dastool_daily/$READS/contigs/metabat_$sample.contigs2bin.tsv
done

