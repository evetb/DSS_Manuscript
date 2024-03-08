#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=10:00:0
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running concoct on the cage F2 co-assembly
# changed READS to "M2" to run the same
# script for cage M2

module load singularity/3.8

cwd=/scratch/etb
READS=F2

cd $cwd/concoct

#mkdir concoct_$READS

CCT=$cwd/concoct/concoct_$READS

cd $CCT


singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img cut_up_fasta.py $cwd/megahit/megahit_$READS/final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img concoct_coverage_table.py $CCT/contigs_10K.bed $cwd/bowtie2/$READS-alignedReads.bam > $CCT/coverage_table.tsv
singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img concoct --composition_file $CCT/contigs_10K.fa --coverage_file $CCT/coverage_table.tsv -t 2 -b $CCT
singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img merge_cutup_clustering.py $CCT/clustering_gt1000.csv > $CCT/clustering_merged.csv

mkdir $CCT/bins
singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img extract_fasta_bins.py $cwd/megahit/megahit_$READS/final.contigs.fa $CCT/clustering_merged.csv --output_path $CCT/bins
