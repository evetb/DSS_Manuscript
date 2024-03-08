#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=10:00:0
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running concoct for the bins
# created for each individual day
# here doing cage F2, but READS
# was changed to "M2" to run
# the same script for cage M2

READS=F2
cwd=/scratch/etb
MEGA=/scratch/etb/megahit_$READS-daily

module load singularity/3.8

cd $cwd/concoct

mkdir concoct_daily_$READS

CCT=$cwd/concoct/concoct_daily_$READS

# check that path to bowtie2 .bam output is correct

for i in $cwd/bowtie2/bowtie_remove_mouse/$READS/*.1.fastq.gz
do
	base="${i##*/}"
	sample=${base%.mouse_removed.1.fastq.gz}
	mkdir $CCT/concoct_$sample
	cd $CCT/concoct_$sample
	singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img cut_up_fasta.py $MEGA/$sample-megahit_out/$sample.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
	singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img concoct_coverage_table.py $CCT/concoct_$sample/contigs_10K.bed $cwd/bowtie2/bowtie_align/$READS/$sample-bams/$sample-alignedReads.bam > $CCT/concoct_$sample/coverage_table.tsv
	singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img concoct --composition_file $CCT/concoct_$sample/contigs_10K.fa --coverage_file $CCT/concoct_$sample/coverage_table.tsv -t 2 -b $CCT/concoct_$sample
	singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img merge_cutup_clustering.py $CCT/concoct_$sample/clustering_gt1000.csv > $CCT/concoct_$sample/clustering_merged.csv
	mkdir $CCT/concoct_$sample/bins
	singularity exec -B /home -B /project -B /scratch $cwd/concoct/concoct-1.1.0.img extract_fasta_bins.py $MEGA/$sample-megahit_out/$sample.contigs.fa $CCT/concoct_$sample/clustering_merged.csv --output_path $CCT/concoct_$sample/bins
done
