#!/bin/bash
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running dastool on all the bins
# created individually (per day)
# for cage F2; READS was changed
# to "M2" to run the same script
# for cage M2

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
        mkdir -p $cwd/dastool/dastool_daily/$READS/dastool_$sample
        singularity exec -B /home -B /project -B /scratch das_tool-1.1.4.img DAS_Tool -i $cwd/dastool/dastool_daily/$READS/contigs/metabat_$sample.contigs2bin.tsv,$cwd/dastool/dastool_daily/$READS/contigs/maxbin_$sample.contigs2bin.tsv,$cwd/dastool/dastool_daily/$READS/contigs/concoct_$sample.contigs2bin.tsv -c $cwd/megahit_$READS-daily/$sample-megahit_out/$sample.contigs.fa -o $cwd/dastool/dastool_daily/$READS/dastool_$sample/$sample --write_bins --write_bin_evals --write_unbinned -t 8
done
