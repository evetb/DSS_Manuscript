#!/bin/bash
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running dastool
# for cage F2;
# READS was changed to 
# "M2" to run the sames script
# for cage M2

module load singularity/3.8

cwd=/scratch/etb
READS=F2

cd $cwd/dastool

singularity exec -B /home -B /project -B /scratch das_tool-1.1.4.img DAS_Tool -i $cwd/dastool/scaffolds_$READS/metabat_$READS.scaffolds2bin.tsv,$cwd/dastool/scaffolds_$READS/maxbin_$READS.scaffolds2bin.tsv,$cwd/dastool/scaffolds_$READS/cleaned_concoct_$READS.contigs2bin.tsv -c $cwd/megahit/megahit_$READS/final.contigs.fa -o DAS --write_bins --write_bin_evals --write_unbinned -t 8
