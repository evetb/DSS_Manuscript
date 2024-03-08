#!/bin/bash

#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running dastool's Fasta_to_Contig2Bin script
# for dastool input
# for all bins created from the cage F2
# coassembly;
# READS was changed to "M2"
# to run the same script 
# for cage M2

module load singularity/3.8

cwd=/scratch/etb
READS=F2

cd $cwd/dastool

#maxbin2
singularity exec -B /home -B /project -B /scratch das_tool-1.1.4.img Fasta_to_Contig2Bin.sh -i /scratch/etb/maxbin2/$READS/bins -e fasta > /scratch/etb/dastool/scaffolds_$READS/maxbin_$READS.scaffolds2bin.tsv

#concoct
singularity exec -B /home -B /project -B /scratch das_tool-1.1.4.img Fasta_to_Contig2Bin.sh -i /scratch/etb/concoct/concoct_$READS/bins -e fa > /scratch/etb/dastool/scaffolds_$READS/concoct_$READS.scaffolds2bin.tsv

#metabat
singularity exec -B /home -B /project -B /scratch das_tool-1.1.4.img Fasta_to_Contig2Bin.sh -i /scratch/etb/metabat2/$READS/bins -e fa > /scratch/etb/dastool/scaffolds_$READS/metabat_$READS.scaffolds2bin.tsv

