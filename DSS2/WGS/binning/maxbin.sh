#!/bin/bash
#SBATCH --account=def-corinnem
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running maxbin2 on coassembled
# reads from F2; READS was changed to
# "M2" to run the same script for
# cage M2

module load StdEnv/2020 gcc/9.3.0 maxbin/2.2.7

cwd=/project/def-corinnem/etb/dss2_WGS
READS=F2

cd $cwd/maxbin2/$READS

run_MaxBin.pl -contig $cwd/megahit/megahit_$READS/final.contigs.fa -abund $cwd/metabat2/$READS/$READS-depth.txt -out $cwd/maxbin2/$READS/bins/$READS-maxbin -min_contig_length 1500 -thread 16
