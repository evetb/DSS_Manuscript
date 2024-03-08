#!/bin/bash
#SBATCH --account=def-corinnem
#SBATCH --mem=50G
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running megahit to co-assemble
# all reads from all days within a single cage
# here doing cage F2

FASTA_DIR='/home/etb/projects/def-corinnem/etb/dss2_WGS/bowtie2/bowtie_remove_mouse/F2'
ASSEM_DIR='/home/etb/projects/def-corinnem/etb/dss2_WGS/megahit'

cd $FASTA_DIR

# from https://merenlab.org/tutorials/assembly-based-metagenomics/#co-assembly
# and https://www.protocols.io/view/mg-hw4-co-assembly-using-megahit-dm6gp3x1vzpn/v3
R1s=$(ls ./*.1.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')
R2s=$(ls ./*.2.fastq.gz | python -c 'import sys; print(",".join([x.strip() for x in sys.stdin.readlines()]))')

module load StdEnv/2020 megahit/1.2.9

megahit -1 $R1s -2 $R2s -o $ASSEM_DIR/megahit_F2 -t 16
