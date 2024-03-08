#!/bin/bash

#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running fastqc from a virtual environment
# on all of the raw reads for DSS #2

cd /home/etb/projects/def-corinnem/etb/dss2_WGS/reads

source /home/etb/projects/def-corinnem/etb/myvenv/bin/activate

module load fastqc/0.11.9

# command from https://sbc.shef.ac.uk/workshops/2020-02-12-command-line/read-quality.nb.html
fastqc *.fastq.gz -o /home/etb/projects/def-corinnem/etb/dss2_WGS/fastqc_output
