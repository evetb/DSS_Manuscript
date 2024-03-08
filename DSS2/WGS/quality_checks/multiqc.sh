#!/bin/bash

#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=5:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running multiqc on the fastqc output
# from the reads of the DSS #2 experiment

cd /home/etb/projects/def-corinnem/etb/dss2_WGS/fastqc_output

source /home/etb/projects/def-corinnem/etb/myvenv/bin/activate

multiqc . -o /home/etb/projects/def-corinnem/etb/dss2_WGS/multiqc_output
