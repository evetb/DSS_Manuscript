#!/bin/bash
#SBATCH --time=2:0:0
#SBATCH --tasks=16
#SBATCH --mem=6G
#SBATCH --array=1-20
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# script to use bowtie2 to align
# reads to mouse genome, and then
# remove those reads which match
# to remove mouse (host) contamination
# of gut bacterial sequences

cwd=/project/def-corinnem/etb/dss2_WGS

module load StdEnv/2020 bowtie2/2.4.4

# NOTE: DIRECTORY HERE MUST BE THE ONE WITH samples.txt FILE
cd /project/def-corinnem/etb/dss2_WGS/bowtie2

mkdir -p bowtie_remove_mouse/logs

# NOTE: IT IS CRUCIAL THAT THE samples.txt FILE HAVE A HEADER "sample", THEN ONE SAMPLE PER LINE, AND NO EXTRA LINES/CONTENT
SAMPLE=( $(tail -n +2 samples.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') )

# Access & export mouse genome bowtie2 index from Compute Canada
index=/cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.GRCm38/genome/bowtie2_index/Mus_musculus.GRCm38
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6

# Note that the % in the bowtie2 output is a command to put the number 1 for the forward reads and the number 2 for the reverse reads
bowtie2 -x $index -1 $cwd/reads-trimmed/${SAMPLE}_forward_paired.fq.gz -2 $cwd/reads-trimmed/${SAMPLE}_reverse_paired.fq.gz --sensitive --un-conc-gz bowtie_remove_mouse/${SAMPLE}.mouse_removed.%.fastq.gz --threads ${SLURM_TASKS_PER_NODE} 1> /dev/null 2> bowtie_remove_mouse/logs/${SAMPLE}.mouse_removed.log
