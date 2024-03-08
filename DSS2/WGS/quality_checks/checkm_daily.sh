#!/bin/bash

#SBATCH --mem=45G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:0
#SBATCH --array=1-10
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running checkm for the bins
# created from each day 
# (not co-asssembled)

# READS was changed to "M2"
# to run the same script for
# cage M2

cwd=/home/etb/projects/def-corinnem/etb
scr=/scratch/etb
READS=F2

SAMPLE=( $(tail -n +2 $cwd/dss2_WGS/bowtie2/samples_$READS.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') )

module load StdEnv/2020 hmmer/3.2.1 prodigal/2.6.3 pplacer/1.1.alpha19

source $cwd/myvenv/bin/activate

mkdir -p $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE}

checkm lineage_wf -t 8 -x fa $scr/dastool/dastool_daily/$READS/dastool_${SAMPLE}/${SAMPLE}-DASTool_bins $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE}

checkm qa $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE}/lineage.ms $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE} -f $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE}/qa_output1_${SAMPLE}.tsv -t $SLURM_CPUS_PER_TASK --tab_table -o 1

checkm qa $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE}/lineage.ms $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE} -f $scr/checkm/checkm_daily/$READS/checkm_${SAMPLE}/qa_output2_${SAMPLE}.tsv -t $SLURM_CPUS_PER_TASK --tab_table -o 2

