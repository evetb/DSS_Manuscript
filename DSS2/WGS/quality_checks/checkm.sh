#!/bin/bash

#SBATCH --mem=5G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:0
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<eve.beauchemin@mail.mcgill.ca>

# running checkm on the bins
# from the co-assemblies per cage

# READS was changed to "M2"
# to run the same script for
# cage M2

cwd=/home/etb/projects/def-corinnem/etb
scr=/scratch/etb
READS=F2

module load StdEnv/2020 hmmer/3.2.1 prodigal/2.6.3 pplacer/1.1.alpha19

source $cwd/myvenv/bin/activate

#checkm lineage_wf -t 8 -x fa $scr/dastool/dastool_$READS/DAS_DASTool_bins $scr/checkm/checkm_$READS

checkm qa $scr/checkm/checkm_$READS/lineage.ms $scr/checkm/checkm_$READS -f $scr/checkm/checkm_$READS/qa_output1.tsv -t $SLURM_CPUS_PER_TASK --tab_table -o 1

checkm qa $scr/checkm/checkm_$READS/lineage.ms $scr/checkm/checkm_$READS -f $scr/checkm/checkm_$READS/qa_output2.tsv -t $SLURM_CPUS_PER_TASK --tab_table -o 2
