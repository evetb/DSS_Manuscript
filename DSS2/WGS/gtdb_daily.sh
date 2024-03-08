#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --tasks=16
#SBATCH --mem=64G

# using GTDB-Tk to taxonmically identify
# all bins created from each time point
# (individually) from cage F2
# the READS variable was changed to "M2"
# to run the same script for cage M2

module load singularity/3.8 python/3.9.6 gcc/9.3.0 prodigal/2.6.3 hmmer/3.2.1 pplacer/1.1.alpha19 fastani/1.32 fasttree/2.1.11 mash/2.3

source /project/def-corinnem/etb/dss2_WGS/gtdbtk-env/bin/activate

BINS=/home/etb/projects/def-corinnem/etb/dss2_WGS/GTDB-Tk_bins
READS=F2
OUT=/home/etb/projects/def-corinnem/etb/dss2_WGS

cd /home/etb/projects/def-corinnem/etb/dss2_WGS

mkdir -p GTDB-Tk_out/$READS-daily

GTDBTK_DATA_PATH=/project/def-corinnem/databases/gtdb-r207v2/; export GTDBTK_DATA_PATH
gtdbtk classify_wf --extension fa --genome_dir $BINS/$READS-daily_bins/ --out_dir $OUT/GTDB-Tk_out/$READS-daily/ --cpus ${SLURM_TASKS_PER_NODE}
deactivate
