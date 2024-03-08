#!/bin/bash

# simple script for converting 
# cage F2 CONCOCT csv output into 
# tab separated contigs2bin file
# READS was changed to "M2"
# to run the same script for
# cage M2

# with the help of 
# https://stackoverflow.com/questions/7050717/bash-find-a-keyword-in-a-file-and-delete-its-line
# https://unix.stackexchange.com/questions/148114/how-to-add-words-to-an-existing-column
# https://stackoverflow.com/questions/3806874/how-to-merge-two-files-line-by-line-in-bash

cwd=/scratch/etb/
READS=F2
FASTA_DIR=$cwd/bowtie2/bowtie_remove_mouse/$READS

cd $FASTA_DIR

# base & sample are to get the appropriate names
# for each read, (which includes cage
# and day) based on their names
# from the Bowtie2-decontaminated reads 

for i in *.1.fastq.gz
do
        base="${i##*/}"
        sample=${base%.mouse_removed.1.fastq.gz}
        cd $cwd/concoct/concoct_daily_$READS/concoct_$sample
        cp clustering_gt1000.csv file
        sed -i '/id/d' file
        cut -f1 -d '.' file > temp_contigs
        cut -f2 -d ',' file > temp_bin
        paste temp_contigs temp_bin > temp_c2b
        awk -F'\t' -vOFS='\t' '{ $2= "concoct." $2}1' < temp_c2b > $cwd/dastool/dastool_daily/$READS/contigs/concoct_$sample.contigs2bin.tsv 
        rm file temp_contigs temp_bin temp_c2b
done
