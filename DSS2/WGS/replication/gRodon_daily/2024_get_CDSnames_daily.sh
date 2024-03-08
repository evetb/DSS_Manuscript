#!/bin/bash

# script to get CDS names
# for all proteins from each bin
# or all bins in a community,
# as needed for gRodon

for i in *.gff3 
do 
	BASE="${i##*/}" 
	BBASE=${BASE%.gff3}
	INPUT="$BBASE"
	BIN=$(echo $INPUT| cut -d'.' -f 1)
	sed -n '/##FASTA/q;p' $i | awk '$3=="CDS"' | awk '{print $9'} | awk 'gsub(";.*","")' | awk 'gsub("ID=","")' > CDS_names_$BIN.txt
done

