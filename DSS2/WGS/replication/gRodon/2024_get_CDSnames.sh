#!/bin/bash

# script to get CDS names
# from names of each bin
# needed for gRodon input

for i in *.gff3 
do 
	BASE="${i##*/}" 
	BBASE=${BASE%.gff3}
	INPUT="$BBASE"
	BIN=$(echo $INPUT| cut -d'_' -f 3)
	sed -n '/##FASTA/q;p' $i | awk '$3=="CDS"' | awk '{print $9'} | awk 'gsub(";.*","")' | awk 'gsub("ID=","")' > CDS_names_$BIN.txt
done
