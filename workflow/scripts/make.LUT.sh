#!/bin/bash

gitax=$1
taxonly=$2
taxname=$3
taxtax=$4

mkdir -p LUT;
cut -f3 $gitax | sort -u >  $taxonly
taxonkit lineage -n $taxonly | taxonkit reformat | \
	csvtk -H -t cut -f 1,3 | \
	csvtk pretty -t -s ","> $taxname;
taxonkit lineage $taxonly | \
	taxonkit reformat -t | \
	csvtk -H -t cut -f 1,4 | \
	csvtk -H -t sep -f 2 -s ';' -R | \
	awk '{ print $1, $7, $8 }' | \
	sort -u > $taxtax;
