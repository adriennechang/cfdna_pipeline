#!/bin/bash

## Save destination
DEST=$1

## Make directory at destination
mkdir ${DEST}/cfdna
mkdir ${DEST}/cfdna/results
mkdir ${DEST}/data

## Copy everything except results, references, software
cp -r prep_tables ${DEST}/cfdna/
cp -r workflow ${DEST}/cfdna/
cp README.txt ${DEST}/cfdna/
cp sequencing_prep.tsv ${DEST}/cfdna/
cp reference_methylomes.tsv ${DEST}/cfdna/
cp Snakefile ${DEST}/cfdna/
cp cfdna_v2.yml ${DEST}/cfdna/

find ${DEST}/cfdna -type d -exec chmod 755 {} \;
find ${DEST}/cfdna -type f -exec chmod 666 {} \;

find ${DEST}/cfdna/workflow/scripts/ -type d -exec chmod 755 {} \;
find ${DEST}/cfdna/workflow/scripts/ -type f -exec chmod 777 {}\;
