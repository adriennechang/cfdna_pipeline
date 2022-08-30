#!/bin/bash

## Save destination
DEST=$1

## Make directory at destination
mkdir ${DEST}/cfDNA_pipeline
mkdir ${DEST}/cfDNA_pipeline/results
mkdir ${DEST}/data

## Copy everything except results, references, software
cp -r prep_tables ${DEST}/cfDNA_pipeline/
cp -r workflow ${DEST}/cfDNA_pipeline/
cp README.md ${DEST}/cfDNA_pipeline/
cp add_assembly_accession.txt ${DEST}/cfDNA_pipeline/
cp sequencing_prep.tsv ${DEST}/cfDNA_pipeline/
cp reference_methylomes.tsv ${DEST}/cfDNA_pipeline/
cp Snakefile ${DEST}/cfDNA_pipeline/

find ${DEST}/cfDNA_pipeline -type d -exec chmod 755 {} \;
find ${DEST}/cfDNA_pipeline -type f -exec chmod 666 {} \;

find ${DEST}/cfDNA_pipeline/workflow/scripts/ -type d -exec chmod 755 {} \;
find ${DEST}/cfDNA_pipeline/workflow/scripts/ -type f -exec chmod 777 {}\;
