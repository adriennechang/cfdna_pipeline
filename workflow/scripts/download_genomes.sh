#!/bin/bash

outdir=$1
outfile=$2
outall=$3

wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt -O $outdir"archaea.txt"
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O $outdir"bacteria.txt"
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt -O $outdir"fungi.txt"
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt -O $outdir"viral.txt"
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/assembly_summary.txt -O $outdir"human.txt"

cat $outdir"archaea.txt" $outdir"bacteria.txt" $outdir"fungi.txt" $outdir"viral.txt" $outdir"human.txt" > $outall


awk -F "\t" '{ if ($5== "representative genome" || $5 == "reference genome" ) print $0}' $outall > $outfile || true;

LC_ALL=C fgrep -f add_assembly_accession.txt $outall >> $outfile || true;
