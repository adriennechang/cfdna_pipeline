#!/usr/bin/env bash

bismark=$1
methref=$2
deduplicate=$3
methylome=$4
out_bam=$5


mkdir -p references/reference_methylomes/sra_processing/
wd=references/reference_methylomes/sra_processing/$methylome
mkdir -p $wd

#if [[ $methylome == "lung_5" ]]; then
#  #prefetch SRR5263257 --max-size 100GB
#  fastq-dump --split-files /home/apc88/ncbi/public/sra/SRR5263257.sra --outdir $wd
#  $bismark --genome $methref -1 $wd/SRR5263257_1.fastq -2 $wd/SRR5263257_2.fastq -o $wd --parallel 20
#  $deduplicate $wd/SRR5263257_1_bismark_bt2_pe.bam --output_dir $wd
#  mv $wd/SRR5263257_1_bismark_bt2_pe.deduplicated.bam $out_bam
#elif [[ $methylome == "lung_6" ]]; then
#  #prefetch SRR5263259 --max-size 100GB
#  #fastq-dump --split-files /home/apc88/ncbi/public/sra/SRR5263259.sra --outdir $wd
#  #$bismark --genome $methref -1 $wd/SRR5263259_1.fastq -2 $wd/SRR5263259_2.fastq -o $wd --parallel 20
#  $deduplicate $wd/SRR5263259_1_bismark_bt2_pe.bam --output_dir $wd
#  mv $wd/SRR5263259_1_bismark_bt2_pe.deduplicated.bam $out_bam
#
#fi
