def get_reference_genome(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver=genome_ver
    if ver == "hg19":
        map = 2948611470
    if ver == "hg38":
        map = 2937639396
    return(map)

rule mapping_stats:
	input: raw = DATA + '{project}/{sample}_R1.fastq.gz',
		trimmed = OUTPUT + '{project}/{sample}/{sample}_R1_trim.fastq',
		raw_bam = OUTPUT + '{project}/{sample}/aligned/raw_aligned/{sample}.bam',
		deduped_bam = OUTPUT + '{project}/{sample}/aligned/all_chr/{sample}_mapped_all_chr.bam',
		deduped_bai = OUTPUT + '{project}/{sample}/aligned/all_chr/{sample}_mapped_all_chr.bam.bai',
	output: stats = OUTPUT + '{project}/{sample}/{sample}_mapping_stats.txt'
	params: mappable_chr21 = '40088623',
		map = get_reference_genome
	shell:
		"""
		original_reads=$(pigz -dc {input.raw} | wc -l)
		original_reads=$((original_reads / 4 ))
		trimmed_reads=$(wc -l {input.trimmed} | cut -d' ' -f1)
		trimmed_reads=$((trimmed_reads/4))
		mapped_reads=$(samtools view -c {input.raw_bam})
		mapped_reads=$((mapped_reads/2))
		deduped_reads=$(samtools view -c {input.deduped_bam})
		deduped_reads=$((deduped_reads/2))
		mapping_eff=$(echo "scale=2;$mapped_reads/$trimmed_reads" | bc)
		deduped_frac=$(echo "scale=2;$deduped_reads/$mapped_reads" | bc)
		
		depth_var=$(samtools depth {input.deduped_bam} | \
			awk '{{ sum+=$3 }} END {{ print sum/{params.map}" "NR/{params.map} }}')
		depth=$(echo $depth_var | cut -f1 -d' ')
		bp_frac=$(echo $depth_var | cut -f2 -d' ')
		echo -e "SAMPLE\tNUM_READS\tREADS_AFTER_TRIM\tALIGNED\tMAPPING_EFFICIENCY\tDEDUPED_READS\tDEDUP_FRAC\tDEPTH\tFRACTION_PB_BP" > {output.stats};
        	echo -e "{wildcards.sample}\t$original_reads\t$trimmed_reads\t$mapped_reads\t$mapping_eff\t$deduped_reads\t$deduped_frac\t$depth\t$bp_frac" >> {output.stats};
		"""
