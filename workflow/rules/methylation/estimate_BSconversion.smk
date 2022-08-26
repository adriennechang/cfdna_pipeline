def get_reference_genome_fasta(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver = genome_ver
                seq = seq_type
    if seq == "standard":
        if ver == "hg38":
            gen = HG38STD + 'hg38.fa'
        if ver == "hg19":
            gen = HG19STD + 'hg19.fa'
    if seq == "bisulfite":
        if ver == "hg38":
            gen = HG38METH + 'hg38.fa'
        if ver == "hg19":
            gen = HG19METH + 'hg19.fa'
    return(gen)


rule estimate_BSconversion:
	input: mapped = OUTPUT + '{project}/{sample}/aligned/all_chr/{sample}_mapped_all_chr.bam',
		genome = get_reference_genome_fasta
	output:
		mr = temp(OUTPUT + '{project}/{sample}/{sample}.mr'),
		bsrate = OUTPUT + '{project}/{sample}/{sample}.bsrate.txt',
		sample_rate = OUTPUT + '{project}/{sample}/{sample}.bsconversion'
	shell:
		"""
		format_reads -o {output.mr} -L 500 {input.mapped};
		bsrate -c {input.genome} -o {output.bsrate} {output.mr};
		X=$(head -1 {output.bsrate})
		echo -e "{wildcards.sample}\t$X" > {output.sample_rate}
		"""
