def get_host_genome(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver = genome_ver
    if ver == "hg38":
        gen = HG38_CHROMO_SIZES
    if ver == "hg19":
        gen = HG19_CHROMO_SIZES
    return(gen)


rule aln_stats:
	input: raw = DATA + '{project}/{sample}_R1.fastq.gz',
		nhpair = OUTPUT + '{project}/{sample}/{sample}.nonhuman_paired.fastq',
		nhunpair = OUTPUT + '{project}/{sample}/{sample}.nonhuman_unpaired.fastq',
		nphix1 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R1.fastq.gz',
		nphunpair = OUTPUT + '{project}/{sample}/{sample}.nonhost_unpaired.fastq.gz',
		normalbam = OUTPUT + '{project}/{sample}/{sample}.host.pos.sorted.bam',
		i1 = OUTPUT + '{project}/{sample}/{sample}.host.pos.sorted.bam.bai',
		i2 = OUTPUT + '{project}/{sample}/{sample}.host.duprmv.sorted.bam.bai',
		rmvdupbam = OUTPUT + '{project}/{sample}/{sample}.host.duprmv.sorted.bam',
		genome = get_host_genome,
	output: stat = OUTPUT + '{project}/{sample}/{sample}_stats.align.tab',
		o1 = temp(OUTPUT + '{project}/{sample}/{sample}.nonhost.R1.fastq'),
		o2 = temp(OUTPUT + '{project}/{sample}/{sample}.nonhost_unpaired.fastq'),
	params: chr = 'chr21'
	threads: 4
	shell: 
		"""
		gunzip -c {input.nphix1} > {output.o1};
		gunzip -c {input.nphunpair} > {output.o2};
		bash {COUNT_FASTQ} {input.nhpair} nonhg19_pairs {wildcards.sample} > {output.stat};
		bash {COUNT_FASTQ} {input.nhunpair} nonhg19_orphans {wildcards.sample} >> {output.stat};
		bash {COUNT_FASTQ} {output.o1} nonphiX_pairs {wildcards.sample} >> {output.stat};
		bash {COUNT_FASTQ} {output.o2} nonphiX_orphans {wildcards.sample} >> {output.stat};
		bash {CALC_COV} {input.rmvdupbam} {params.chr} {wildcards.sample} {input.genome} >> {output.stat};
		bash {NON_DUP} {input.normalbam} {input.rmvdupbam} {wildcards.sample} >> {output.stat};
		totalreads=($( gunzip -c {input.raw} | LC_ALL=C grep "@" | wc -l ));
		echo "total_reads" {wildcards.sample} $totalreads >> {output.stat};
		"""
