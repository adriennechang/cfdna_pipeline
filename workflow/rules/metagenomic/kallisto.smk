rule kallisto:
	input: nonhumanfa_r1 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R1.fastq.gz',
		nonhumanfa_r2 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R2.fastq.gz',
		nonhuman_fa = OUTPUT + '{project}/{sample}/{sample}.nonhost_unpaired.fastq.gz',
	output: tab = OUTPUT + '{project}/{sample}/abundance.tsv',
		bam = OUTPUT + '{project}/{sample}/{sample}.bam'
	threads: 80
	shell: 
		"""
		/programs/kallisto/kallisto quant -t 1 --pseudobam --single -i {KALIDX} -l 50 -s 20 -o results/{wildcards.project}/{wildcards.sample}/ {input.nonhumanfa_r1} {input.nonhumanfa_r2} {input.nonhuman_fa} > {output.bam};
		"""
