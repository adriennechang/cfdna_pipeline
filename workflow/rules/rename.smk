rule rename:
	input: i1 = DATA + '{project}/{sample}.R1.fastq.gz', i2 = DATA + '{project}/{sample}.R2.fastq.gz'
	output: o1 = DATA + '{project}/{sample}_R1.fastq.gz', o2 = DATA + '{project}/{sample}_R2.fastq.gz'
	shell:
		"""
		cp {input.i1} {output.o1};
		cp {input.i2} {output.o2};
		"""
rule touch:
	input: o1 = DATA + '{project}/{sample}_R1.fastq.gz', o2 = DATA + '{project}/{sample}_R2.fastq.gz'
	shell: "touch {input.o1}; touch {input.o2}"
