rule merge:
	input: i1 = OUTPUT + '{project}/{sample}/{sample}_R1.paired.fastq',
		i2 = OUTPUT + '{project}/{sample}/{sample}_R2.paired.fastq',
		i3 = OUTPUT + '{project}/{sample}/{sample}_R1.unpaired.fastq',
		i4 = OUTPUT + '{project}/{sample}/{sample}_R2.unpaired.fastq'
	output: o1 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}_merged.extendedFrags.fastq'),
		o2 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}_merged.notCombined_1.fastq'),
		o3 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}_merged.notCombined_2.fastq'),
		o4 = temp(OUTPUT + '{project}/{sample}/{sample}_unpaired.merged.fastq'),
		o5 = temp(OUTPUT + '{project}/{sample}/{sample}_R1.merged.fastq'),
		o6 = temp(OUTPUT + '{project}/{sample}/{sample}_R2.merged.fastq'),
	threads: 4
	log: 'logs/{project}/{sample}.merge.log'
	shell:
		"""
		(flash -m 10 {input.i1} {input.i2} -d {OUTPUT}/{wildcards.project}/{wildcards.sample} -o tmp_{wildcards.sample}_merged) &>{log};
		sed -e '/@/ s/\ .*$/-C/g' {output.o1} > {output.o4};
		cat {output.o2} > {output.o5};
		cat {output.o3} > {output.o6};
		"""
