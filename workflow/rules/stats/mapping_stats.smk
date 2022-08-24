rule mapping_stats:
	input: raw = DATA + '{project}/{sample}_R1.fastq.gz',
`		trimmed = OUTPUT + '{project}/{sample}/{sample}_R1_trim.fastq'),
		raw_bam = OUTPUT + '{project}/{sample}/aligned/raw_aligned/{sample}.bam',
		deduped_bam = OUTPUT + '{project}/{sample}/aligned/all_chr/{sample}_mapped_all_chr.bam',
		deduped_bai = OUTPUT + '{project}/{sample}/aligned/all_chr/{sample}_mapped_all_chr.bam.bai',
		unmapped_r1 = OUTPUT + '{project}/{sample}/aligned/unmapped/{sample}_unmapped_R1.fastq.gz',
		decon = OUTPUT


stats = OUTPUT + '{project}/{sample}/{sample}_mapping_stats.txt'

