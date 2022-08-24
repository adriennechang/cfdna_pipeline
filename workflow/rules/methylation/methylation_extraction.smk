rule methylation_extraction:
	input: name_sorted = OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/{sample}.deduplicated.bam'
	output: CpG_bg = OUTPUT + '{project}/{sample}/methylation_extraction/bedGraph/{sample}.bedGraph.gz',
		CpG_bg_gz = temp(OUTPUT + '{project}/{sample}/methylation_extraction/{sample}.deduplicated.bedGraph.gz'),
		CpG_bismark = OUTPUT + '{project}/{sample}/methylation_extraction/bismark/{sample}.bismark.cov.gz',
		mbias = OUTPUT + '{project}/{sample}/methylation_extraction/mbias/{sample}.M-bias.txt',
		CHGOB = OUTPUT + '{project}/{sample}/methylation_extraction/CHG/{sample}_CHG_OB.txt.gz',
		CHGOT = OUTPUT + '{project}/{sample}/methylation_extraction/CHG/{sample}_CHG_OT.txt.gz',
		CHHOB = OUTPUT + '{project}/{sample}/methylation_extraction/CHH/{sample}_CHH_OB.txt.gz',
		CHHOT = OUTPUT + '{project}/{sample}/methylation_extraction/CHH/{sample}_CHH_OT.txt.gz',
		CpGOB = OUTPUT + '{project}/{sample}/methylation_extraction/CpG/{sample}_CpG_OB.txt.gz',
		CpGOT = OUTPUT + '{project}/{sample}/methylation_extraction/CpG/{sample}_CpG_OT.txt.gz',
		log = 'logs/{project}/{sample}.methylation_extraction.log'
	threads: 6
	params: outdir = OUTPUT + '{project}/{sample}/methylation_extraction/'
	shell:
		"""
		{METHEXT} --parallel {threads} \
			--bedGraph \
			-p \ 
			-o {params.outdir} \
			--gzip \
			{input.name_sorted};
		gunzip -c {params.outdir}{wildcards.sample}.deduplicated.bedGraph.gz | \
			sort-bed - | gzip > {output.CpG_bg};
		mv {params.outdir}{wildcards.sample}.deduplicated.bismark.cov.gz {output.CpG_bismark};
		mv {params.outdir}{wildcards.sample}.deduplicated_splitting_report.txt {output.log};
		mv {params.outdir}{wildcards.sample}.deduplicated.M-bias.txt {output.mbias};
		mv {params.outdir}CHG_OT_{wildcards.sample}.deduplicated.txt.gz {output.CHGOT};
		mv {params.outdir}CHG_OB_{wildcards.sample}.deduplicated.txt.gz {output.CHGOB};
		mv {params.outdir}CHH_OT_{wildcards.sample}.deduplicated.txt.gz {output.CHHOT};
		mv {params.outdir}CHH_OB_{wildcards.sample}.deduplicated.txt.gz {output.CHHOB};
		mv {params.outdir}CpG_OT_{wildcards.sample}.deduplicated.txt.gz {output.CpGOT};
		mv {params.outdir}CpG_OB_{wildcards.sample}.deduplicated.txt.gz {output.CpGOB};
		"""
