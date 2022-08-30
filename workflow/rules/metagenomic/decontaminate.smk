def get_adapter_file(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver,seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                prep = prep_type
    if prep == "MEYER":
        adapt = MEYER
    if prep == "SRSLY":
        adapt = SRSLY
    if prep == "NEXTERA":
        adapt = NEXTERA
    return(adapt)


rule subset_FQ:
	input: i1 = OUTPUT + '{project}/{sample}/phix/{sample}_unmapped_R1.fastq',
                i2 = OUTPUT + '{project}/{sample}/phix/{sample}_unmapped_R2.fastq',
		raw1 = DATA + '{project}/{sample}_R1.fastq.gz',
		raw2 = DATA + '{project}/{sample}_R2.fastq.gz'
	output: temp1 = temp(OUTPUT + '{project}/{sample}/{sample}_R1.fq'),
		temp2 = temp(OUTPUT + '{project}/{sample}/{sample}_R2.fq'),
		o1 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_R1.fastq',
		o2 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_R2.fastq',
	shell:
		"""
		zcat {input.raw1} | sed 's/ 1.*$/-1/g' > {output.temp1};
		zcat {input.raw2} | sed 's/ 2.*$/-2/g' > {output.temp2};
		LC_ALL=C grep "@" {input.i1} | sed 's/^@//g' | seqtk subseq {output.temp1} - > {output.o1};
		LC_ALL=C grep "@" {input.i2} | sed 's/^@//g' | seqtk subseq {output.temp2} - > {output.o2};
		"""

rule adapt_qual:
	input: o1 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_R1.fastq',
                o2 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_R2.fastq',
		adapter = get_adapter_file
	output: r1 =  OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_R1.fastq',
		r2 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_R2.fastq'
	shell:
		"""
		{BBDUK} in1={input.o1} in2={input.o2} \
			out1={output.r1} out2={output.r2} -Xmx1g \
			ref={input.adapter} \
			threads=1 \
			maq=32;
		"""

rule dedupe:
	input: r1 =  OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_R1.fastq',
                r2 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_R2.fastq'
	output: r1 =  OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_dedup_R1.fastq',
                r2 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_dedup_R2.fastq'
	resources: mem_mb = 20000
	log: 'logs/{project}/{sample}.dedupe.log'
	shell:
		"""
		{DEDUP} -Xmx65000m in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} dedupe &> {log};
		"""

rule flash:
	input: r1 =  OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_dedup_R1.fastq',
                r2 = OUTPUT + '{project}/{sample}/{sample}_original_unmapped_trim_dedup_R2.fastq'
	output: combined = OUTPUT + '{project}/{sample}/{sample}_flash.extendedFrags.fastq',
		notcombined = OUTPUT + '{project}/{sample}/{sample}_flash.notCombined_1.fastq'
	resources: mem_mb=20000
	log: 'logs/{project}/{sample}.flash.log'
	shell:
		"""
		flash -q -M75 -O -o {wildcards.sample}_flash -d {OUTPUT}{wildcards.project}/{wildcards.sample} {input.r1} {input.r2} &> {log};
		"""

rule combine:
	input: combined = OUTPUT + '{project}/{sample}/{sample}_flash.extendedFrags.fastq',
                notcombined = OUTPUT + '{project}/{sample}/{sample}_flash.notCombined_1.fastq'
	output: o1 = OUTPUT + '{project}/{sample}/{sample}.flash.cat.fastq'
	shell:
		"""
		sed '1~4s/@.*/&-COMBINED/' {input.combined} > {output.o1};
		sed '1~4s/@.*/&-RISINGLEEND/' {input.notcombined} >> {output.o1};
		"""

rule decontaminate:
	input: cat =OUTPUT + '{project}/{sample}/{sample}.flash.cat.fastq',
		CT = CT_DECON,
		GA = GA_DECON
	output:
		r1_CT = temp(OUTPUT + '{project}/{sample}/decontaminate/{sample}.CT.fastq'),
		r1_pass1 = temp(OUTPUT + '{project}/{sample}/decontaminate/{sample}.pass1.fastq'),
		r1_pass2 = temp(OUTPUT + '{project}/{sample}/decontaminate/{sample}.pass2.fastq'),
		decon_r1 = temp(OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.fa'),
		nonhu_ct = temp(OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.CT.fa'),
	resources: mem_mb=80000
	shell:
		"""
		python {INSIL_CONV} -i {input.cat} -o {output.r1_CT} -r R1;
		{BBDUK} in={output.r1_CT} \
			out={output.r1_pass1} \
			-Xmx{resources.mem_mb}m \
			prealloc=t rcomp=f \
			ref={input.GA} k=50;
		{BBDUK} in={output.r1_pass1} \
			out={output.r1_pass2} \
			-Xmx{resources.mem_mb}m \
			prealloc=t rcomp=f \
			ref={input.CT} k=50;
		grep "@" {output.r1_pass2} | sed 's/^@//g' | seqtk subseq {input.cat} - | \
			fastq_to_fasta -Q33 -o {output.decon_r1};
		fastq_to_fasta -Q33 -i {output.r1_pass2} -o {output.nonhu_ct};
		"""
