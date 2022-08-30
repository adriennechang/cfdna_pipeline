def get_adapter_file(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                prep = prep_type
    if prep == "MEYER":
        adapt = MEYER
    if prep == "SRSLY":
        adapt = SRSLY
    return(adapt)
    

rule trimmomatic:
	input: i1 = OUTPUT + '{project}/{sample}/{sample}_R1.fastq',
		i2 = OUTPUT + '{project}/{sample}/{sample}_R2.fastq',
		adapter_file = get_adapter_file,
	output: o1 = temp(OUTPUT + '{project}/{sample}/{sample}_R1.paired.fastq'),
		o2 = temp(OUTPUT + '{project}/{sample}/{sample}_R1.unpaired.fastq'),
		o3 = temp(OUTPUT + '{project}/{sample}/{sample}_R2.paired.fastq'),
		o4 = temp(OUTPUT + '{project}/{sample}/{sample}_R2.unpaired.fastq'),
	threads: 4
	log: 'logs/{project}/{sample}.trimmomatic.log'
	shell:
		"""
		trimmomatic PE -phred33 \
			-threads {threads} \
			{input.i1} {input.i2} \
			{output.o1} {output.o2} {output.o3} {output.o4} \
			ILLUMINACLIP:{input.adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 &>{log};
		"""
