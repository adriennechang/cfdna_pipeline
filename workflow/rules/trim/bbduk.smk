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

def get_sample_info(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
        if str(wcrds.sample) == sample_id:
            return[sample_id, project_id, prep_type, path, genome_ver, seq_type]

rule bbduk:
	input: adapter_file = get_adapter_file,
		i1 = OUTPUT + '{project}/{sample}/{sample}_R1.fastq',
		i2 = OUTPUT + '{project}/{sample}/{sample}_R2.fastq'
	output: o1 = temp(OUTPUT + '{project}/{sample}/{sample}_R1_trim.fastq'),
		o2 = temp(OUTPUT + '{project}/{sample}/{sample}_R2_trim.fastq')
	threads: 6
	log: 'logs/{project}/{sample}/{sample}/bbduk.log'
	params:
		min_avg_phred='10',
		min_entropy='0.25',
		prep_type = get_sample_info,
		mem_mb= 1000
	shell:
		"""
		{BBDUK} in1={input.i1} in2={input.i2} \
			out1={output.o1} out2={output.o2} -Xmx1g \
			threads={threads} ref={input.adapter_file} \
			maq={params.min_avg_phred} entropy={params.min_entropy} \
			tbo tpe &>{log};
		"""
