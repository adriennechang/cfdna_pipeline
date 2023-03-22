rule fastqc:
	input: i1 =DATA + '{project}/{sample}_R1.fastq.gz',
		i2 = DATA + '{project}/{sample}_R2.fastq.gz'
	output: o1 = temp(OUTPUT + '{project}/{sample}/{sample}_R1.fastq'),
            o2 = temp(OUTPUT + '{project}/{sample}/{sample}_R2.fastq'),
            o3 = OUTPUT + '{project}/{sample}/{sample}_R1_fastqc.html',
            o4 = OUTPUT + '{project}/{sample}/{sample}_R2_fastqc.html',
	threads: 2
	shell:
		"""
 		zcat {input.i1} | {ADD_STR} -s "-1" > {output.o1} ;
        	zcat {input.i2} | {ADD_STR} -s "-2" > {output.o2} ;
        	fastqc -t {threads} --extract {output.o1} ;
        	fastqc -t {threads} --extract {output.o2} ;
        	"""
