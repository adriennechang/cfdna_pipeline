def get_reference_genome(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver = genome_ver
                seq = seq_type
    if seq == "bowtie":
        if ver == "hg38":
            gen = HG38STD
        if ver == "hg19":
            gen = HG19STD
    if seq == "bisulfite":
        if ver == "hg38":
            gen = HG38METH
        if ver == "hg19":
            gen = HG19METH
    return(gen)

def get_reference_genome_fasta(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver=genome_ver
                seq=seq_type
    if seq == "bowtie":
        if ver == "hg38":
            gen = HG38STD + 'hg38.fa'
        if ver == "hg19":
            gen = HG19STD + 'hg19.fa'
    if seq == "bisulfite":
        if ver == "hg38":
            gen = HG38METH + 'hg38.fa'
        if ver == "hg19":
            gen = HG19METH + 'hg19.fa'
    return(gen)


rule host_bismark:
	input: genome = get_reference_genome,
		i1 = OUTPUT + '{project}/{sample}/{sample}_R1_trim.fastq',
		i2 = OUTPUT + '{project}/{sample}/{sample}_R2_trim.fastq'
	output: bam = temp(OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/{sample}.bam'),
		o1 = temp(OUTPUT + '{project}/{sample}/bismark_aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz'),
		o2 = temp(OUTPUT + '{project}/{sample}/bismark_aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'),
	log: 'logs/{project}/{sample}.bismark.log'
	threads: 4
	params: outdir = OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/'
	shell:
		"""
		{BISMARK} --genome {input.genome} \
			--parallel {threads} \
			--quiet \
			--unmapped \
			-o {params.outdir} \
			-1 {input.i1} \
			-2 {input.i2};
		mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_pe.bam {output.bam};
		mv {params.outdir}{wildcards.sample}_R1_trim_bismark_bt2_PE_report.txt {log};
		mv {params.outdir}{wildcards.sample}_R1_trim.fastq_unmapped_reads_1.fq.gz {output.o1};
		mv {params.outdir}{wildcards.sample}_R2_trim.fastq_unmapped_reads_2.fq.gz {output.o2};
		"""

rule filter_bismark:
	input: bam = OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/{sample}.bam'
	output: dup = temp(OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/{sample}.deduplicated.bam'),
		mapped_all_chr = temp(OUTPUT + '{project}/{sample}/bismark_aligned/all_chr/{sample}_mapped_all_chr.bam'),
		mapped_all_chr_bai = temp(OUTPUT + '{project}/{sample}/bismark_aligned/all_chr/{sample}_mapped_all_chr.bam.bai'),
		dedup_file = temp(OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/{sample}.deduplication_report.txt')
	params:
		mapQ='15',
		outdir = OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/'
	log: 'logs/{project}/{sample}.deduplication.log'
	threads: 1
	shell:
		"""
		{RMDUPS} -o {wildcards.sample} -p --output_dir {params.outdir} --bam {input.bam} &> {log};
		samtools sort -@ {threads} {output.dup} | samtools view -h -q {params.mapQ} -o {output.mapped_all_chr};
		samtools index {output.mapped_all_chr}
		"""

rule bowtie2:
	input:
		genome = get_reference_genome_fasta,
		 i1 = OUTPUT + '{project}/{sample}/{sample}_R1_trim.fastq',
                i2 = OUTPUT + '{project}/{sample}/{sample}_R2_trim.fastq'
	output: bam = OUTPUT + '{project}/{sample}/bowtie_aligned/raw_aligned/{sample}.bam',
		o1 = OUTPUT + '{project}/{sample}/bowtie_aligned/unmapped/{sample}_unmapped_R1.fastq.gz',
		o2 = OUTPUT + '{project}/{sample}/bowtie_aligned/unmapped/{sample}_unmapped_R2.fastq.gz'
	log: 'logs/{project}/{sample}.bowtie.log'
	params: unmapped = OUTPUT + '{project}/{sample}/bowtie_aligned/unmapped/{sample}_unmapped_R%.fastq.gz'
	threads: 4
	shell:
		"""
		genome_prefix=$(echo {input.genome} | cut -f1 -d'.')
		bowtie2 -x $genome_prefix \
			-1 {input.i1} -2 {input.i2} \
			-p {threads} -q --score-min L,0,-0.2 --ignore-quals --no-mixed --no-discordant --dovetail --maxins 500 \
			--un-conc-gz {params.unmapped} | samtools view -f3 -bh - > {output.bam}
		"""

rule filter_bowtie2:
	input: bam= OUTPUT + '{project}/{sample}/bowtie_aligned/raw_aligned/{sample}.bam',
	output: mapped_all_chr = temp(OUTPUT + '{project}/{sample}/bowtie_aligned/all_chr/{sample}_mapped_all_chr.bam'),
		mapped_all_chr_bai = temp(OUTPUT + '{project}/{sample}/bowtie_aligned/all_chr/{sample}_mapped_all_chr.bam.bai')
	params:
		mapQ='15'
	threads: 1
	shell:
		"""
		samtools sort -n -@ {threads} {input.bam} | \
			samtools fixmate -m -@ {threads} - - | \
			samtools sort -@ {threads} | \
			samtools markdup -r -@ {threads} - - | samtools sort -@ {threads} -o {output.mapped_all_chr};
		samtools index {output.mapped_all_chr};
		"""


def aggregate_bam_files(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                seq = seq_type
    if seq == "bisulfite":
        outbam = OUTPUT + '{project}/{sample}/bismark_aligned/raw_aligned/{sample}.bam'
        mapped_all_chr = OUTPUT + '{project}/{sample}/bismark_aligned/all_chr/{sample}_mapped_all_chr.bam'
        mapped_all_chr_bai = OUTPUT + '{project}/{sample}/bismark_aligned/all_chr/{sample}_mapped_all_chr.bam.bai'
        unmapped_R1 = OUTPUT + '{project}/{sample}/bismark_aligned/unmapped/{sample}_pe_unmapped_R1.fastq.gz'
        unmapped_R2 = OUTPUT + '{project}/{sample}/bismark_aligned/unmapped/{sample}_pe_unmapped_R2.fastq.gz'
    if seq == "bowtie":
        outbam = OUTPUT + '{project}/{sample}/bowtie_aligned/raw_aligned/{sample}.bam'
        mapped_all_chr = OUTPUT + '{project}/{sample}/bowtie_aligned/all_chr/{sample}_mapped_all_chr.bam'
        mapped_all_chr_bai = OUTPUT + '{project}/{sample}/bowtie_aligned/all_chr/{sample}_mapped_all_chr.bam.bai'
        unmapped_R1 =  OUTPUT + '{project}/{sample}/bowtie_aligned/unmapped/{sample}_unmapped_R1.fastq.gz'
        unmapped_R2 =  OUTPUT + '{project}/{sample}/bowtie_aligned/unmapped/{sample}_unmapped_R2.fastq.gz'
    return[outbam, mapped_all_chr, mapped_all_chr_bai, unmapped_R1, unmapped_R2]

rule aggregate_bams:
	input: aggregate_bam_files
	output: bamout = OUTPUT + '{project}/{sample}/aligned/raw_aligned/{sample}.bam',
		mapped_all_chr = OUTPUT + '{project}/{sample}/aligned/all_chr/{sample}_mapped_all_chr.bam',
		mapped_all_chr_bai = OUTPUT + '{project}/{sample}/aligned/all_chr/{sample}_mapped_all_chr.bam.bai',
		unmapped_R1 = OUTPUT + '{project}/{sample}/aligned/unmapped/{sample}_unmapped_R1.fastq.gz',
		unmapped_R2 = OUTPUT + '{project}/{sample}/aligned/unmapped/{sample}_unmapped_R2.fastq.gz'
	shell:
		"""
		cp {input[0]} {output.bamout};
		cp {input[1]} {output.mapped_all_chr};
		cp {input[2]} {output.mapped_all_chr_bai};
		cp {input[3]} {output.unmapped_R1};
		cp {input[4]} {output.unmapped_R2};
		"""

rule phix_alignment:
	input: i1 = OUTPUT + '{project}/{sample}/aligned/unmapped/{sample}_unmapped_R1.fastq.gz',
		i2 = OUTPUT + '{project}/{sample}/aligned/unmapped/{sample}_unmapped_R2.fastq.gz',
	output: 
		bam = temp(OUTPUT + '{project}/{sample}/phix/{sample}.phix.bam'),
		o1 = OUTPUT + '{project}/{sample}/phix/{sample}_unmapped_R1.fastq',
		o2 = OUTPUT +  '{project}/{sample}/phix/{sample}_unmapped_R2.fastq'
	params:
		genome_prefix = PHIXSTD + '/phix',
		unmapped_pe = OUTPUT + '{project}/{sample}/phix/{sample}_unmapped_R%.fastq',
	shell:
		"""
		bowtie2 -x {params.genome_prefix} \
			-1 {input.i1} -2 {input.i2} \
			--local --very-sensitive-local \
			--un-conc {params.unmapped_pe} | \
			samtools view -bh - > {output.bam};
		"""	



