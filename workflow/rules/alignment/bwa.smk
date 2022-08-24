def get_reference_genome(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver = genome_ver
    if ver == "hg38":
        gen = HG38
    if ver == "hg19":
        gen = HG19
    return(gen)

rule host_align_bwa:
	input: i1 = OUTPUT + '{project}/{sample}/{sample}_R1.merged.fastq',
		i2 = OUTPUT + '{project}/{sample}/{sample}_R2.merged.fastq',
		i3 = OUTPUT + '{project}/{sample}/{sample}_unpaired.merged.fastq',
		genome = get_reference_genome,
	output: prepairfq = temp(OUTPUT + '{project}/{sample}/{sample}.host.merged.inter.fastq'),
		preorphanfq = temp(OUTPUT + '{project}/{sample}/{sample}.host.merged.empty.fastq'),
		prepairfq2 = temp(OUTPUT + '{project}/{sample}/{sample}.host.merged.inter2.fastq'),
		o1 = OUTPUT + '{project}/{sample}/{sample}.host.bam',
		o2h = temp(OUTPUT + '{project}/{sample}/{sample}.host.name.sorted.bam'),
		o3h = temp(OUTPUT + '{project}/{sample}/{sample}.host.fixmate.bam'),
		o4h = temp(OUTPUT + '{project}/{sample}/{sample}.host.pos.sorted.bam'),
		o5 = temp(OUTPUT + '{project}/{sample}/{sample}.host.pos.sorted.bam.bai'),
		o6h = OUTPUT+'{project}/{sample}/{sample}.host.duprmv.sorted.bam',
		o6 = OUTPUT + '{project}/{sample}/{sample}.host.duprmv.sorted.bam.bai',
		o3 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.nonhuman_paired.bam'),
		o4 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.nonhuman_unpaired.bam'),
		fq11 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.nonhuman_paired.R1.fastq'),
		fq12 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.nonhuman_paired.R2.fastq'),
		fq2 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.nonhuman_unpaired.fastq'),
		tmpfqtemp1 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.host.tmp_1.t.fastq'),
		tmpfqtemp2 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.host.tmp_2.t.fastq'),
		tmpfq1 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.host.tmp_1.fastq'),
		tmpfq2 = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.host.tmp_2.fastq'),
		tmpfqC = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.host.tmp_C.fastq'),
		pairfq = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.host.tmp_pair.fastq'),
		orphanfq = temp(OUTPUT + '{project}/{sample}/tmp_{sample}.host.tmp_orphan.fastq'),
		unpfq = temp(OUTPUT + '{project}/{sample}/{sample}.nonhuman_unpaired.fastq'),
		finalfq1 = temp(OUTPUT + '{project}/{sample}/{sample}.nonhuman.R1.fastq'),
		finalfq2 = temp(OUTPUT + '{project}/{sample}/{sample}.nonhuman.R2.fastq'),
		finalpair = temp(OUTPUT + '{project}/{sample}/{sample}.nonhuman_paired.fastq'),
	threads: 4
	log: 'logs/{project}/{sample}.host.log'
	shell:
		"""
		sed -i 's/\ 1.*$/-1/g' {input.i1};
		sed -i 's/\ 2.*$/-2/g' {input.i2};
		sed -i 's/\-[1-2]/\-C/g' {input.i3};
		python2.7 {INTER} {input.i1} {input.i2} {output.prepairfq} {output.preorphanfq};
		sed 's/-1//g; s/-2//g' {output.prepairfq} | cat - {input.i3} > {output.prepairfq2};
		(bwa mem -p -t {threads} {input.genome} {output.prepairfq2} | \
			samtools view -F 1024,2048 -bS - > {output.o1}) &>{log};
		samtools sort -n -@ {threads} {output.o1} -o {output.o2h};
		samtools fixmate -m {output.o2h} {output.o3h};
		samtools sort -@ {threads} {output.o3h} -o {output.o4h};
		samtools markdup -r {output.o4h} {output.o6h};
		samtools index {output.o6h} > {output.o6};
		samtools index {output.o4h} > {output.o5};
		samtools view -b -f 0x000D {output.o1} -o {output.o3};
		samtools view -b -F 0x0001 -f 0x0004 {output.o1} -o {output.o4};
		bedtools bamtofastq -i {output.o3} -fq {output.fq11} -fq2 {output.fq12};
		bedtools bamtofastq -i {output.o4} -fq {output.fq2};
		sed -e '/@/ s/\/1/\-1/' {output.fq11} > {output.tmpfqtemp1};
		sed -e '/@/ s/\/2/\-2/' {output.fq12} > {output.tmpfqtemp2};
		cat {output.tmpfqtemp1} {output.tmpfqtemp2} {output.fq2} | \
			{F2S} | \
			grep -P "\-1\t" | \
			{S2F} > {output.tmpfq1};
		cat {output.tmpfqtemp1} {output.tmpfqtemp2} {output.fq2} | \
                        {F2S} | \
			grep -P "\-2\t" | \
			{S2F} > {output.tmpfq2};
		cat {output.tmpfqtemp1} {output.tmpfqtemp2} {output.fq2} | \
                        {F2S} | \
			grep -P "\-C\t" | \
			{S2F} > {output.tmpfqC};
		python2.7 {INTER} {output.tmpfq1} {output.tmpfq2} {output.pairfq} {output.orphanfq};
		cat {output.orphanfq} {output.tmpfqC} > {output.unpfq};
		cat {output.pairfq} | {F2S} | \
			grep -P "\-1\t" | {S2F} > {output.finalfq1};
		cat {output.pairfq} | {F2S} | \
			grep -P "\-2\t" | {S2F} > {output.finalfq2};
		cp {output.pairfq} {output.finalpair};
		"""

rule phix_align_bwa:
	input: i1 = OUTPUT+'{project}/{sample}/{sample}.nonhuman.R1.fastq',
		i2 = OUTPUT +'{project}/{sample}/{sample}.nonhuman.R2.fastq',
		i3 = OUTPUT + '{project}/{sample}/{sample}.nonhuman_unpaired.fastq',
	output: prepairfq = temp(OUTPUT +'{project}/{sample}/{sample}.merged.inter.phix.fastq'),
		preorphanfq = temp(OUTPUT+'{project}/{sample}/{sample}.merged.empty.phix.fastq'),
		prepairfq2 = temp(OUTPUT+'{project}/{sample}/{sample}.merged.inter2.phix.fastq'),
		o1 = temp(OUTPUT+'{project}/{sample}/{sample}.nonhuman_phiX.bam'),
		o3 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.nonhuman_nonphiX_paired.bam'),
		o4 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.nonhuman_nonphiX_unpaired.bam'),
		fq11 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.nonhuman_nonphiX_paired.R1.fastq'),
		fq12 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.nonhuman_nonphiX_paired.R2.fastq'),
		fq2 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.nonhuman_nonphpiX_unpaired.fastq'),
		tmpfqtemp1 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.tmp_1.t.phix.fastq'),
		tmpfqtemp2 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.tmp_2.t.phix.fastq'),
		tmpfq1 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.tmp_1.phix.fastq'),
		tmpfq2 = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.tmp_2.phix.fastq'),
		tmpfqC = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.tmp_C.phix.fastq'),
		pairfq = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.tmp_pair.phix.fastq'),
		orphanfq = temp(OUTPUT+'{project}/{sample}/tmp_{sample}.tmp_orphan.phix.fastq'),
		nh1 = OUTPUT+'{project}/{sample}/{sample}.nonhost.R1.fastq.gz',
		nh2 = OUTPUT+'{project}/{sample}/{sample}.nonhost.R2.fastq.gz',
		nhu = OUTPUT+'{project}/{sample}/{sample}.nonhost_unpaired.fastq.gz',
	threads: 4
	shell:
		"""
		python2.7 {INTER} {input.i1} {input.i2} {output.prepairfq} {output.preorphanfq};
		sed 's/-1//g; s/-2//g' {output.prepairfq} | cat - {input.i3} > {output.prepairfq2};
		bwa mem -p -t {threads} {PHIX} {output.prepairfq2} | \
			samtools view -F 1024,2048 -bS - > {output.o1};
		samtools view -b -f 0x000D {output.o1} -o {output.o3};
		samtools view -b -F 0x0001 -f 0x0004 {output.o1} -o {output.o4};
		bedtools bamtofastq -i {output.o3} -fq {output.fq11} -fq2 {output.fq12};
		bedtools bamtofastq -i {output.o4} -fq {output.fq2};
		sed -e '/@/ s/\/1/\-1/' {output.fq11} > {output.tmpfqtemp1};
		sed -e '/@/ s/\/2/\-2/' {output.fq12} > {output.tmpfqtemp2};
		cat {output.tmpfqtemp1} {output.tmpfqtemp2} {output.fq2} | \
			{F2S} | \
			grep -P "\-1\t" | {S2F} > {output.tmpfq1};
		cat {output.tmpfqtemp1} {output.tmpfqtemp2} {output.fq2} | \
                        {F2S} | \
			grep -P "\-2\t" | {S2F} > {output.tmpfq2};
		cat {output.tmpfqtemp1} {output.tmpfqtemp2} {output.fq2} | \
                        {F2S} | \
			grep -P "\-C\t" | {S2F} > {output.tmpfqC};
		python2.7 {INTER} {output.tmpfq1} {output.tmpfq2} {output.pairfq} {output.orphanfq};
		cat {output.orphanfq} {output.tmpfqC} | pigz > {output.nhu};
		cat {output.pairfq} | {F2S} | \
			grep -P "\-1\t" | {S2F} | \
			pigz > {output.nh1};
		cat {output.pairfq} | {F2S} | \
			grep -P "\-2\t" | {S2F} | \
			pigz > {output.nh2};
		"""
