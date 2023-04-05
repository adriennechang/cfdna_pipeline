if OLD_MET_REF == "FALSE":
    if MAKE_NEWDB == "FALSE":
        use_db = NCBI_22
        use_gi = GI_TAXID_22
        use_length = GI_LENGTH_22
        use_gdt =  GDT_22
        use_lut = LUT_22
        filter_strand = FILTER_STRAND
        grammy_pre = GRAMMY_PRE_ACC
        human_gis = HUMAN_22
    else:
        use_db = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std.fna'
        use_gi = 'references/' + NEWDB_NAME +  '/' + NEWDB_NAME + '.acc.taxids'
        use_length = 'references/' + NEWDB_NAME + '/taxids_lengths.txt'
        use_gdt = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std'
        use_lut = 'references/' + NEWDB_NAME + '/LUT/taxids_names_lengths_tax.tab'
        grammy_pre = GRAMMY_PRE_ACC
        filter_strand = FILTER_STRAND
        human_gis = 'references/' + NEWDB_NAME + '/' + NEWDB_NAME + '.human.gis'

if OLD_MET_REF == "TRUE":
    use_db = NCBI_06
    use_gi = GI_TAXID_06
    use_length = GI_LENGTH_06
    use_gdt = GDT_06
    use_lut = LUT_06
    filter_strand = FILTER_STRAND_GI
    grammy_pre = GRAMMY_PRE_GI


def get_reference_genome(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver,seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver = genome_ver
    if ver == "hg19":
        gen = HG19STD + 'hg19.fa'
    if ver == "hg38":
        gen = HG38STD + 'hg38.fa'
    return(gen)

rule decontaminate_stds:
	input: i1 =  OUTPUT + '{project}/{sample}/{sample}.flash.cat.fastq',
		ref = get_reference_genome
	output: pass0 = OUTPUT + '{project}/{sample}/decontaminate/{sample}.pass.fastq',
		pass1 = OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.stds.fa'
	resources: mem_mb=80000
	shell:
		"""
		{BBDUK} in={input.i1} out={output.pass0} -Xmx{resources.mem_mb}m prealloc=t ref={input.ref} k=50;
		grep "@" {output.pass0} | sed 's/^@//g' | seqtk subseq {input.i1} - | fastq_to_fasta -Q33 -o {output.pass1};
		"""

rule hs_blastn_stds:
	input: i1 = OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.stds.fa',
	output: o1 = OUTPUT + '{project}/{sample}/blast/{sample}.outfmt6'
	threads: 8
	resources: mem_mb = 80000
	shell:
		"""
		{HSBLASTN} align -query {input.i1} \
			-db {use_db} \
			-evalue 0.0001 \
			-perc_identity 95 \
			-num_threads {threads} \
			-outfmt 6 > {output.o1};
		"""

rule filter_blast_stds:
	input: i1 = OUTPUT + '{project}/{sample}/blast/{sample}.outfmt6',
		clean = OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.stds.fa'
	output: hblast = temp(OUTPUT + '{project}/{sample}/blast/{sample}.human'),
		blast_filt = OUTPUT + '{project}/{sample}/blast/{sample}.sorted.outfmt6'
	run:
		import os
		import dnaio
		import os.path
		
		if os.path.exists(human_gis):
		    cmd = f"grep -f {human_gis} {input.i1} | cut -f1 | sort -u > {output.hblast}"
		    os.system(cmd)
		    human_readIDs = [line.rstrip('\n') for line in open(output.hblast)]
		else:
		    cmd = f"echo 'empty' > {output.hblast}"
		    os.system(cmd)
		    human_readIDs = []

		with open(input.i1) as f, open(output.blast_filt, 'w') as w:
		    for hit in f:
		        read_id = hit.strip().split('\t')[0]
		        mol_len = float(hit.strip().split('\t')[3])
		        mapped_len = float(hit.strip().split('\t')[12])
		        ratio = mol_len/mapped_len
		        if (mol_len >=50) and (ratio > 0.9):
		            if human_readIDs != [] and read_id not in human_readIDs:
		                w.write(hit)
		            if human_readIDs == []:
		                w.write(hit)

rule get_relevant_stds:
	input: blast_filt = OUTPUT + '{project}/{sample}/blast/{sample}.sorted.outfmt6',
		clean = OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.stds.fa'
	output: cleantblat = OUTPUT + '{project}/{sample}/blast/{sample}.tblat.1'
	shell:
		"""
		grep "Plus/Plus" {input.blast_filt} > {output.cleantblat};
		grep "Plus/Minus" {input.blast_filt} >> {output.cleantblat};
		"""

rule grammy_clean_stds:
	input: clean = OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.stds.fa',
		cleantblat = OUTPUT + '{project}/{sample}/blast/{sample}.tblat.1'
	output: fa_gz = OUTPUT + '{project}/{sample}/blast/{sample}.fa.gz',
		fasta_gz = temp(OUTPUT + '{project}/{sample}/blast/{sample}.fasta.gz'),
		rdt = temp(OUTPUT + '{project}/{sample}/blast/{sample}.rdt'),
		mtx = temp(OUTPUT + '{project}/{sample}/blast/{sample}.mtx'),
		lld = temp(OUTPUT + '{project}/{sample}/blast/{sample}.lld'),
		btp = temp(OUTPUT + '{project}/{sample}/blast/{sample}.btp'),
		est = temp(OUTPUT + '{project}/{sample}/blast/{sample}.est'),
		gra = temp(OUTPUT + '{project}/{sample}/blast/{sample}.gra'),
		avl = temp(OUTPUT + '{project}/{sample}/blast/{sample}.avl')
	resources: mem_mb =1
	shell:
		"""
		if [ $(wc -l {input.cleantblat} | cut -d' ' -f1) -gt 1 ]
		then
			cut -f1 {input.cleantblat} | sort -u | seqtk subseq {input.clean} - | gzip -1 > {output.fa_gz};

			if [[ {OLD_MET_REF} == "FALSE" ]]
			then
				docker1 run -v `pwd`/{OUTPUT}{wildcards.project}/{wildcards.sample}/blast:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy \
					bash -c "cd /data && \
						grammy_rdt -t illumina . . && \
						grammy_pre_acc -q "40,75,-5" {wildcards.sample} /{use_gdt} && \
						grammy_em -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx && \
						grammy_post {wildcards.sample}.est /{use_gdt} {wildcards.sample}.btp"
			else
				docker1 run -v `pwd`/{OUTPUT}/{wildcards.project}/{wildcards.sample}/blast:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy \
					bash -c "cd /data && \
						grammy_rdt -t illumina . . && \
						grammy_pre -q "40,75,-5" {wildcards.sample} /{use_gdt} && \
						grammy_em -c L -b 5 -t .00001 -n 100 {wildcards.sample.mtx && \
						grammy_post {wildcards.sample}.est /{use_gdt} {wildcards.sample}.bt["
			fi

			cp {output.fa_gz} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.fa_gz};
			cp {output.fasta_gz} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.fasta_gz};
			cp {output.rdt} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.rdt};
			cp {output.mtx} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.mtx};
			cp {output.lld} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.lld};
			cp {output.btp} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.btp};
			cp {output.est} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.est};
			cp {output.gra} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.gra};
			cp {output.avl} {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.avl};

			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.fa_gz} {output.fa_gz};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.fasta_gz} {output.fasta_gz};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.rdt} {output.rdt};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.mtx} {output.mtx};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.lld} {output.lld};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.btp} {output.btp};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.est} {output.est};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.gra} {output.gra};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/blast/tmp_{output.avl} {output.avl};



		else
			touch {output.fa_gz} {output.fasta_gz} {output.rdt} {output.mtx} {output.lld} {output.btp} {output.est} {output.gra} {output.avl};
		fi
		"""


rule annotate_grammy_stds:
	input: gra = OUTPUT + '{project}/{sample}/blast/{sample}.gra',
		cleantblat = OUTPUT + '{project}/{sample}/blast/{sample}.tblat.1',
		stats = OUTPUT + '{project}/{sample}/{sample}_mapping_stats.txt'
	output: tab = temp(OUTPUT + '{project}/{sample}/blast/{sample}.tab'),
		anno = OUTPUT + '{project}/{sample}/blast/{sample}.grammy.tab'
	shell:
		"""
		if [ $(cut -f1 {input.cleantblat} | sort -u | wc -l) -gt 2 ]
		then
			Rscript {FILT_GRAM_BS} {input.gra} {output.tab} {wildcards.sample};
			Rscript {ANN_GRAM_BS} {output.tab} {input.cleantblat} {input.stats} {use_lut} {output.anno};
		else
			touch {output.tab}
			head -n1 {use_lut} > {output.anno}
		fi
		"""
