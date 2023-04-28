if OLD_MET_REF == "FALSE":
    if MAKE_NEWDB == "FALSE":
        use_db = NCBI_22
        use_gi = GI_TAXID_22
        use_length = GI_LENGTH_22
        use_gdt =  GDT_22
        use_lut = LUT_22
        filter_strand = FILTER_STRAND
	grammy_pre = GRAMMY_PRE_ACC
    else:
        use_db = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std.fna'
        use_gi = 'references/' + NEWDB_NAME +  '/' + NEWDB_NAME + '.acc.taxids'
        use_length = 'references/' + NEWDB_NAME + '/taxids_lengths.txt'
        use_gdt = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std'
        use_lut = 'references/' + NEWDB_NAME + '/LUT/taxids_names_lengths_tax.tab'
        grammy_pre = GRAMMY_PRE_ACC
        filter_strand = FILTER_STRAND

if OLD_MET_REF == "TRUE":
    use_db = NCBI_06
    use_gi = GI_TAXID_06
    use_length = GI_LENGTH_06
    use_gdt = GDT_06
    use_lut = LUT_06
    filter_strand = FILTER_STRAND_GI
    grammy_pre = GRAMMY_PRE_GI


rule hsblastn:
	input:
		i1 = OUTPUT+ '{project}/{sample}/{sample}.nonhost.R1.fastq.gz',
		i2 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R2.fastq.gz',
		i3 = OUTPUT +'{project}/{sample}/{sample}.nonhost_unpaired.fastq.gz',
		db = use_db,
		gi_to_taxid = use_gi
	output:
		o1 = temp(OUTPUT+'{project}/{sample}/{sample}.nonhost.R1.fasta'),
		o2 = temp(OUTPUT+'{project}/{sample}/{sample}.nonhost.R2.fasta'),
		o3 = temp(OUTPUT+'{project}/{sample}/{sample}.nonhost_unpaired.fasta'),
		oall = temp(OUTPUT + '{project}/{sample}/{sample}.nonhost.fasta'),
		blast1 = temp(OUTPUT+'{project}/{sample}/{sample}.nonhost.outfmt6'),
		human_like1 = temp(OUTPUT + '{project}/{sample}/{sample}.humanlike'),
		rejected1 = temp(OUTPUT + '{project}/{sample}/{sample}.rejected'),
		tblat0 = temp(OUTPUT + '{project}/{sample}/{sample}.tblat.0'),
		tblat1 = OUTPUT + '{project}/{sample}/{sample}.tblat.1',
		orphans = temp(OUTPUT + '{project}/{sample}/{sample}_orphans.txt'),
	threads: 12
	resources:
		mem_mb=80000
	shell:
		"""
		zcat {input.i1} | fastq_to_fasta -Q33 -i - -o {output.o1};
		zcat {input.i2} | fastq_to_fasta -Q33 -i - -o {output.o2};
		zcat {input.i3} | fastq_to_fasta -Q33 -i - -o {output.o3};
		cat {output.o1} {output.o2} {output.o3} > {output.oall};
		{HSBLASTN} align -query {output.oall} \
                    	-db {input.db} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 | \
		python2 {filter_strand} --filename_out {output.blast1} \
			--acc_to_tax {input.gi_to_taxid} \
			--conversion std \
			--human_like {output.human_like1} \
			--rejected_hits {output.rejected1};
		cat {output.blast1} | \
			awk '{{ print $0, $4/$13 }} ' | \
			awk '{{ if ($NF > 0.90) print $0 }}' | \
			awk '{{ if ($3 > 90) print $0 }}' | \
			awk '{{ $13=""; print $0 }}' > {output.tblat0};
		zcat {input.i1} {input.i2} {input.i3} | {F2S} | \
			cut -f1 | grep -v -e 'C$$' > {output.orphans};

		python2.7 {RMDBL} {output.tblat0} {output.orphans} > {output.tblat1};
		sed -i -e 's/ /\t/g' {output.tblat1};
		"""


rule grammy:
	input:
		nonhumanfa_r1 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R1.fastq.gz',
		nonhumanfa_r2 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R2.fastq.gz',
		nonhuman_fa = OUTPUT + '{project}/{sample}/{sample}.nonhost_unpaired.fastq.gz',
		tblat1 = OUTPUT + '{project}/{sample}/{sample}.tblat.1'
	output:
		nonhumanfa_gz = temp(OUTPUT + '{project}/{sample}/{sample}.fa.gz'),
		nonhumanfasta_gz = temp(OUTPUT + '{project}/{sample}/{sample}.fasta.gz'),
		rdt = temp(OUTPUT + '{project}/{sample}/{sample}.rdt'),
		mtx = temp(OUTPUT + '{project}/{sample}/{sample}.mtx'),
		lld = temp(OUTPUT + '{project}/{sample}/{sample}.lld'),
		btp = temp(OUTPUT + '{project}/{sample}/{sample}.btp'),
		est = temp(OUTPUT + '{project}/{sample}/{sample}.est'),
		gra = temp(OUTPUT + '{project}/{sample}/{sample}.gra'),
		avl = temp(OUTPUT + '{project}/{sample}/{sample}.avl'),
		tab = temp(OUTPUT + '{project}/{sample}/{sample}.tab'),
	resources: mem_mb=1
	shell:
		"""
		zcat {input.nonhumanfa_r1} {input.nonhumanfa_r2} {input.nonhuman_fa} | \
			fastq_to_fasta -Q33 -i -  | sed 's/_1.*/-1/g' |  sed 's/_2.*/-2/g' | gzip -1 > {output.nonhumanfa_gz}

		if [[ {OLD_MET_REF} == "FALSE" ]]
                then
                        docker1 run -v `pwd`/results/{wildcards.project}/{wildcards.sample}:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy  \
                                bash -c "cd /data  &&
                                        grammy_rdt -t illumina . . && \
                                        grammy_pre_acc -q "40,75,-5" {wildcards.sample} /{use_gdt} && \
                                        grammy_em -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx && \
                                        grammy_post {wildcards.sample}.est /{use_gdt} {wildcards.sample}.btp"
		else
			docker1 run -v `pwd`/results/{wildcards.project}/{wildcards.sample}:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy  \
                                bash -c "cd /data  &&
                                        grammy_rdt -t illumina . . && \
                                        grammy_pre -q "40,75,-5" {wildcards.sample} /{use_gdt} && \
                                        grammy_em -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx && \
                                        grammy_post {wildcards.sample}.est /{use_gdt} {wildcards.sample}.btp" 
               fi
		mkdir -p results/{wildcards.project}/{wildcards.sample}/tmp_results/{wildcards.project}/{wildcards.sample};

		cp {output.nonhumanfa_gz} results/{wildcards.project}/{wildcards.sample}/tmp_{output.nonhumanfa_gz};
		cp {output.nonhumanfasta_gz} results/{wildcards.project}/{wildcards.sample}/tmp_{output.nonhumanfasta_gz};
		cp {output.rdt} results/{wildcards.project}/{wildcards.sample}/tmp_{output.rdt};
		cp {output.mtx} results/{wildcards.project}/{wildcards.sample}/tmp_{output.mtx};
		cp {output.lld} results/{wildcards.project}/{wildcards.sample}/tmp_{output.lld};
		cp {output.btp} results/{wildcards.project}/{wildcards.sample}/tmp_{output.btp};
		cp {output.est} results/{wildcards.project}/{wildcards.sample}/tmp_{output.est};
		cp {output.gra} results/{wildcards.project}/{wildcards.sample}/tmp_{output.gra};
		cp {output.avl} results/{wildcards.project}/{wildcards.sample}/tmp_{output.avl};

		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.nonhumanfa_gz} {output.nonhumanfa_gz};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.nonhumanfasta_gz} {output.nonhumanfasta_gz};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.rdt} {output.rdt};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.mtx} {output.mtx};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.lld} {output.lld};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.btp} {output.btp};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.est} {output.est};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.gra} {output.gra};
		mv -f results/{wildcards.project}/{wildcards.sample}/tmp_{output.avl} {output.avl};

		rm -r results/{wildcards.project}/{wildcards.sample}/tmp_results/

		cat {output.gra} | {TRANSPOSE} | \
                        {FILTER} -c 1 -mins 0 | \
                        {ADD_COL} -b -s "{wildcards.sample}" > {output.tab};

		"""

rule grammy_annotate:
	input: tab = OUTPUT + '{project}/{sample}/{sample}.tab',
                stat = OUTPUT + '{project}/{sample}/{sample}_stats.align.tab'
	output: o1 = OUTPUT + '{project}/{sample}/{sample}.grammy.tab'
	threads: 4
	shell:
		"""
		if [[ {OLD_MET_REF} == "FALSE" ]]
		then
			Rscript {ANN_GRAMMY_ACC} {OUTPUT}{wildcards.project}/{wildcards.sample}/ {wildcards.sample} {use_lut} {OUTPUT}{wildcards.project}/{wildcards.sample}/ {output.o1}
		else
			Rscript {ANN_GRAMMY_GI} {OUTPUT}{wildcards.project}/{wildcards.sample}/ {wildcards.sample} {use_lut}
		fi
		"""
