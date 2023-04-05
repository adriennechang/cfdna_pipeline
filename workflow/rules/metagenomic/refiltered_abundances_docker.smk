if OLD_MET_REF == "FALSE":
    if MAKE_NEWDB == "FALSE":
        use_db = NCBI_22
	use_gdt = GDT_22
        use_gi = GI_TAXID_22
        use_lut = LUT_22
        grammy_pre = GRAMMY_PRE_ACC
        refilt = REFILT
    else:
        use_db = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std.fna'
        use_gi = 'references/' + NEWDB_NAME +  '/' + NEWDB_NAME + '.acc.taxids'
        use_length = 'references/' + NEWDB_NAME + '/taxids_lengths.txt'
        use_gdt = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std'
        use_lut = 'references/' + NEWDB_NAME + '/LUT/taxids_names_lengths_tax.tab'
        grammy_pre = GRAMMY_PRE_ACC
        refilt = REFILT

if OLD_MET_REF == "TRUE":
    use_db_CT = CT_06
    use_db_GA = GA_06
    use_gdt_CT = GDT06_CT
    use_gdt_GA = GDT06_GA
    human_gis = HUMAN_06
    grammy_pre = GRAMMY_PRE_GI
    use_lut = LUT_06
    refilt = REFILT_GI



rule refilter:
	input: unfilteredgrammy = OUTPUT + '{project}/{sample}/unfiltered/{sample}.grammy.tab',
		unfilteredtblat = OUTPUT + '{project}/{sample}/unfiltered/{sample}.tblat.1',
		filteredgrammy= OUTPUT + '{project}/{sample}/C_poor/{sample}.grammy.tab',
		filteredtblat = OUTPUT + '{project}/{sample}/C_poor/{sample}.tblat.1',
	output: filtcomp = OUTPUT + '{project}/{sample}/refiltered/figures/unfilt_vs_filt_AdjBlast.png',
		refilt = OUTPUT + '{project}/{sample}/refiltered/{sample}.tblat.1'
	params: outdir = OUTPUT + '{project}/{sample}/refiltered/'
	shell:
		"""
		mkdir -p {params.outdir};
		mkdir -p {params.outdir}figures;
		if [ $(wc -l {input.filteredtblat} | cut -d' ' -f1) -gt 2 ]
		then
			python {refilt}  {use_db} {input.unfilteredgrammy} {input.unfilteredtblat} {input.filteredgrammy} {input.filteredtblat} {use_gi} {params.outdir} {output.refilt};
		else
			touch {output.filtcomp};
			touch {output.refilt}
		fi
		"""

rule grammy_clean_refiltered:
	input: cleanfa = OUTPUT + '{project}/{sample}/C_poor/{sample}.cpoor.fa',
		cleantblat = OUTPUT + '{project}/{sample}/refiltered/{sample}.tblat.1',
	output: fasta_gz = OUTPUT + '{project}/{sample}/refiltered/{sample}.fa.gz',
		fa_gz = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.fasta.gz'),
		rdt = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.rdt'),
		mtx = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.mtx'),
		lld = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.lld'),
		btp = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.btp'),
		est = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.est'),
		gra = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.gra'),
		avl = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.avl')
	resources: mem_mb = 1
	shell:
		"""
		if [ $(wc -l {input.cleantblat} | cut -d' ' -f1) -gt 1 ]
		then
			cut -f1 {input.cleantblat} | sort -u | seqtk subseq {input.cleanfa} - | gzip -1 > {output.fasta_gz};

			if [[ {OLD_MET_REF} == "FALSE" ]]
			then
				docker1 run -v `pwd`{OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy \
					bash -c "cd /data && \
						grammy_rdt -t illumina . . && \
						grammy_pre_acc -q "40,75,-5" {wildcards.sample} /{use_gdt} && \
						grammy_em -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx && \
						grammy_post {wildcards.sample}.est /{use_gdt} {wildcards.sample}.btp"
			else
				docker1 run -v `pwd`{OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy \
					bash -c "cd /data && \
							grammy_rdt -t illumina . . && \
							grammy_pre -q "40,75,-5" {wildcards.sample} /{use_gdt} && \
							grammy_em -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx && \
							grammy_post {wildcards.sample}.est /{use_gdt} {wildcards.sample}.btp"
			fi


		cp {output.fasta_gz} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.fasta_gz};
		cp {output.fa_gz} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.fa_gz};
		cp {output.rdt} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.rdt};
		cp {output.mtx} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.mtx};
		cp {output.lld} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.lld};
		cp {output.btp} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.btp};
		cp {output.est} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.est};
		cp {output.gra} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.gra};
		cp {output.avl} {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.avl};

		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.fasta_gz} {output.fasta_gz};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.fa_gz} {output.fa_gz};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.rdt} {output.rdt};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.mtx} {output.mtx};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.lld} {output.lld};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.btp} {output.btp};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.est} {output.est};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.gra} {output.gra};
		mv {OUTPUT}{wildcards.project}/{wildcards.sample}/refiltered/tmp_{output.avl} {output.avl};


		else
			touch {output.fasta_gz};
			touch {output.rdt} {output.mtx} {output.lld} {output.btp} {output.est} {output.gra} {output.avl};
		fi
		"""

rule annotate_grammy_refiltered:
	input: gra = OUTPUT + '{project}/{sample}/refiltered/{sample}.gra',
		tblat = OUTPUT + '{project}/{sample}/refiltered/{sample}.tblat.1',
		stats = OUTPUT + '{project}/{sample}/{sample}_mapping_stats.txt',
	output: tab = temp(OUTPUT + '{project}/{sample}/refiltered/{sample}.tab'),
		anno = OUTPUT + '{project}/{sample}/refiltered/{sample}.grammy.tab'
	shell:
		"""
		if [ $(wc -l {input.tblat} | cut -d' ' -f1) -gt 1 ]
		then
			Rscript {FILT_GRAM_BS} {input.gra} {output.tab} {wildcards.sample};
			Rscript {ANN_GRAM_BS} {output.tab} {input.tblat} {input.stats} {use_lut} {output.anno}
		else
			touch {output.tab}
			head -n1 {use_lut} > {output.anno}
		fi
		"""
