if OLD_MET_REF == "FALSE":
    if MAKE_NEWDB == "FALSE":
        use_db_CT = NCBI_CT
        use_db_GA = NCBI_GA
        use_gdt_CT = GDT_CT
        use_gdt_GA = GDT_GA
        human_gis = HUMAN_22
        use_gi = GI_TAXID_22
        grammy_pre = GRAMMY_PRE_ACC
        use_lut = LUT_22
    else:
        use_db_CT = 'references/' + NEWDB_NAME + '/CT/' + NEWDB_NAME + '_CT.fna'
        use_db_GA = 'references/' + NEWDB_NAME + '/GA/' + NEWDB_NAME + '_GA.fna'
        use_gdt_CT = 'references/' + NEWDB_NAME + '/CT/' + NEWDB_NAME + '_CT'
        use_gdt_GA = 'references/' + NEWDB_NAME + '/GA/' + NEWDB_NAME + '_GA'
        use_gi = 'references/' + NEWDB_NAME + '/' + NEWDB_NAME + '.acc.taxids'
        human_gis = 'references/' + NEWDB_NAME + '/' + NEWDB_NAME + '.human.gis'
        grammy_pre = GRAMMY_PRE_ACC
        use_lut = 'references/' + NEWDB_NAME + '/LUT/taxids_names_lengths_tax.tab'

if OLD_MET_REF == "TRUE":
    use_db_CT = CT_06
    use_db_GA = GA_06
    use_gdt_CT = GDT06_CT
    use_gdt_GA = GDT06_GA
    human_gis = HUMAN_06
    grammy_pre = GRAMMY_PRE_GI
    use_lut = LUT_06



rule hs_blastn_wgbs:
	input: fa = OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.CT.fa',
		db_CT = use_db_CT,
		db_GA = use_db_GA,
		gi_tax = use_gi,
	output: out_GA = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.outfmt6',
		out_CT = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.outfmt6'	
	threads: 8
	resources: mem_mb=80000
	shell:
		"""
                {HSBLASTN} align -query {input.fa} \
                        -db {input.db_CT} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 > {output.out_CT};
                {HSBLASTN} align -query {input.fa} \
                        -db {input.db_GA} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 > {output.out_GA};
                """

rule filter_blast_CT:
	input: blast = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.outfmt6',
		clean = OUTPUT + '{project}/{sample}/C_poor/{sample}.cpoor.fa',
	output: human_blast = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.human'),
		blast_filt = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.sorted.outfmt6',

	run:
		import os
                import dnaio
		import os.path
	
		if os.path.exists(human_gis):
		    cmd = f"grep -f {human_gis} {input.blast} | cut -f1 | sort -u > {output.human_blast}"
		    os.system(cmd)
		    human_readIDs = [line.rstrip('\n') for line in open(output.human_blast)]
		else:
		    cmd = f" echo 'empty' > {output.human_blast}"
		    os.system(cmd)
		    human_readIDs=[]

		with open(input.blast) as f, open(output.blast_filt, 'w') as w:
		    for hit in f:
		        read_id = hit.strip().split('\t')[0]
		        mol_len = float(hit.strip().split('\t')[3])
		        mapped_len = float(hit.strip().split('\t')[12])
		        ratio = mol_len/mapped_len
		        if (mol_len >= 50) and (ratio > 0.9):
		            if human_readIDs != [] and read_id not in human_readIDs:
		                w.write(hit)
		            if human_readIDs == []:
		                w.write(hit)

rule filter_blast_GA:
	input: blast = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.outfmt6',
		clean = OUTPUT + '{project}/{sample}/C_poor/{sample}.cpoor.fa',
	output: human_blast = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.human'),
		blast_filt = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.sorted.outfmt6',
	run:
		import os
                import dnaio
		import os.path

		if os.path.exists(human_gis):
		    cmd = f"grep -f {human_gis} {input.blast} | cut -f1 | sort -u > {output.human_blast}"
		    os.system(cmd)
		    human_readIDs = [line.rstrip('\n') for line in open(output.human_blast)]
		else:
		    cmd = f"echo 'empty' > {output.human_blast}"
		    os.system(cmd)
		    human_readIDs = []

		with open(input.blast) as f, open(output.blast_filt, 'w') as w:
		    for hit in f:
		        read_id = hit.strip().split('\t')[0]
		        mol_len = float(hit.strip().split('\t')[3])
		        mapped_len = float(hit.strip().split('\t')[12])
		        ratio = mol_len/mapped_len
		        if (mol_len >= 50) and (ratio > 0.9):
		            if human_readIDs != [] and read_id not in human_readIDs:
		                w.write(hit)
		            if human_readIDs == []:
		                w.write(hit)

rule get_relevant:
        input: CT = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.sorted.outfmt6',
                GA = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.sorted.outfmt6',
        output: clean = OUTPUT + '{project}/{sample}/unfiltered/{sample}.tblat.1'
        shell:
                """
                grep "Plus/Plus" {input.CT} > {output.clean} || true;
                grep "Plus/Minus" {input.GA} >> {output.clean} || true;
                """
rule grammy_clean:
	input: fasta =  OUTPUT + '{project}/{sample}/decontaminate/{sample}.decon.CT.fa',
                tblat = OUTPUT + '{project}/{sample}/unfiltered/{sample}.tblat.1'
	output: fa = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.fa.gz'),
		fasta = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.fasta.gz'),
                rdt = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.rdt'),
                mtx = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.mtx'),
                lld = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.lld'),
                btp = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.btp'),
                est = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.est'),
                gra = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.gra'),
                avl = temp(OUTPUT + '{project}/{sample}/unfilered/{sample}.avl')
	resources: mem_mb = 1
	shell:
                """
                if [ $(wc -l {input.tblat} | cut -d' ' -f1) -gt 1 ]
                then
                        cut -f1 {input.tblat} | sort -u | seqtk subseq {input.fasta} - | gzip -1 > {output.fa};
    
			if [[ {OLD_MET_REF} == "FALSE" ]]
			then
				docker1 run -v `pwd`/{OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy \
					bash -c "cd /data &&  \
						grammy_rdt -t illumina . . && \
						grammy_pre_acc -q "40,75,-5" {wildcards.sample} /{use_gdt_GA} \
						grammy_em -c L -b 5 -t .00001 -n 100 {wildcarsd.sample}.mtx && \
						grammy_post {wildcards.sample}.est /{use_gdt_GA} {wildcards.sample}.btp"
			else
				docker1 run -v `pwd`/{OUTPUT}/{wildcards.project}/{wildcards.sample}/unfiltered:/data -v `pwd`/main_reference:/references biohpc_ac2763/grammy \
					bash -c " cd /data && \
						grammy_rdt -t illumina . . && \
						grammy_pre -q "40,75,-5" {wildcards.sample} /{use_gdt_GA} \
						grammy_em -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx && \
						grammy_post {wildcards.sample}.est /{use_gdt_GA} {wildcards.sample}.btp"
			fi
		
			cp {output.fa} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.fa};
			cp {output.fasta} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.fasta};
			cp {output.rdt} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.rdt};
			cp {output.mtx} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.mtx};
			cp {output.est} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.est};
			cp {output.lld} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.lld};
			cp {output.btp} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.btp};
			cp {output.gra} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.gra};
			cp {output.avl} {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.avl};

			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.fa} {output.fa};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.fasta} {output.fasta};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.rdt} {output.rdt};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.lld} {output.lld};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.mtx} {output..mtx};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.btp} {output.btp};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.est} {output.est};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.gra} {output.gra};
			mv {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/tmp_{output.avl} {output.avl};

                else
                       touch {output.fa} {output.rdt} {output.mtx} {output.lld} {output.btp} {output.est} {output.gra} {output.avl};
                fi
                """

rule annotate_grammy_bs:
	input: gra = OUTPUT + '{project}/{sample}/unfiltered/{sample}.gra',
                tblat = OUTPUT + '{project}/{sample}/unfiltered/{sample}.tblat.1',
                stats = OUTPUT + '{project}/{sample}/{sample}_mapping_stats.txt'
	output: tab = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.tab'),
                anno = OUTPUT + '{project}/{sample}/unfiltered/{sample}.grammy.tab'
	shell:
                """
                if [ $(wc -l {input.tblat} | cut -d' ' -f1) -gt 2 ]
                then
                        Rscript {FILT_GRAM_BS} {input.gra} {output.tab} {wildcards.sample};
                        Rscript {ANN_GRAM_BS} {output.tab} {input.tblat} {input.stats} {use_lut} {output.anno}
                else
                        touch {output.tab}
                        head -n1 {use_lut} > {output.anno}
                fi
                """
