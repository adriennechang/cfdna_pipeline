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
                {HSBLASTN} align -query {input.clean} \
                        -db {input.db_CT} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 > {output.out_CT};
                {HSBLASTN} align -query {input.clean} \
                        -db {input.db_GA} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 > {output.out_GA};
                """

rule filter_blast_CT:
	input: blast = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.outfmt6',
		human = human_gis,
		clean = OUTPUT + '{project}/{sample}/unfiltered/{sample}.cpoor.fa',
	output: human_blast = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.human',
		blast_filt = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.sorted.outfmt6',

	run:
		import os
                import dnaio

                cmd = f"grep -f {input.human} {input.blast} | cut -f1 | sort -u > {output.human_blast}"
                os.system(cmd)
                human_readIDS = [line.rstrip('\n') for line in open(output.human_blast)]
                with open(input.blast) as f, open(output.blast_filt, 'w') as w:
                    for hit in f:
                        read_id = hit.strip().split('\t')[0]
                        mol_len = float(hit.strip().split('\t')[3])
                        mapped_len = float(hit.strip().split('\t')[12])
                        ratio = mol_len/mapped_len
                        if (mol_len >= 50) and (ratio > 0.9):
                            if read_id not in human_readIDs:
                                    w.write(hit)

rule filter_blast_GA:
	input: blast = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.outfmt6',
		human = human_gis,
		clean = OUTPUT + '{project}/{sample}/unfiltered/{sample}.cpoor.fa',
	output: human_blast = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.human',
		blast_filt = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.sorted.outfmt6',
	run:
		import os
                import dnaio

                cmd = f"grep -f {input.human} {input.blast} | cut -f1 | sort -u > {output.human_blast}"
                os.system(cmd)
                human_readIDS = [line.rstrip('\n') for line in open(output.human_blast)]
                with open(input.blast) as f, open(output.blast_filt, 'w') as w:
                    for hit in f:
                        read_id = hit.strip().split('\t')[0]
                        mol_len = float(hit.strip().split('\t')[3])
                        mapped_len = float(hit.strip().split('\t')[12])
                        ratio = mol_len/mapped_len
                        if (mol_len >= 50) and (ratio > 0.9):
                            if read_id not in human_readIDs:
                                    w.write(hit)

rule get_relevant:
        input: CT = OUTPUT + '{project}/{sample}/unfiltered/{sample}_CT.sorted.outfmt6',
                GA = OUTPUT + '{project}/{sample}/unfiltered/{sample}_GA.sorted.outfmt6',
                clean = OUTPUT + '{project}/{sample}/unfiltered/{sample}.cpoor.fa',
                human = human_gis,
        output: clean = OUTPUT + '{project}/{sample}/unfiltered/{sample}.tblat.1'
        shell:
                """
                grep "Plus/Plus" {input.CT} > {output.clean} || true;
                grep "Plus/Minus" {input.GA} >> {output.clean} || true;
                """
rule grammy_clean:
	input: fasta =  OUTPUT + '{project}/{sample}/unfiltered/{sample}.cpoor.fa',
                tblat = OUTPUT + '{project}/{sample}/unfiltered/{sample}.tblat.1'
	output: fa = OUTPUT + '{project}/{sample}/unfiltered/{sample}.fa.gz',
                rdt = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.rdt'),
                mtx = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.mtx'),
                lld = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.lld'),
                btp = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.btp'),
                est = temp(OUTPUT + '{project}/{sample}/unfiltered/{sample}.est'),
                gra = OUTPUT + '{project}/{sample}/unfiltered/{sample}.gra',
                avl = temp(OUTPUT + '{project}/{sample}/unfilered/{sample}.avl')
	resources: mem_mb = 1
	shell:
                """
                if [ $(wc -l {input.tblat) | cut -d' ' -f1) -gt 1 ]
                then
                        cut -f1 {input.tblat} | sort -u | seqtk subseq {input.fasta} - | gzip -1 > {output.fa};
                        cd {OUTPUT}{wildcards.project}/{wildcards.sample}/unfiltered/;
                        python2.7 {GRAMMY_RDT} -t illumina . .;
                        python2.7 {grammy_pre} -q "40,75,-5" {wildcards.sample} {use_gdt_GA};
                        python2.7 {GREMMY_EM} -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx;
                        python2.7 {GRAMMY_POST} {wildcards.sample}.est {use_gdt_GA} {wildcards.sample}.btp;
                        cd ../../../../;
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
