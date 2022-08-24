if OLD_MET_REF == "FALSE":
    if MAKE_NEWDB == "FALSE":
        use_db = NCBI_22
        use_gi = GI_TAXID_22
        use_length = GI_LENGTH_22
        use_gdt = GDT_22
        use_lut = LUT_22
        filter_strand = FILTER_STRAND
        grammy_pre = GRAMMY_PRE_ACC
        use_db_CT = NCBI_CT
        use_db_GA = NCBI_GA
        use_gdt_CT = GDT_CT
        use_gdt_GA = GDT_GA
        human_gis = HUMAN_22
    else:
        use_db = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std.fna'
        use_gi = 'references/' + NEWDB_NAME +  '/' + NEWDB_NAME + '.acc.taxids'
        use_length = 'references/' + NEWDB_NAME + '/taxids_lengths.txt'
        use_gdt = 'references/' + NEWDB_NAME + '/std/' + NEWDB_NAME + '_std.gdt'
        use_lut = 'references/' + NEWDB_NAME + '/LUT/taxids_names_lengths_tax.tab'
        grammy_pre = GRAMMY_PRE_ACC
        human_gis = 'references/' + NEWDB_NAME + '/' + NEWDB_NAME + '.human.gis'

def get_conversion(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                seqt = seq_type
    if seqt == "standard":
        con = ['std']
    if seqt == "bisulfite":
        con = ['GA','CT']
    return(con)


rule hsblastn_bisulfite:
	input: i1 = OUTPUT + '{project}/{sample}/{sample}.cpoor.fa',
		gi_to_taxid = use_gi
	output: blast_GA = temp(OUTPUT + '{project}/{sample}/{sample}_GA.outfmt6'),
		blast_CT = temp(OUTPUT + '{project}/{sample}/{sample}_CT.outfmt6'),
	threads: 12
	resources: mem_mb=80000
	shell:
		"""
		{HSBLASTN} align -query {input.i1} \
				-db {use_db_GA} \
				-evalue 0.0001 \
				-perc_identity 95 \
				-num_threads {threads} \
				-outfmt 6 > {output.blast_GA};

		{HSBLASTN} align -query {input.i1} \
				-db {use_db_CT} \
				-evalue 0.0001 \
				-perc_identity 95 \
				-num_threads {threads} \
				-outfmt 6 > {output.blast_CT};
		"""



rule filter_hsblastn_bisulfite:
	input: blast_GA = temp(OUTPUT + '{project}/{sample}/{sample}_GA.outfmt6'),
                blast_CT = temp(OUTPUT + '{project}/{sample}/{sample}_CT.outfmt6'),

	
output:  tblatGA = temp(OUTPUT + '{project}/{sample}/{sample}.GA.pe'),
		tblatCT = temp(OUTPUT + '{project}/{sample}/{sample}.CT.pe'),
		tblatpe = temp(OUTPUT + '{project}/{sample}/{sample}.tblat.pe'),
		tblat1 = temp(OUTPUT + '{project}/{sample}/{sample}.tblat.1')
	threads: 10
	shell:
		"""
		filesize=$(du {input.GA_gz_R1} | cut -f1)
		if [ "$filesize" -gt 3296400 ]
		then
			bash {FILTA} {input.GA_gz_R1} {input.GA_gz_R2} {wildcards.sample}_GA {input.taxid_lengths} {threads} {output.tblatGA}
			bash {FILTA} {input.CT_gz_R1} {input.CT_gz_R2} {wildcards.sample}_CT {input.taxid_lengths} {threads} {output.tblatCT}
		else
			Rscript {FILTB} {input.GA_gz_R1} {input.GA_gz_R2} {output.tblatGA} {input.taxid_lengths}
			Rscript {FILTB} {input.CT_gz_R1} {input.CT_gz_R2} {output.tblatCT} {input.taxid_lengths}
		fi
		Rscript {CONSOLIDATE} {output.tblatGA} {output.tblatCT} {output.tblatpe} {output.tblat1}
		"""


rule hsblastn_std:
	input: i1 = OUTPUT + '{project}/{sample}/nonhuman_std/{sample}_R1.fa',
		i2 = OUTPUT + '{project}/{sample}/nonhuman_std/{sample}_R2.fa',
		gi_to_taxid = use_gi,
		db = use_db
	output: blastR1 = temp(OUTPUT + '{project}/{sample}/{sample}_std_R1.outfmt6'),
		blastR2 = temp(OUTPUT + '{project}/{sample}/{sample}_std_R2.outfmt6'),
                human_R1 = temp(OUTPUT + '{project}/{sample}/{sample}_std_R1.humanlike'),
		human_R2 = temp(OUTPUT + '{project}/{sample}/{sample}_std_R2.humanlike'),
                rejected_R1 = temp(OUTPUT + '{project}/{sample}/{sample}_std_R1.rejected'),
		rejected_R2 = temp(OUTPUT + '{project}/{sample}/{sample}_std_R2.rejected'),
		gz_R1 = OUTPUT + '{project}/{sample}/{sample}_std_R1.outfmt6.gz',
		gz_R2 = OUTPUT + '{project}/{sample}/{sample}_std_R2.outfmt6.gz'
	threads: 12
	resources: mem_mb=80000
	shell:
		"""
		{HSBLASTN} align -query {input.i1} \
			-db {input.db} \
			-evalue 0.0001 \
			-perc_identity 95 \
			-num_threads {threads} \
			-outfmt 6 | \
		python2 {filter_strand} --filename_out {output.blastR1} \
			--acc_to_tax {input.gi_to_taxid} \
		 	--conversion std \
			--human_like {output.human_R1} \
			--rejected_hits {output.rejected_R1};
		{HSBLASTN} align -query {input.i2} \
                        -db {input.db} \
                        -evalue 0.0001 \
                        -perc_identity 95 \
                        -num_threads {threads} \
                        -outfmt 6 | \
                python2 {filter_strand} --filename_out {output.blastR2} \
                        --acc_to_tax {input.gi_to_taxid} \
                        --conversion std \
                        --human_like {output.human_R2} \
                        --rejected_hits {output.rejected_R2};
		gzip -c {output.blastR1} > {output.gz_R1};
		gzip -c {output.blastR2} > {output.gz_R2};
		"""

rule filter_blastn_stds:
	input: gz_R1 = OUTPUT + '{project}/{sample}/{sample}_std_R1.outfmt6.gz',
                gz_R2 = OUTPUT + '{project}/{sample}/{sample}_std_R2.outfmt6.gz',
		taxid_lengths = use_length
	output: tblatpe = temp(OUTPUT + '{project}/{sample}/{sample}.std.tblat.pe'),
                tblat1 = temp(OUTPUT + '{project}/{sample}/{sample}.std.tblat.1')
	threads: 10
	shell:
		"""
		filesize=$(du {input.gz_R1} | cut -f1)
		if [ $filesize -gt 3696400 ]
		then
			bash {FILTA} {input.gz_R1} {input.gz_R2} {wildcards.sample]_std {input.taxid_lengths} {threads} {output.tblatpe}
		else
			Rscript {FILTB} {input.gz_R1} {input.gz_R2} {output.tblatpe} {input.taxid_lengths}
		fi
		cat {output.tblatpe} | grep -v qseqid | awk '{{ print $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $3, $1 }}' | tr ' ' '\t' > {output.tblat1};
		"""

def aggregate_tblat_files(wildcards):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                seq = seq_type
    if seq_type == "bisulfite":
        nonhu_r1 = OUTPUT + '{project}/{sample}/nonhuman/{sample}_R1.fa'
        nonhu_r2 = OUTPUT + '{project}/{sample}/nonhuman/{sample}_R2.fa'
        tblatpe = OUTPUT + '{project}/{sample}/{sample}.tblat.pe'
        tblat1 = OUTPUT + '{project}/{sample}/{sample}.tblat.1'
    if seq_type == "standard":
        nonhu_r1 = OUTPUT + '{project}/{sample}/nonhuman_std/{sample}_R1.fa'
        nonhu_r2 = OUTPUT + '{project}/{sample}/nonhuman_std/{sample}_R2.fa'
        tblatpe = OUTPUT + '{project}/{sample}/{sample}.std.tblat.pe'
        tblat1 = OUTPUT + '{project}/{sample}/{sample}.std.tblat.1'
    return[nonhu_r1, nonhu_r2, tblatpe, tblat1]


rule aggregate_tblats:
	input: aggregate_tblat_files
	output: nonhu_r1 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R1.fa',
		nonhu_r2 = OUTPUT + '{project}/{sample}/{sample}.nonhost.R2.fa',
		tblat1 = OUTPUT + '{project}/{sample}/grammy/{sample}.tblat.1',
		tblatpe = OUTPUT + '{project}/{sample}/grammy/{sample}.tlbat.pe'
	shell:
		"""
		cp {input[0]} {output.nonhu_r1};
		cp {input[1]} {output.nonhu_r2};
		cp {input[2]} {output.tblatpe};
		cp {input[3]} {output.tblat1};
		"""
rule filter_contam:
	input: tblatpe = OUTPUT + '{project}/{sample}/grammy/{sample}.tblat.pe',
		reference_fasta = use_db,
		LUT = use_lut,
		unmap_R1 = OUTPUT + '{project}/{sample}/aligned/unmapped/{sample}_unmapped_R1.fastq.gz',
		unmap_R2 = OUTPUT + '{project}/{sample}/aligned/unmapped/{sample}_unmapped_R2.fastq.gz',
		model1 = CSFMOD1,
		model2 = CFSMOD2
	output: names=temp(OUTPUT + '{project}/{sample}/{sample}.reads'),
		R1_sort = temp(OUTPUT + '{project}/{sample}/{sample}_temp_sort_R1.fastq'),
		R2_sort = temp(OUTPUT + '{project}/{sample}/{sample}_temp_sort_R2.fastq'),
		tlbat = temp(OUTPUT + '{project}/{sample}/grammy/{sample}.sort.tblat.pe'),
		removed = OUTPUT + '{project}/{sample}/filt_grammy/{sample}.removed.blast',
		keep = OUTPUT + '{project}/{sample}/filt_grammy/{sample}.blast',
		tblat = OUTPUT + '{project}/{sample}/filt_grammy/{sample}.tblat.1'
	params: dir = OUTPUT + '{project}/{sample}/filt/grammy/'
	log: 'logs/{project}/{sample}.contamination_filtering.log'
	threads: 1
	shell:
		"""
		cut -f2 {input.tblatpe} | tail -n+2 | cut -f2 -d'_' > {output.names};
		zcat {input.unmap_R1} | cut -f1 -d'_' | cut -f1 -d' ' | \
			seqtk subseq - {output.names} | \
			paste - - - - | sort -k1,1 -S 3G -T ./ | tr '\t' '\n' > {output.R1_sort};
		zcat {input.unmap_R2} | cut -f1 -d'_' | cut -f1 -d' ' | \
			seqtk subseq - {output.names} | \
			paste - - - - | sort -k1,1 -S 3G -T ./ | tr '\t' '\n' > {output.R2_sort};

		(head -n1 {input.tblatpe} && tail -n+2 {input.tblatpe} | sort -k2 -T ./) > {output.tlbat};
		python
