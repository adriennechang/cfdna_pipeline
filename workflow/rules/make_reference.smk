REF_DIR = 'references/' + NEWDB_NAME + '/'

rule find_genomes:
	input: insc = DL_GENOMES,
	output: summary = REF_DIR + 'assembly_summary.txt',
		allsum = REF_DIR + 'all_assembly_summmary.txt'
	shell:
		"""
		mkdir -p references/;
		mkdir -p {REF_DIR};
		mkdir -p {REF_DIR}genomes;
		sh {input.insc} {REF_DIR} {output.summary} {output.allsum};
		"""

rule add_genomes:
	input: newfile = 'add_assembly_acession.txt',
		allsum = REF_DIR + 'all_assembly_summmary.txt' 
	output: summary = REF_DIR + 'assembly_summary.txt'
	shell:
		"""
		LC_ALL=C fgrep -f {input.newfile} {input.allsum} >> {output.summary};
		"""
	

rule download_genomes:
	input: summary = REF_DIR + 'assembly_summary.txt',
	output: dl = temp(REF_DIR + 'genomes/download.sh'),
                tempin = temp(REF_DIR + 'tmp_' + NEWDB_NAME + '.fna'),
                dupids = temp(REF_DIR + 'genomes/duplicate_ids.txt'),
                infa = REF_DIR + NEWDB_NAME + '.fna'
	shell:
		"""
		awk -F "\t" '{{ print $20 }}' {input.summary} | \
			awk -F "/" '{{ print "wget -P genomes/ "$0"/"$NF"_genomic.fna.gz" }}' > {REF_DIR}genomes/download.sh;
		cd {REF_DIR};
		sh genomes/download.sh;
		cd ../../
		find {REF_DIR}genomes/ -name "*.gz" | xargs zcat > {output.tempin};
		seqkit rmdup --id-ncbi -D {output.dupids} < {output.tempin} > {output.infa};
		touch {output.dupids};
		find {REF_DIR}genomes/ -name "*.gz" -exec rm {{}} \;
		"""

rule download_acc:
	input: infa = REF_DIR + NEWDB_NAME + '.fna'
	output: gz1 = temp(REF_DIR + 'nucl_wgs_accession2taxid.gz'),
		gz2 = temp(REF_DIR + 'nucl_gb.accession2taxid.gz'),
		acc = REF_DIR + 'accession2taxid.txt',
		acc2 = REF_DIR + NEWDB_NAME + '.acc',
		acctax = REF_DIR + NEWDB_NAME + '.acc.taxids',
		human = REF_DIR + NEWDB_NAME + '.human.gis'
	shell:
		"""
		wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz -O {output.gz1};
		wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz -O {output.gz2};
		gunzip -c {output.gz1} > {output.acc};
		gunzip -c {output.gz2} >> {output.acc};
		LC_ALL=C fgrep ">" {input.infa} | \
			awk '{{ print $1 }}' | \
			sed 's/>//g' > {output.acc2};
		LC_ALL=C fgrep -f {output.acc2} {output.acc} > {output.acctax};
		awk '$3=="9606" {{print $0}}' {output.acctax} > {output.human};
		"""

rule convert_fasta:
	input: infa = REF_DIR + NEWDB_NAME + '.fna'
	output: CT_conv = REF_DIR + 'CT/' + NEWDB_NAME + '_CT.fna',
		GA_conv = REF_DIR + 'GA/' + NEWDB_NAME + '_GA.fna'
	shell:
		"""
		python {REF_CONV} {input.infa} {output.CT_conv} {output.GA_conv};
		"""

rule copy_std_fasta:
	input: infa = REF_DIR + NEWDB_NAME + '.fna'
	output: out = REF_DIR + 'std/' + NEWDB_NAME + '_std.fna',
	shell:
		"""
		cp {input.infa} {output.out};
		"""

rule make_blastdb:
	input: acctax = REF_DIR + NEWDB_NAME + '.acc.taxids',
		acc = REF_DIR + NEWDB_NAME + '.acc',
		infa = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.fna',
	output:	counts = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.fna.counts',
		obinary = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.fna.counts.obinary',
	params: genome = NEWDB_NAME+'_{conversion}',
		acctax = NEWDB_NAME+'.acc.taxids',
		acc = NEWDB_NAME+'.acc'
	conda: 'workflow/envs/python2.yml'
	shell:
		"""
		cd {REF_DIR}{wildcards.conversion};
		makeblastdb -in {params.genome}.fna -out {params.genome}.fna -dbtype nucl -parse_seqids -taxid_map ../{params.acctax};
		blastdb_aliastool -db {params.genome}.fna -seqidlist ../{params.acc} -dbtype nucl -out {params.genome}.curated;
		windowmasker -in {params.genome}.fna -infmt blastdb -mk_counts -out {params.genome}.fna.counts;
		windowmasker -in {params.genome}.fna.counts -sformat obinary -out {params.genome}.fna.counts.obinary -convert;
		{HSBLASTN} index {params.genome}.fna;
		cd ../../../;
		"""

rule index_reference:
	input: infa = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.fna',
	output: fai = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.fna.fai'
	shell:
		"""
		samtools faidx {input.infa};
		"""

rule get_genome_sizes:
	input: fai = REF_DIR + 'std/' + NEWDB_NAME + '_std.fna.fai',
		acctax = REF_DIR + NEWDB_NAME + '.acc.taxids'
	output: lengths = REF_DIR + 'taxids_lengths.txt'
	shell:
		"""
		Rscript {GENOMESIZE} {input.fai} {input.acctax} {output.lengths};
		"""

rule make_grammydb:
	input: genome = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.fna',
		obinary = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.fna.counts.obinary',
		acctax = REF_DIR + NEWDB_NAME + '.acc.taxids',
	output: prompts = REF_DIR + '{conversion}/grammy_prompts_{conversion}',
		taxids = REF_DIR + 'grefs/{conversion}/gid_tid.dmp',
		gdt = REF_DIR + '{conversion}/' + NEWDB_NAME + '_{conversion}.gdt'
	params: genome = NEWDB_NAME +'_{conversion}'
	threads: 40
	conda: 'workflow/envs/python2.yml'
	shell:
		"""
		rm -f {output.prompts};
		mkdir -p {REF_DIR}grefs/;
		mkdir -p {REF_DIR}grefs/{wildcards.conversion};
		blastdbcmd -db {REF_DIR}{wildcards.conversion}/{params.genome}.curated -entry all -outfmt "%a" | sort | uniq | LC_ALL=C fgrep -f - {input.acctax} | cut -f3 | sort | uniq | while read taxid; do echo "sh workflow/scripts/build.grefs.sh $taxid {input.acctax} {wildcards.conversion} {REF_DIR}{wildcards.conversion}/{params.genome}.curated {REF_DIR}grefs/{wildcards.conversion}" >> {output.prompts}; done
		perl_fork_univ.pl {output.prompts} {threads};
		cat {input.acctax} | awk '{{ print $2 "\t" $3 }}' > {output.taxids};
		TIDS=$(cat {input.acctax} | cut -f3 | sort -u | tr '\n' ',');
		TAXIDS=${{TIDS::-1}};
		python2.7 {GRAMMY_GDT} -p 200000 -d {output.taxids} -r {REF_DIR}grefs/{wildcards.conversion} {params.genome} $TAXIDS
		mv {params.genome}.gdt {output.gdt};
		"""
 
rule LUTfile:
	input: acctax = REF_DIR + NEWDB_NAME + '.acc.taxids',
		lengths = REF_DIR + 'taxids_lengths.txt'
	output: taxonly = temp(REF_DIR + 'LUT/taxids.txt'),
		tax_name = REF_DIR + 'LUT/taxid_names.tab',
		taxgenus = temp( REF_DIR + 'LUT/taxids_genus.tab'),
		taxnamelengthtax = REF_DIR + 'LUT/taxids_names_lengths_tax.tab'
	conda: 'workflow/envs/python3.yml'
	shell:
		"""
		echo "sh {MAKELUT} {input.acctax} {output.taxonly} {output.tax_name} {output.taxgenus}" | bash;
		Rscript {MERGELUT} {output.tax_name} {input.lengths} {output.taxgenus} {output.taxnamelengthtax};
		"""
