rule blast:
	input: nh1 = 'results/{project}/{sample}/{sample}.nonhost.R1.fastq.gz',
		nh2 = 'results/{project}/{sample}/{sample}.nonhost.R2.fastq.gz',
		nhu = 'results/{project}/{sample}/{sample}.nonhost_unpaired.fastq.gz',
	output: nonhostfastq = temp('results/{project}/{sample}/{sample}.nonhost.fastq'),
		nonhostfasta = temp('results/{project}/{sample}/{sample}.nonhost.fasta'),
		ob = 'results/{project}/{sample}/{sample}.blast',
		dbl = temp('results/{project}/{sample}/{sample}.doublecounts.txt'),
		tblat = 'results/{project}/{sample}/{sample}.tblat.1',
	threads: 4
	shell:
		"""
		zcat {input.nh1} {input.nh2} {input.nhu} > {output.nonhostfastq};
		fastq_to_fasta -Q33 -i {output.nonhostfastq} -o {output.nonhostfasta};
		BLASTDB=/home/id93_0001/shared/Infections/Data/NCBI/nt blastn -query {output.nonhostfasta} \
			-db {NCBI} \
			-out {output.ob} \
			-evalue 0.0001 \
			-perc_identity 90 \
			-culling_limit 1 \
			-num_threads {threads} \
			-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids scomnames sskingdoms qlen' \
			-task megablast;
		sed -i 's/\ /_/g' {output.ob};
		cat {output.ob}  | awk '{{ print $1 }}' | grep -v "\-C" | awk -F "-" '{{ print $1 }}' | \
			sort | uniq -c | awk '{{ if ($1 ==2) print $2 }}' | sed 's/$/\-2/g' > {output.dbl};
		LC_ALL=C fgrep -v -f {output.dbl} {output.ob} | sed 's/\ /\t/g' > {output.tblat};

		"""

rule grammy:
	input: nonhostfasta = 'results/{project}/{sample}/{sample}.nonhost.fasta',
		tblat1 = 'results/{project}/{sample}/{sample}.tblat.1'
	output: nhfastazip = temp('results/{project}/{sample}/{sample}.fa.gz'),
		tab = 'results/{project}/{sample}/{sample}.tab',
		mtx = 'results/{project}/{sample}/{sample}.mtx',
		est = 'results/{project}/{sample}/{sample}.est',
		btp = 'results/{project}/{sample}/{sample}.btp',
		gra = 'results/{project}/{sample}/{sample}.gra'
	threads: 4
	shell:
		"""
		cat {input.nonhostfasta} | pigz > {output.nhfastazip};
		cd results/{wildcards.project}/{wildcards.sample};
		python2.7 ../../../software/GRAMMy/grammy/grammy_rdt.py -t illumina . .;
		python2.7 ../../../software/GRAMMy/grammy/grammy_pre.py -q "40,75,-5" {wildcards.sample} {GDT};
		python2.7 ../../../software/GRAMMy/grammy/grammy_em.py -c L -b 5 -t .00001 -n 100 {wildcards.sample}.mtx;
		python2.7 ../../../software/GRAMMy/grammy/grammy_post.py {wildcards.sample}.est {GDT} {wildcards.sample}.btp;
		cd ../../../;
		cat {output.gra} | software/perl/t/transpose.pl | \
			software/perl/t/filter.pl -c 1 -mins 0 | \
			software/perl/t/add_column.pl -b -s "{wildcards.sample}" > {output.tab};
		"""
	
rule grammy_annotate:
	input: tab = 'results/{project}/{sample}/{sample}.tab',
		stat = 'results/{project}/{sample}/{sample}_stats.align.tab'
	output: 'results/{project}/{sample}/{sample}.grammy.tab'
	threads: 4
	shell: 
		"""
		Rscript scripts/annotate_grammy_main.R results/{wildcards.project}/{wildcards.sample}/ {wildcards.sample} references/NCBIGenomes06/LUTGrammy/taxids_names_lengths_tax.tab;
		"""
