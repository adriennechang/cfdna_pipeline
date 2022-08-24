rule ToBedGraph:
    input:
        'references/reference_methylomes_' + NEWMETH_NAME + '/downloaded/{methylome}'
    output:
        bedGraph='references/reference_methylomes_' + NEWMETH_NAME + '/ToBedGraph/{methylome}.bedGraph'
    threads: 20
    params:
        outdir='references/reference_methylomes_' + NEWMETH_NAME + '/ToBedGraph/'
    shell:
        """
        file_type=$(grep -P "{wildcards.methylome}\t" {METHYLOME_TABLE} | cut -f6)
        if [[ $file_type == *bw ]] || [[ $file_type == *bigWig ]]; then
            bigWigToBedGraph {input[0]} {output.bedGraph}
        fi
        if [[ $file_type == *bed.gz ]]; then
            gunzip -c {input[0]} > {output.bedGraph}
        fi
        if [[ $file_type == *bam ]] || [[ $file_type == *sra ]]; then
            rm -fr {wildcards.methylome}_temp_bam_processing
            mkdir -p {wildcards.methylome}_temp_bam_processing
            {METHEXT} --parallel 20 --bedGraph {input[0]} --output {wildcards.methylome}_temp_bam_processing --mbias_off --no_header
            zcat {wildcards.methylome}_temp_bam_processing/*bedGraph.gz > {output.bedGraph}
            rm -fr {wildcards.methylome}_temp_bam_processing
        fi
        """

# Files under build hg38 need to be scaled to hg19. files already in hg19 are
# copied over
rule liftover_human_builds:
    input:
        bedGraph='references/reference_methylomes_' + NEWMETH_NAME + '/ToBedGraph/{methylome}.bedGraph',
        tohg19='references/hg38ToHg19.over.chain.gz',
	tohg38='references/hg19ToHg38.over.chain.gz'
    output:
        hg19='references/reference_methylomes_' + NEWMETH_NAME + '/liftover_hg19/{methylome}.bedGraph.lifted',
	hg38='references/reference_methylomes_' + NEWMETH_NAME + '/liftover_hg38/{methylome}.bedGraph.lifted'
    shell:
        """
        build=$(grep -P "{wildcards.methylome}\t" {METHYLOME_TABLE} | cut -f4)
        if [[ $build == hg38 ]]
        then
            python2.7 {CROSSMAP} bed {input.tohg19} {input.bedGraph} {output.hg19}
	    cp {input.bedGraph} {output.hg38}
        else
		python2.7 {CROSSMAP} bed {input.tohg38} {input.bedGraph} {output.hg19};
		cp {input.bedGraph} {output.hg19};
	fi
			
        """

rule normalize_tissues:
    input:
        lifted = 'references/reference_methylomes_' + NEWMETH_NAME + '/liftover_{genome}/{methylome}.bedGraph.lifted'
    output:
        renamed_file = 'references/reference_methylomes_' + NEWMETH_NAME + '/normalized_{genome}/{methylome}.bedGraph'
    resources:
        mem_mb = 100
    run: #need to scale between 0 and 1... can take a while for big files so we cheat a bit for these large files
        import pandas as pd
        import numpy as np
        import os

        chunksize = 10 ** 6
        for chunk in pd.read_table(input.lifted, chunksize=chunksize, names=('chr', 'start', 'end', 'meth')):
            chunk['meth'] = chunk['meth']/np.max(chunk['meth'])
            chunk.to_csv(output.renamed_file, sep='\t', float_format='%.8f', index=False, header=False, mode='a')

rule singleBP_BedGraph: #scaleback
    input:
        normalized_file='references/reference_methylomes_' + NEWMETH_NAME + '/normalized_{genome}/{methylome}.bedGraph',
        hg19 = HG19METH + '/hg19.fa',
        hg19_lengths = HG19_CHROMO_SIZES,
	hg38 = HG38METH + '/hg38.fa',
	hg38_lengths = HG38_CHROMO_SIZES
    output:
        single_tmp = temp('references/reference_methylomes_' + NEWMETH_NAME + '/singleBP_bedGraph_{genome}/{methylome}.tmp'),
        singleBP_file ='references/reference_methylomes_' + NEWMETH_NAME + '/singleBP_bedGraph_{genome}/{methylome}.singlebp.bedGraph'
    shell:
        """
	if [[ {wildcards.genome} == hg19 ]]
	then
        	python {BPBED} {input.normalized_file} {input.hg19_lengths} {input.hg19} {output.single_tmp}
        	bedtools sort -i {output.single_tmp} | bedtools merge -i - -c 4 -o mean > {output.singleBP_file}
	else
		python {BPBED} {input.normalized_file} {input.hg38_lengths} {input.hg38} {output.single_tmp}
		bedtools sort -i {output.single_tmp} | bedtools merge -i - -c 4 -o mean > {output.singleBP_file}
	fi
        """
