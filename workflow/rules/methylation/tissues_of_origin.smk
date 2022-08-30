def get_reference(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path,genome_ver, seq_type = line.strip().split('\t')
            if str(wcrds.sample) == sample_id:
                ver = genome_ver
    if ver == "hg38" and NEW_METHYL == "FALSE":
        bed = GOLDENBED_HG38
        mm = METHYLMATRIX_HG38
    if ver == "hg19" and NEW_METHYL == "FALSE":
        bed = GOLDENBED_HG19
        mm = METHYLMATRIX_HG19
    if ver == "hg38" and NEW_METHYL == "TRUE":
        bed = 'references/reference_methylomes_' + NEWMETH_NAME + '/processed_refs_hg38/golden_markers.bed'
        mm = 'references/reference_methylomes_' + NEWMETH_NAME + '/MethylMatrix_hg38/MethylMatrix.txt'
    if ver == "hg19" and NEW_METHYL == "TRUE":
        bed = 'references/reference_methylomes_' + NEWMETH_NAME + '/processed_refs_hg19/golden_markers.bed'
        mm = 'references/reference_methylomes_' + NEWMETH_NAME + '/MethylMatrix_hg19/MethylMatrix.txt'
    return[bed,mm]

rule binned_methylation:
	input: CpG_bg = OUTPUT + '{project}/{sample}/methylation_extraction/bedGraph/{sample}.bedGraph.gz',
	output: CpG_bg_tmp = temp(OUTPUT + '{project}/{sample}/tissues_of_origin/{sample}_mapped_autosomal_CpG.bedGraph.tmp'),
		intersected = temp(OUTPUT + '{project}/{sample}/tissues_of_origin/{sample}.intersect.golden.tmp'),
		binned = OUTPUT+ '{project}/{sample}/tissues_of_origin/golden_markers/{sample}'
	params: ref = get_reference
	shell:
		"""
		zcat {input.CpG_bg} | tail -n+2 | bedtools sort -i - > {output.CpG_bg_tmp};
		bedtools intersect -wo -a {params.ref[0]} -b {output.CpG_bg_tmp} -sorted | \
			awk '$6-$5==1 {{ print $0 }}' | \
			awk 'NF{{NF-=1}};1' > {output.intersected};
		Rscript {TOFAGG} {output.intersected} {output.binned};
		"""


rule tissue_of_origin:
	input: binned = OUTPUT + '{project}/{sample}/tissues_of_origin/golden_markers/{sample}',
	output: tof = OUTPUT + '{project}/{sample}/tissues_of_origin/{sample}.tsv'
	params: sum_to_one = 'FALSE',
		other = 'TRUE',
		group_by_celltype = 'FALSE',
		removals='colon_5/hsc_5/hsc_2',
		ref = get_reference
	shell:
		"""
		Rscript {TOF} \
			{input.binned} \
			{params.ref[1]} \
			{output.tof} \
			{METHYLOME_TABLE} \
			{wildcards.sample} \
			{params.sum_to_one} \
			{params.other} \
			{params.group_by_celltype} \
			{params.removals}
		"""
