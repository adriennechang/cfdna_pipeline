rule download_metholomes:
    output:
        'references/reference_methylomes_' + NEWMETH_NAME + '/downloaded/{methylome}'
    threads: 1
    shell:
        """
        sleep $((RANDOM % 100))
        filetype=$(grep -P "{wildcards.methylome}\t" {METHYLOME_TABLE} | cut -f6)
        if [[ $filetype == sra ]]; then
            bash {SRA} {BISMARK} {HG19METH} {BISDEDUP} {wildcards.methylome} {output[0]}
        else
            url=$(grep -P "{wildcards.methylome}\t" {METHYLOME_TABLE} | cut -f5)
            wget -O {output[0]} $url --no-check-certificate
        fi
        """

rule download_liftover_chain:
    output:
        'references/hg38ToHg19.over.chain.gz', 'references/hg19ToHg38.over.chain.gz'
    shell:
        """
        wget -O {output[0]} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz;
	wget -O {output[1]} http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz;
        """
