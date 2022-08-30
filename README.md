# cfDNA_pipeline
Author: Adrienne Chang  
Created: 8/22/2022

## Setup  
Before running the cfDNA pipeline:  
1. Install miniconda
2. Copy the pipeline 
3. Create the conda environment

#### *Install miniconda*
1. Change to your home directory  
`cd ~`  
2. Download and run the installation script  
`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`  
`chmod u+x Miniconda3-latest-Linux-x84_64.sh`  
`./Miniconda3-latest-Linux-x86_64.sh`

#### *Copy the pipeline*
1. Change into the shared directory and execute the copy script  
`cd /workdir/cfDNA_pipeline/`  
`sh copy_pipeline.sh <pipeline destination>`

#### *Create the conda environment*
1. Activate conda  
`source /workdir/<netid>/miniconda3/bin/activate`  
2. Create the conda environment  
`conda create -c bioconda -n cfdna_pipeline -y python=3.6 csvtk=0.24.0 fastqc=0.11.9 flash=1.2.11 seqkit=2.3.0 seqtk=1.3 taxonkit=0.12.0 trimmomatic=0.39`  
3. Activate the conda environment  
`conda activate cfdna_pipeline`  
4. Install additional packages:  
`conda install -c conda-forge -y ratelimiter r-limsolve r-matrixcalc r-matrix r-data.table zipp importlib-metadata configargparse appdirs r-dplyr`  
`conda install -c bioconda -y ucsc-bigwigtobedgraph abismal htslibb methpipe`  
`pip3 install pandas datrie pyparsing perl`  


## Running the pipeline  
*Notes:*  
 *1. Currently, the pipeline only supports analysis of standard sequenced (not bisulfite) samples. The bisulfite pipeline is coming soon!*   
 *2. The pipeline only works for 2x75 sequencing after SRSLY, MEYER, or NEXTERA library preparation.*    

The plug and play cfDNA pipeline follows the general steps:  
1. Copy data in `data/`  
2. Create a sequecing prep table  
3. Choose the metagenomic reference   
4. Execute

#### *Copy data to `data/`*  
Copy raw data to `data/<project>/` as `<sample_id>_R[1-2].fastq.gz`. It is recommended to rename the files to remove adapter information.  

#### *Create a sequencing prep table*  
 Create a new sequencing prep table (tab-delimited) using the columns below and save a copy to the `prep_tables/` folder with the file structure `sequencing_prep_<project>.tsv`:  
 - **sample_id**: Sample ID that will be used for naming (e.g., `<sample_id>_R[1-2].fastq.gz`)  
 - **project_id**: Project name  
 - **prep_type**: The library preparation type used [MEYER or SRSLY or NEXTERA]  
 - **path**: The path to the fastq location (e.g., `data/`)  
 - **genome_ver**: The genome version used for the analysis [hg19 or hg38]  
 - **seq_type**:  The sequencing type [standard or bisulfite (wgbs or SIFT-seq) or bowtie (standard-seq for comparison with bisulfite analysis)]   
    *Note: currently only standard and bowtie are supported*  
    
#### *Choose the metagenomic reference*  
The default metagenomic reference is NCBIGenomes22, which contains all reference and representative genomes (archaea, bacteria, fungi, and viral) and the human genome. There are two additional options for metagenomic reference currently supported:  
1. Using the old metagenomic reference, NCBIGenomes06:  
     To use NCBIGenomes06, open the Snakefile and change `OLD_MET_REF = "TRUE"`. **For this to work, `MAKE_NEWDB` must be set to `"FALSE"`**.  
2. Adding additional genomes to the metagenomic database:
     1. If additional genomes needed for analysis are missing from the database, set `MAKE_NEWDB = "TRUE"` and give your new database a name `NEWDB_NAME = <new_db_name>`  
     2. Create a new file `"add_assembly_accession.txt"` that contains a single column of assembly accessions (e.g., GCF_002287175.1). The assembly accessions can be obtained from the first column of an `assembly summary.txt` file downloaded from the [NCBI FTP site](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/).  
3. If you only want to make a new metagenomic reference (and not analyze any samples), change `REF_ONLY = "TRUE"`.  

#### Execute  
Execute the pipeline using snakemake: `snakemake --cores <cores>`. Additional snakemake options can be found [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options). Results are found in `results/<project><sample>`.  

------------------------------------------------------------------------
# TO DO:   
1. Test analysis pipelines
2. Incorporate rules needed to generate a MethylMatrix (testing)
