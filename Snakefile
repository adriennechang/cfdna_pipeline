##############################
###   cfDNA pipeline v2    ###
### Author: Adrienne Chang ###
###   Created: 8/22/2022   ###
##############################

import os
import pandas as pd
from natsort import natsorted
import itertools


### Updating the metagenomic database ###
# A new metagenomic database will only be created if
# 	(1) MAKE_NEWDB is set to "TRUE", and
# 	(2) A new file "add_assembly_accession.txt" is created.
# 	    The new file should contain a single column of assembly accessions (i.e. GCF_002287175.1).
# 	    This can be obtained from the first column of the assembly_summary.txt file downloaded
# 	    from ftp://ftp.ncbi.nlm.nih.gov/, and
# 	(3) A NEWDB_NAME is set
MAKE_NEWDB = "FALSE"
NEWDB_NAME = "NCBIGenomes22_human_virus"


### CHANGE ONLY IF USING NCBIGENOMES06 ###
### CANNOT BE USED WITH MAKE_NEWDB ###
OLD_MET_REF = "FALSE"

### CHANGE ONLY IF MAKING A NEW METHYLMATRIX ###
NEW_METHYL = "FALSE"
NEWMETH_NAME = "meth_test"

### CHANGE ONLY IF YOU ONLY WANT TO MAKE A NEW DB"
REF_ONLY = "FALSE"

## METAGENOMIC RE-ESTIMATION USING KALLISTO? ##
KALLISTO = "FALSE"

#####################################
### DO NOT MODIFY BELOW THIS LINE ###
#####################################
if MAKE_NEWDB == "TRUE" and OLD_MET_REF == "TRUE":
    sys.exit("Error: cannot make a new database and analyze samples using an old reference.")

### LOAD CONFIG ###
configfile: 'workflow/main.yaml'
DATA=config['DATA']
MP = config['MASTER_PATH']
SEQ_PREP_TABLE = config['SEQ_PREP_TABLE']
METHYLOME_TABLE = config['METHYLOME_TABLE']
OUTPUT = config['OUTPUT']
GENOMES = config['GENOMES']
CHR = config['CHR']

# References
KALIDX = MP + config['KALIDX']
HG19 = MP + config['HG19']
HG19_CHROMO_SIZES = MP + config['HG19_CHROMO_SIZES']
HG38 = MP + config['HG38']
HG38_CHROMO_SIZES = MP + config['HG38_CHROMO_SIZES']
PHIX = MP + config['PHIX']
HG19STD = MP + config['HG19STD']
HG19METH = MP + config['HG19METH']
HG38STD = MP + config['HG38STD']
HG38METH = MP + config['HG38METH']
PHIXSTD = MP + config['PHIXSTD']
PHIXMETH = MP + config['PHIXMETH']

REFERENCE_METHYLOMES = config['REFERENCE_METHYLOMES']

NCBI_06 = MP + config['NCBI_06']
GI_TAXID_06 = MP + config['GI_TAXID_06']
GI_LENGTH_06 = MP + config['GI_LENGTH_06']
GDT_06 =  config['GDT_06']
LUT_06 = MP + config['LUT_06']
GDT06_CT =  config['GDT06_CT']
GDT06_GA =  config['GDT06_GA']
HUMAN_06 = MP +config['HUMAN_06']
NCBI_22 = MP + config['NCBI_22']
GI_TAXID_22 = MP + config['GI_TAXID_22']
GI_LENGTH_22 = MP + config['GI_LENGTH_22']
GDT_22 = config['GDT_22']
LUT_22 = MP + config['LUT_22']
NCBI_CT = MP + config['NCBI_CT']
NCBI_GA = MP + config['NCBI_GA']
CT_06 = MP + config['CT_06']
GA_06 = MP + config['GA_06']
GDT_CT = config['GDT_CT']
GDT_GA = config['GDT_GA']
HUMAN_22 = MP + config['HUMAN_22']
NCBI_06_CT = MP + config['NCBI_06_CT']
NCBI_06_GA = MP + config['NCBI_06_GA']

GOLDENBED_HG19 = MP + config['GOLDENBED_HG19']
GOLDENBED_HG38 = MP + config['GOLDENBED_HG38']
METHYLMATRIX_HG19 = MP + config['METHYLMATRIX_HG19']
METHYLMATRIX_HG38 = MP + config['METHYLMATRIX_HG38']

CT_DECON = MP + config['CT_DECON']
GA_DECON = MP + config['GA_DECON']
BOTH_DECON = MP + config['BOTH_DECON']

SRSLY = MP + config['SRSLY']
MEYER = MP + config['MEYER']
NEXTERA = MP + config['NEXTERA']

CFSMOD1 = MP + config['CFSMOD1']
CFSMOD2 = MP + config['CFSMOD2']

# Scripts
REFILT = config['REFILT']
REFILT_GI = config['REFILT_GI']
FILT_GRAM_BS = config['FILT_GRAM_BS']
ANN_GRAM_BS = config['ANN_GRAM_BS']
AGGMETH = config['AGGMETH']
BPBED = config['BPBED']
RMDBL = config['RMDBL']
DL_GENOMES = config['DL_GENOMES']
REF_CONV = config['REF_CONV']
GENOMESIZE = config['GENOMESIZE']
MAKELUT = config['MAKELUT']
MERGELUT = config['MERGELUT']
FILTER_STRAND = config['FILTER_STRAND']
FILTER_STRAND_GI = config['FILTER_STRAND_GI']
COUNT_FASTQ = config['COUNT_FASTQ']
CALC_COV = config['CALC_COV']
NON_DUP = config['NON_DUP']
INTER = config['INTER']
ANN_GRAMMY_GI = config['ANN_GRAMMY_GI']
ANN_GRAMMY_ACC = config['ANN_GRAMMY_ACC']
SRA = config['SRA']
TOFAGG = config['TOFAGG']
TOF = config['TOF']
INSIL_CONV = config['INSIL_CONV']
FILTA = config['FILTA']
FILTB = config['FILTB']
CONSOLIDATE = config['CONSOLIDATE']

# Programs
BBDUK = MP + config['BBDUK']
DEDUP2 = MP + config['DEDUP2']
ADD_STR = MP + config['ADD_STR']
F2S = MP + config['F2S']
S2F = MP + config['S2F']
TRANSPOSE = MP + config['TRANSPOSE']
FILTER = MP + config['FILTER']
ADD_COL = MP + config['ADD_COL']
HSBLASTN = MP + config['HSBLASTN']
GRAMMY_GDT = MP + config['GRAMMY_GDT']
GRAMMY_EM = MP + config['GRAMMY_EM']
GRAMMY_POST = MP + config['GRAMMY_POST']
GRAMMY_RDT = MP + config['GRAMMY_RDT']
GRAMMY_PRE_ACC = MP + config['GRAMMY_PRE_ACC']
GRAMMY_PRE_GI = MP + config['GRAMMY_PRE_GI']
BISMARK = MP + config['BISMARK']
BISDEDUP = MP + config['BISDEDUP']
METHEXT = MP + config['METHEXT']
RMDUPS = MP + config['RMDUPS']
DEDUP = MP + config['DEDUP']
CROSSMAP = MP + config['CROSSMAP']
METILENE = MP + config['METILENE']

### LOAD UNIVERSAL VARIABLES ###
sample_information  = pd.read_csv("sequencing_prep.tsv", index_col = 0, sep = "\t")

def get_project(sample): return(sample_information.loc[sample, "project_id"])
def get_fastq_path(sample): return(sample_information.loc[sample, "path"])
def get_prep_type(sample): return(sample_information.loc[sample,"seq_type"])

SAMPLES = sample_information.index.values
PROJECTS = list(set([get_project(s) for s in SAMPLES]))
samplesInProjects = {}
for p in PROJECTS:
    samplesInProjects[p] = sample_information.loc[sample_information.project_id == p].index.tolist()
SEQTYPES = list(set([get_prep_type(s) for s in SAMPLES]))
samplesInSeqs = {}
for t in SEQTYPES:
    samplesInSeqs[t] = sample_information.loc[sample_information.seq_type== t].index.tolist()

def get_sample_info(wcrds):
    with open(SEQ_PREP_TABLE) as f:
        next(f)
        for line in f:
            sample_id, project_id, prep_type, path, genome_ver, seq_type = line.strip().split('\t')
        if str(wcrds.sample) == sample_id:
            return[sample_id, project_id, prep_type, path, genome_ver, seq_type]

### FUNCTIONS ###
def string_split(string,delimiter, number_of_splits):
    s=string.split(delimiter, number_of_splits)
    return(s)

tissues = natsorted(set(x.split('_')[0] for x in config['REFERENCE_METHYLOMES']))
comparing_groups=list(itertools.combinations(tissues,2))
comparing_groups=[x[0]+'_'+x[1] for x in comparing_groups]

### RULE SETS ###
MAKEREF = [expand('references/' + NEWDB_NAME + '/' + 'LUT/taxids_names_lengths_tax.tab', conversion=config['CONVERSION']),expand('references/' + NEWDB_NAME + '/{conversion}/' + NEWDB_NAME + '_{conversion}.gdt', conversion = config['CONVERSION'])]

MAKE_METHYL = [expand('references/reference_methylomes_' + NEWMETH_NAME + '/singleBP_bedGraph_{genome}/{methylome}.singlebp.bedGraph', comp_group = comparing_groups, methylome = config['REFERENCE_METHYLOMES'], genome=config['GENOMES'])]


# Pick pipeline based on sample prep
if 'standard' in samplesInSeqs.keys(): 
	if KALLISTO == "FALSE":
		STD_ANALYSIS = [ OUTPUT + '{project}/{sample}/{sample}.grammy.tab'.format(sample=sample, project=get_project(sample)) for sample in samplesInSeqs['standard']]
	if KALLISTO == "TRUE":
		STD_ANALYSIS = [ OUTPUT + '{project}/{sample}/{sample}.bam'.format(sample=sample, project=get_project(sample)) for sample in samplesInSeqs['standard']]
else:
	STD_ANALYSIS = []

if 'bisulfite' in samplesInSeqs.keys():
	METH_ANALYSIS = [OUTPUT + '{project}/{sample}/refiltered/{sample}.grammy.tab'.format(sample=sample,project=get_project(sample)) for sample in samplesInSeqs['bisulfite']]
else:
	METH_ANALYSIS = []

if 'bowtie' in samplesInSeqs.keys():
	BIS_ANALYSIS = [OUTPUT + '{project}/{sample}/blast/{sample}.grammy.tab'.format(sample=sample, project=get_project(sample)) for sample in samplesInSeqs['bowtie']]
else:
	BIS_ANALYSIS = []

ANALYSIS = STD_ANALYSIS + METH_ANALYSIS + BIS_ANALYSIS

# Pick pipeline based on reference

if MAKE_NEWDB == "FALSE" and NEW_METHYL == "FALSE":
    inputs = ANALYSIS
elif MAKE_NEWDB == "TRUE" and NEW_METHYL == "FALSE" and REF_ONLY == "FALSE":
    inputs = MAKEREF + ANALYSIS 
elif MAKE_NEWDB == "TRUE" and NEW_METHYL == "FALSE" and REF_ONLY == "TRUE":
    inputs = MAKEREF
elif MAKE_NEWDB == "FALSE" and NEW_METHYL == "TRUE":
    inputs = MAKE_METHYL 
elif MAKE_NEWDB == "TRUE" and NEW_METHYL == "TRUE":
    inputs = MAKEREF + MAKE_METHYL
else:
    inputs = MAKEREF + MAKE_METHYL + ANALYSIS
    
    
rule all:
	input: inputs


### CONDITION TO CONVERT .fastq TO _fastq ###
if (os.path.exists(DATA + '{project}/{sample}_R1.fastq.gz')):
	ruleorder: touch > fastqc
elif (os.path.exists(DATA + '{project}/{sample}.R1.fastq.gz')):
	ruleorder: rename > fastqc

### CONDITION TO MAKE NEW REF BEFORE BLAST ###
if MAKE_NEWDB == "TRUE":
    ruleorder: LUTfile > make_grammydb > hsblastn

### RULES TO GENERATE NEW METAGENOMIC DB ###
include: 'workflow/rules/make_reference.smk'

### RULES TO GENERATE NEW TOF DB ###
include: 'workflow/rules/reference_methylomes/download_methylomes.smk'
include: 'workflow/rules/reference_methylomes/format_methylomes.smk'
include: 'workflow/rules/reference_methylomes/create_methylmatrix.smk'
#include: 'workflow/rules/reference_methylomes/get_DMRs2.smk'

### RULES TO RUN SAMPLE ANALYSIS ###
include: 'workflow/rules/rename.smk'
include: 'workflow/rules/fastqc.smk'
include: 'workflow/rules/trim/trimmomatic.smk'
include: 'workflow/rules/merge.smk'
include: 'workflow/rules/alignment/bwa.smk'
include: 'workflow/rules/stats/align_stats.smk'
#include: 'workflow/rules/metagenomic/hsblastn.smk'
include: 'workflow/rules/metagenomic/hsblastn_docker.smk'
include: 'workflow/rules/metagenomic/kallisto.smk'

### BISULFITE ONLY RULES ###
include: 'workflow/rules/trim/bbduk.smk'
include: 'workflow/rules/alignment/bismark.smk'
include: 'workflow/rules/methylation/estimate_BSconversion.smk'
include: 'workflow/rules/stats/mapping_stats.smk'
include: 'workflow/rules/methylation/methylation_extraction.smk'
include: 'workflow/rules/methylation/tissues_of_origin.smk'
include: 'workflow/rules/metagenomic/decontaminate.smk'
#include: 'workflow/rules/metagenomic/C_poor_abundances.smk'
include: 'workflow/rules/metagenomic/C_poor_abundances_docker.smk'
#include: 'workflow/rules/metagenomic/unfiltered_abundances.smk'
include: 'workflow/rules/metagenomic/unfiltered_abundances_docker.smk'
#include: 'workflow/rules/metagenomic/refiltered_abundances.smk'
include: 'workflow/rules/metagenomic/refiltered_abundances_docker.smk'
#include: 'workflow/rules/metagenomic/bowtie.smk'
include: 'workflow/rules/metagenomic/bowtie_docker.smk'
