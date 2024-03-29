#!/usr/bin/python3
"""
Gets taxid and filters based on logical strands
File must come as stdin
"""
import sys
import subprocess
import argparse

def argument_parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o',"--filename_out", type=str,help="Name of output file")
    parser.add_argument("--acc_to_tax", type=str,help="acc to taxid file")
    parser.add_argument("--conversion", type=str, help="Conversion", choices= ['CT','GA','std'])
    parser.add_argument("--human_like", type=str, help="Human-like file")
    parser.add_argument("--rejected_hits", type=str, help="Rejected hits")

    args = parser.parse_args()


    outfile = args.filename_out
    acc_to_tax = args.acc_to_tax
    conversion = args.conversion
    human_like = args.human_like
    rejected = args.rejected_hits

    return(outfile, acc_to_tax, conversion, human_like, rejected)

def get_gi(sseqid):
    gi = sseqid.split('|')[1]
    return(gi)

def get_taxid(gi, acc_to_tax):
    if gi not in acc_to_tax:
        #print(gi)
        taxid = 'GI_NOT_FOUND'
    else:
        taxid = acc_to_tax[gi]
    return(taxid)


def main():
    outfile, acc_to_tax, conversion, human_like, rejected = argument_parsing()
    with open(acc_to_tax) as f:
        acc_to_tax = dict(x.rstrip().split('\t') for x in f)

    with open(outfile, 'w') as w1, open(human_like,'w') as w2, open(rejected, 'w') as w3:
        for line in sys.stdin:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart,qend, sstart, send, evalue, bitscore, qlen, strand = line.strip().split('\t')
            gi = sseqid
            taxid = get_taxid(gi, acc_to_tax)

            if taxid == "9606":
                w2.write('\t'.join([qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, strand, taxid]) + '\n')
            else:
                if conversion == 'CT' and strand == 'Plus/Plus':
                    w1.write('\t'.join([qseqid, sseqid, pident, length, mismatch,  gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, strand, taxid]) + '\n')
                elif conversion == 'GA' and strand == 'Plus/Minus':
                    w1.write('\t'.join([qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, strand, taxid]) + '\n')
                elif conversion == 'std':
                    w1.write('\t'.join([qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, strand, taxid]) + '\n')
                else:
                    w3.write('\t'.join([qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, strand, taxid]) + '\n')

if __name__ == '__main__':
    main()
