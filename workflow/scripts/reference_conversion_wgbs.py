#!/usr/bin/env python
import sys

def parse_arguments():
    input_fasta = sys.argv[1]
    CT = sys.argv[2]
    GA = sys.argv[3]
    return(input_fasta, CT, GA)

def convert_database(input_fasta, CT, GA):
    with open(input_fasta) as f, open(CT, 'w') as ct, open(GA, 'w') as ga:
        for line in f:
            if line[0] == '>':
                ct.write(line)
                ga.write(line)
            else:
                ct.write(line.upper().replace('C', 'T'))
                ga.write(line.upper().replace('G', 'A'))

def main():
    input_fasta, CT, GA = parse_arguments()
    convert_database(input_fasta, CT, GA)

if __name__ == '__main__':
    main()
