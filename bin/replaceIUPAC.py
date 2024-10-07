#!/usr/bin/env python3
from Bio import SeqIO
import sys
from argparse import ArgumentParser

def _getSeqs(fa):
    for seq in SeqIO.parse(open(fa,'rt'), 'fasta'):
        for b in ['R', 'Y', 'S', 'W', 'K', 'M', 'B','M','D','H', 'V']:
            seq.seq=seq.seq.replace(b,'N')
        yield seq


def run(opts):
    outSeqs=_getSeqs(opts.input_fasta)
    outSeqs=list(outSeqs)
    with open(opts.output_fasta, 'wt') as outf:
        SeqIO.write(outSeqs, outf, 'fasta')


if __name__ == '__main__':
    parser = ArgumentParser(description='replace IUPAC codes with N from fasta seqs')
    parser.add_argument('-i', '--input_fasta', required=True, 
	    help='input fasta file')
    parser.add_argument('-o', '--output_fasta', required=True, 
	    help='output fasta file')

opts, unknown_args = parser.parse_known_args()

run(opts)