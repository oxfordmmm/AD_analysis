#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO


def run(ref, mask, output):
    df = pd.read_csv(mask)
    maskDict = {}
    for i, row in df.iterrows():
        maskDict.setdefault(row['chrom'], []).append({'start': row['start'], 'end': row['end']})
    print(maskDict)
    seqs=[]
    for seq in SeqIO.parse(ref, 'fasta'):
        if seq.id in maskDict:
            for mask in maskDict[seq.id]:
                print(f'Masking {seq.id}')
                start=mask['start']
                end=mask['end']
                seq.seq = seq.seq[:start] + 'N'*(end-start) + seq.seq[end:]

                
        seqs.append(seq)

    SeqIO.write(seqs, output, 'fasta')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mask reference sequence')
    parser.add_argument('-r', '--reference', required=True, 
                        help='Reference sequence fasta file')
    parser.add_argument('-m', '--mask', required=True,
                        help='Masked regions csv file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output masked reference fasta file')
    args = parser.parse_args()

    run(args.reference, args.mask, args.output)

