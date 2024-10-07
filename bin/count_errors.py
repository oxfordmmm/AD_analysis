#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def count_Ns(fasta_file):
    records = SeqIO.parse(fasta_file, 'fasta')
    rows = []
    for record in records:
        row = {}
        row['header'] = record.id
        sequence = str(record.seq)
        total_bases = len(sequence)
        n_count = sequence.count('N')

        row['total_bases'] = total_bases
        row['n_count'] = n_count
        rows.append(row)
    
    return pd.DataFrame(rows)

def count_differences(vcf_file):
    header = []
    data = []

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header = line.strip('#').strip().split('\t')
            else:
                data.append(line.strip().split('\t'))
    
    df = pd.DataFrame(data, columns=header)

    df = df.groupby('CHROM', as_index=False) \
            .agg(
                ref_diffs=('POS', 'count')
            ) \
            .rename(columns={'CHROM': 'header'}) \
            .reset_index(drop=True) \
            .copy()
    
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("consensus", help="fasta file of consensus")
    parser.add_argument("vcf", help="vcf file of differences from ref")
    parser.add_argument("outfile_csv", help="outfile with N count")
    args = parser.parse_args()

    ns_df = count_Ns(args.consensus)
    diff_df = count_differences(args.vcf)

    df = ns_df.merge(diff_df, on='header', how='left').fillna(0)

    df['header'] = df['header'].str.split('_').str[0]
    df = df.groupby('header').sum().reset_index()
    df['n_percent'] = round(100 * df['n_count'] / df['total_bases'], 2)
    df['diff_percent'] = round(100 * df['ref_diffs'] / df['total_bases'], 2)
    df.to_csv(args.outfile_csv, index=False)