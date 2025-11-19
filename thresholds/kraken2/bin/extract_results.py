#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
from argparse import ArgumentParser

def get_kreport(f):
    df=pd.read_csv(f,sep='\t')
    return df

def filter_report(df, taxids):
    taxids=[int(x) for x in taxids]
    df['Taxonomy ID']=df['Taxonomy ID'].astype(int)
    df2=df[df['Taxonomy ID'].isin(taxids)]
    df2['Scientific Name']=df2['Scientific Name'].str.lstrip()
    return df2

def getSeqkitstats(f):
    df=pd.read_csv(f, sep='\t')
    df['run']=df['file'].str.split('_barcode').str[0]
    df['barcode']=df['file'].str.split('_').str[-1]
    #df['barcode']=df['barcode'].str.replace('barcode','')
    df['barcode']=df['barcode'].str.replace('.gz','')
    #df['barcode']=df['barcode'].str.replace('.fastq','')
    return df

def run(args):
    df=get_kreport(args.input)
    df2=filter_report(df, args.taxids)
    df2['barcode']=args.barcode
    df2['run']=args.run
    df3=getSeqkitstats(args.seqkit_stats)
    df2=df2.merge(df3, on=['run','barcode'], how='left')
    df2.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    parser = ArgumentParser(description='Filter kraken2 report by taxids')
    parser.add_argument('-i','--input', help='Input kraken2 report')
    parser.add_argument('-o','--output', help='Output filtered kraken2 report')
    parser.add_argument('-t','--taxids', nargs='+', 
                        help='Taxids to filter')
    parser.add_argument('-b','--barcode', help='Barcode number')
    parser.add_argument('-r','--run', help='Run name')
    parser.add_argument('-s','--seqkit_stats', help='Seqkit stats file')
    args = parser.parse_args()
    run(args)

