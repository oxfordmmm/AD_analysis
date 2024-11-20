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

def run(args):
    df=get_kreport(args.input)
    df2=filter_report(df, args.taxids)
    df2['barcode']=args.barcode
    df2['run']=args.run
    df2.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    parser = ArgumentParser(description='Filter kraken2 report by taxids')
    parser.add_argument('-i','--input', help='Input kraken2 report')
    parser.add_argument('-o','--output', help='Output filtered kraken2 report')
    parser.add_argument('-t','--taxids', nargs='+', 
                        help='Taxids to filter')
    parser.add_argument('-b','--barcode', help='Barcode number')
    parser.add_argument('-r','--run', help='Run name')
    args = parser.parse_args()
    run(args)

