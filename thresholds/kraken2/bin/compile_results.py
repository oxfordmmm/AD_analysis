#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
from argparse import ArgumentParser

pathogens=['Mastadenovirus', 
    'Coronaviridae', 'Human coronavirus 229E', 'Human coronavirus HKU1', 'Human coronavirus NL63', 'Human coronavirus OC43',
    'Influenza A virus', 'Influenza A virus (A/California/07/2009(H1N1))', 'Influenza A virus (A/Puerto Rico/8/1934(H1N1))', 'Influenza A virus (A/New York/392/2004(H3N2))',
    'Influenza B virus', 'Influenza B virus (B/Lee/1940)',
    'Enterovirus', 
    'Human metapneumovirus',
    'Human respirovirus 1','Human orthorubulavirus 2','Human respirovirus 3','Human parainfluenza virus 4a',
    'Respiratory syncytial virus',
    'Severe acute respiratory syndrome coronavirus 2',
    'Bordetella parapertussis', 'Bordetella pertussis', 'Chlamydia pneumoniae', 'Mycoplasmoides pneumoniae',
    'Escherichia phage MS2', 'Mammalian orthoreovirus','Respirovirus muris', 'Zika virus' ]

index_cols=['run','barcode']
meta_cols=['sample_name', 'seq_name', 'pathogen 1', 'pathogen 2', 'pathogen 3']

def read_reports(l):
    dfs=[]
    for f in l:
        df=pd.read_csv(f,sep='\t')
        dfs.append(df)
    df=pd.concat(dfs)
    return df

def clean_df(df):
    df.reset_index(inplace=True)
    df['run']=df['run'].str.replace('.sup_barcodes_sup_fastq','')
    df['barcode']=df['barcode'].str.replace('SQK-RBK114-96_barcode','')
    df['barcode']=df['barcode'].str.replace('unclassified','0')
    df=df[df['run']!='basecalled_fastq']
    df['barcode']=df['barcode'].astype(int)

    
    for p in pathogens:
        if p not in df.columns:
            df[p]=None
    df=df[index_cols+pathogens]
    return df

def run(args):
    df=read_reports(args.input)
    cols=['run','barcode','Scientific Name','Clades']
    df2=df[cols]
    df2=df2.pivot(index=['run','barcode'],columns='Scientific Name',values='Clades')
    df2=clean_df(df2)

    if args.meta:
        meta=pd.read_csv(args.meta)
        meta['barcode']=meta['barcode'].astype(int)
        df2=df2.merge(meta, left_on=['run','barcode'], right_on=['Run','barcode'])
        df2=df2[index_cols+meta_cols+pathogens]
    df2.to_csv(args.output)

if __name__ == '__main__':
    parser = ArgumentParser(description='Combine kraken2 reports')
    parser.add_argument('-i','--input', nargs='+', help='Input kraken2 reports')
    parser.add_argument('-o','--output', help='Output combined kraken2 report')
    parser.add_argument('-m','--meta', 
                        help='metadata file')
    args = parser.parse_args()
    run(args)