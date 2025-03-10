#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
from argparse import ArgumentParser

pathogens=['Acanthamoeba' ]

index_cols=['run','barcode']
meta_cols=['sample_name', 'seq_name', 'pathogen 1', 'pathogen 2', 'pathogen 3']
names=['percent','clade_reads','taxon_reads','minizers','distinct minizers','taxid','Scientific Name']
#names=['percent','clade_reads','taxon_reads','taxid','Scientific Name']
def read_reports(l):
    dfs=[]
    for f in l:
        try:
            df=pd.read_csv(f,sep='\t', names=names)
            df['run']=f.split('/')[-3]
            df['barcode']=f.split('/')[-2]
            dfs.append(df)
        except:
            pass
    
    df=pd.concat(dfs)
    return df

def clean_df(df):
    df.reset_index(inplace=True)
    print(df.columns)
   
    for p in pathogens:
        if p not in df.columns:
            df[p]=None
    df=df[index_cols+pathogens]
    return df

def run(args):
    df=read_reports(args.input)
    print(df)
    cols=['run','barcode','Scientific Name','Clades']
    #df2=df[cols]
    # left strip the Scientific Name column
    df['Scientific Name']=df['Scientific Name'].str.lstrip()
    df=df[df['Scientific Name'].isin(pathogens)]
    df2=df.pivot(index=['run','barcode'], columns='Scientific Name', values='percent')
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
    parser.add_argument('-m','--meta', required=False,
                        help='metadata file')
    args = parser.parse_args()
    run(args)