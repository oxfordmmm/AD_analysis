#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os
import argparse 
import seaborn as sns
import matplotlib.pyplot as plt


def plot_blast_results(df):
    g = sns.FacetGrid(df, col="barcode", hue='Front/Rear', col_wrap=4, height=4, aspect=1)
    #g.map_dataframe(sns.boxenplot, x='query',y='length')
    g.map(plt.hist, 'query', bins=20)
    g.set_xticklabels(rotation=90)
    #g.set_xlabels('Barcode found from blast')
    #g.set_ylabels('Length of the blast hit')
    g.savefig('blast_results_hist.pdf',format='pdf',dpi=300)




def main(args):
    dfs=[]
    for barcode in args.input:
        df=pd.read_csv(barcode,sep='\t',names=['query','subject','identity','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore', 'qlen', 'slen'])
        df['barcode']=barcode.split('/')[-1].split('.')[0]
        df=df[df['length']>23]
        df['identity']=df['identity'].map(int)
        df=df[df['identity']>=100]
        df['Half']=df['slen']/2
        df['5prime/3prime']=np.where(df['sstart']<df['Half'],'5prime','3prime')
        dfs.append(df)
    df=pd.concat(dfs)
    df.to_csv('blast_results.tsv',sep='\t',index=False)

    plot_blast_results(df)

# main entry 
if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Plot the blast results of the barcode sequences')
    argparser.add_argument('-i', '--input', help='Input files', required=True, nargs='+')

    args = argparser.parse_args()
    main(args)