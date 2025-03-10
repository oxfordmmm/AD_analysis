#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def plot_bacteria(df2, out_file):
    #bacterial_refs=['CP085971.1','NZ_CP025371.1','NC_005043.1','NZ_LR214945.1']
    #df2=df[df['chrom'].isin(bacterial_refs)]
    mins=df2.groupby('chrom')['position'].min()
    maxs=df2.groupby('chrom')['position'].max()

    maxLen=df2['position'].max()

    df2=df2[df2['depth']>0]
    print(df2)

    df3=df2#.groupby(['chrom']).sample(n=1000)

    if len(df3)==0:
        return

    ax=sns.FacetGrid(df3, col='chrom', col_wrap=1, height=5)
    ax.map(sns.lineplot, 'position', 'depth' )
    ax.set(xlim=(0,maxLen))
    plt.savefig(f'{out_file}.pdf')

    df2=df2[df2['depth']>0]
    df2.to_csv(f'{out_file}.csv',index=False)

if __name__ == '__main__':
    df=pd.read_csv(sys.argv[1], sep='\t', names=['chrom','position','depth'])
    plot_bacteria(df, sys.argv[2])