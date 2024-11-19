#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def plot_bacteria(df, out_file):
    bacterial_refs=['CP085971.1','NZ_CP025371.1','NC_005043.1','NZ_LR214945.1']
    df2=df[df['chrom'].isin(bacterial_refs)]
    df3=df2.groupby(['chrom']).sample(n=1000)

    ax=sns.FacetGrid(df3, col='chrom', col_wrap=1, height=5)
    ax.map(sns.lineplot, 'position', 'depth' )
    plt.savefig(f'{out_file}.pdf')

    df2=df2[df2['depth']>0]
    df2.to_csv(f'{out_file}.csv',index=False)

if __name__ == '__main__':
    df=pd.read_csv(sys.argv[1], sep='\t', names=['chrom','position','depth'])
    plot_bacteria(df, sys.argv[2])