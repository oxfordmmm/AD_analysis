#!/usr/bin/env python3
import pandas as pd
import sys

## Analysis of workflow results for derivation set using the HAC basecall model


def run():
    df=pd.read_csv(sys.argv[1])

    # combine Run and barcode to make unique sample ID
    df['barcode']=df['barcode'].astype(str)
    df['Run']=df['Run'].astype(str)
    df['sample_id']=df['Run']+'_'+df['barcode']

    # remove umapped ref from df
    df2=df[df['pathogen']!='unmapped']
    # remove negative controls
    df2=df2[df2['test'].isin(['BIOFIRE', 'ALINITY', 'ALINTY'])]

    # count number of samples
    nsamples=df2['sample_id'].nunique()
    print(f"Number of samples y: {nsamples}")

    # remove gold standard and count number of samples with reads mapping to pathogen
    df2=df2[df2['gold_standard']!=1]
    df2=df2.copy()
    df2=df2[df2['sample num reads']>0]
    nsamples_False_reads=df2['sample_id'].nunique()
    print(f"Number of samples with reads mapping to pathogen X: {nsamples_False_reads}")
    print(f'X/y: {nsamples_False_reads}/{nsamples} ({100*(nsamples_False_reads/nsamples):.0f}%)')

    # A/B (%) true negative samples
    df3=df[df['pathogen']!='unmapped']
    df3=df3[~df3['test'].isin(['BIOFIRE', 'ALINITY'])]
    nsamples_true_neg=df3['sample_id'].nunique()
    print(f"Number of true negative samples B: {nsamples_true_neg}")
    df3=df3[df3['sample num reads']>0]
    nsamples_true_neg_reads=df3['sample_id'].nunique()
    print(f"Number of true negative samples with reads mapping to pathogen A: {nsamples_true_neg_reads}")
    print(f'A/B: {nsamples_true_neg_reads}/{nsamples_true_neg} ({100*(nsamples_true_neg_reads/nsamples_true_neg):.0f}%)')

    # how many actually mapped to the pathogen identified by gold standard
    df4=df[df['gold_standard']==1]
    nsamples_gold_standard=df4['sample_id'].nunique()
    print(f"Number of samples with gold standard pathogen X: {nsamples_gold_standard}")
    df4=df4[df4['sample num reads']>0]
    df4.to_csv('gold_standard_reads.csv', index=False)
    nsamples_gold_standard_reads=df4['sample_id'].nunique()
    print(f"Number of samples with reads mapping to gold standard pathogen Z: {nsamples_gold_standard_reads}")
    print(f'Z/X: {nsamples_gold_standard_reads}/{nsamples_gold_standard} ({100*(nsamples_gold_standard_reads/nsamples_gold_standard):.0f}%)')

    # calculate median read length and IQR
    print("Median read length and IQR for reads mapping to pathogen")
    df5=df[df['pathogen']!='unmapped']
    df5=df5[df5['sample num reads']>0]
    median_read_length=df5['median_read_length'].median()
    q1=df5['median_read_length'].quantile(0.25)
    q3=df5['median_read_length'].quantile(0.75)
    iqr=q3-q1
    print(f"median_read_length: {median_read_length:.0f} ({q1:.0f}-{q3:.0f})")
   # print(f"IQR: {iqr}")

if __name__ == "__main__":
    run()