#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

def calculateThresholds(df, percentage_type_reads, percentage_run_reads):
    total_reads=df['counts'].sum()
        
    g=df.groupby(['chrom'])[['counts']].sum()
    g['type numreads threshold']=(g['counts']/100)*float(percentage_type_reads)
    g['run numreads threshold']=(total_reads/100)*float(percentage_run_reads)
    g.rename(columns={'counts':'run num reads'},inplace=True)
    return g

def getSamples(files):
    dfs=[]
    for f in files:
        df=pd.read_csv(f)
        dfs.append(df)
    df=pd.concat(dfs)
    return df

def applyThresholds(df):
    df['type pass']=np.where(df['sample num reads']>=df['type numreads threshold'],True,False)
    df['run pass']=np.where(df['sample num reads']>=df['run numreads threshold'],True,False)
    df['pass']=np.where(df['type pass'] & df['run pass'],True,False)
    return df

def run(args):
    df=getSamples(args.input)
    thresholds=calculateThresholds(df, args.percentage_type_reads, args.percentage_run_reads)

    df=pd.merge(df,thresholds,on='chrom',how='left')
    df.rename(columns={'counts':'sample num reads'},inplace=True)

    df=applyThresholds(df)
    df.to_csv(args.output,index=False)
    df2=df[df['pass']==True]
    df2.to_csv(args.output.replace('.csv','_pass.csv'),index=False)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='combine reports from all samples for overall run report and stats')

    argparser.add_argument('-i', '--input', help='Input files', required=True, nargs='+')
    argparser.add_argument('-o', '--output', help='Output file', required=True)
    argparser.add_argument('-p', '--percentage_type_reads', required=False, default=1.0,
                        help='Percentage of total segment reads required to pass, default=1%')
    argparser.add_argument('-r', '--percentage_run_reads', required=False, default=0.5,
                        help='Percentage of total run reads required to pass, default=0.5%')

    args=argparser.parse_args()

    run(args)
    