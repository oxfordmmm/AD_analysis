#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

def calculateThresholds(df, percentage_type_reads, percentage_run_reads):

    total_reads_inc_unmapped=df['counts'].sum()
    df=df[df['chrom']!='unmapped']
    total_reads=df['counts'].sum()
        
    g=df.groupby(['chrom'])[['counts']].sum()
    g['type numreads threshold']=(g['counts']/100)*float(percentage_type_reads)
    g['run numreads threshold']=(total_reads/100)*float(percentage_run_reads)
    g['total run reads mapped']=total_reads
    g['total run reads inc unmapped']=total_reads_inc_unmapped
    g.rename(columns={'counts':'run num reads'},inplace=True)
    return g

def getSamples(files):
    dfs=[]
    for f in files:
        df=pd.read_csv(f)
        dfs.append(df)
    df=pd.concat(dfs)
    df['bases_perc']=(df['bases']/df['length'])*100
    return df

def getMeta(args):
    df=pd.read_csv(args.meta_pathogens)
    
    path_dict={}
    for i,r in df.iterrows():
        path_dict.setdefault(r['pathogen'],[]).append(r['reference'])

    # reverse path_dict
    path_dict_rev={v:k for k,values in path_dict.items() for v in values}

    df2=pd.read_csv(args.pathogen_reduced)
    df2.set_index('pathogen',inplace=True)
    path_dict_reduced=df2.to_dict()['pathogen_reduced']

    return path_dict_rev, path_dict_reduced

def applyThresholds(df):
    df['Sample_reads_percent_of_run']=(df['sample num reads']/df['total run reads inc unmapped'])*100
    df['Sample_reads_percent_of_refs']=(df['sample num reads']/df['total run reads mapped'])*100
    df['pathogen reads sample']=df.groupby(['Sample name','pathogen_reduced'])['sample num reads'].transform('sum')
    df['pathogen reads run']=df.groupby(['pathogen_reduced'])['sample num reads'].transform('sum')
    df['Sample_reads_percent_of_type_sample']=(df['sample num reads']/df['pathogen reads sample'])*100
    df['Sample_reads_percent_of_type_run']=(df['sample num reads']/df['pathogen reads run'])*100
    df['type pass']=np.where(df['sample num reads']>=df['type numreads threshold'],True,False)
    df['run pass']=np.where(df['sample num reads']>=df['run numreads threshold'],True,False)
    df['pass']=np.where(df['type pass'] & df['run pass'],True,False)
    return df

def run(args):
    df=getSamples(args.input)
    ref_map, pathogen_reduced = getMeta(args)
    df['pathogen']=df['chrom'].map(ref_map)
    df['pathogen_reduced']=df['pathogen'].map(pathogen_reduced)

    thresholds=calculateThresholds(df, args.percentage_type_reads, args.percentage_run_reads)

    df=pd.merge(df,thresholds,on='chrom',how='left')
    df.rename(columns={'counts':'sample num reads'},inplace=True)

    df=applyThresholds(df)
    df.sort_values(by=['Sample name','chrom'], ascending=[True,True], inplace=True)
    df['batch']=args.batch
    rename_dict={'position cov1':'Cov1','position cov3':'Cov3','position cov5':'Cov5','position cov10':'Cov10',
                 'covBreadth1x': 'Cov1_perc', 'covBreadth3x':'Cov3_perc','covBreadth5x':'Cov5_perc','covBreadth10x':'Cov10_perc'}
    df.rename(columns=rename_dict,inplace=True)
    cols=['batch','Sample name','chrom', 'length','pathogen', 'pathogen_reduced',
          'meanDepth', 'meanDepth_trunc5', 'meanDepth_trunc10','AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc', 
          'Cov1', 'Cov3', 'Cov5', 'Cov10', 
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc',
          'sample num reads','total run reads mapped', 'total run reads inc unmapped', 
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample']
    #df.to_csv('test.csv')
    df=df[cols]
    df.to_csv(args.output,index=False)
    #df2=df[df['pass']==True]
    #df2.sort_values(by=['Sample name','chrom'], ascending=[True,True], inplace=True)
    #df2.to_csv(args.output.replace('.csv','_pass.csv'),index=False)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='combine reports from all samples for overall run report and stats')

    argparser.add_argument('-i', '--input', help='Input files', required=True, nargs='+')
    argparser.add_argument('-o', '--output', help='Output file', required=True)
    argparser.add_argument('-b', '--batch', required=False, default=1,
                        help='Batch name or run name')
    argparser.add_argument('-p', '--percentage_type_reads', required=False, default=1.0,
                        help='Percentage of total segment reads required to pass, default=1%')
    argparser.add_argument('-r', '--percentage_run_reads', required=False, default=0.5,
                        help='Percentage of total run reads required to pass, default=0.5%')
    argparser.add_argument('-mp', '--meta_pathogens', required=True, default=None,
                        help='Path to meta pathogens file containing pathgen:chromosome mapping')
    argparser.add_argument('-pr', '--pathogen_reduced', required=True, default=None,
                           help='Path to reduced pathogens file containing pathgen:reduced(biofire) mapping')

    args=argparser.parse_args()

    run(args)
    