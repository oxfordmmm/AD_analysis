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
    df['Sample_reads_percent_of_refs_AuG_truc10']=df['Sample_reads_percent_of_refs']/df['AuG_trunc10']

    
    # thresholds
    ## run sequencing pass if bases > 400 Mb
    ## sample sequencing pass if reads > 30,000. 
    ## pathogen detected if Sample_num_reads >=2 AND
    ## Aug_trunc10 > 0.003 OR
    ## Cov1_perc > 0.25 OR
    ## Sample_reads_percent_of_refs > 0.007
    ## AND sample_reads_percent_of_refs/aug_trunc10>0.1
    #df['pass']=np.where(df['bases']>400000000,True,False)
    #df['pass']=np.where(df['sample num reads']>30000,True,False)
    
    df['pass']=np.where((df['AuG_trunc10']>0.003) | (df['covBreadth1x']>0.25) | (df['Sample_reads_percent_of_refs']>0.007),True,False)
    
    df['pass']=np.where(df['Sample_reads_percent_of_refs_AuG_truc10']>0.1,df['pass'],False)
    df['pass']=np.where(df['sample num reads']>=2,df['pass'],False)
    
    return df

def getSpikes(df):
    spikes=['orthoreovirus', 'zika', 'MS2', 'murine_respirovirus']
    dfSpikes=df[df['chrom'].isin(spikes)]
    spike_reads=dfSpikes[['Sample name', 'pathogen_reduced', 'sample num reads']]
    df_sr_T=spike_reads.pivot(index=['Sample name'], columns='pathogen_reduced', values='sample num reads')
    print(df_sr_T)

    df=df.merge(df_sr_T, on=['Sample name'], how='left')
    spikes=['orthoreovirus', 'zika', 'MS2', 'murine_respirovirus']
    for spike in spikes:
        df[f'{spike} passed']=np.where(df[spike]>=2, 1, 0)
    df['PC_PASSES']=df['orthoreovirus passed']+df['MS2 passed']+df['zika passed']+df['murine_respirovirus passed']
    df['PCs_passed']=np.where(df['PC_PASSES']>1, 1, 0)
    df['isSpike']=np.where(df['chrom'].isin(spikes),1,0)
    return df

def getSeqkits(files):
    l=[]
    for f in files:
        df=pd.read_csv(f, sep='\t')
        l.append(df)
    df=pd.concat(l)
    print(df)

    total_bases=df['sum_len'].sum()
    total_reads=df['num_seqs'].sum()

    df['total run bases']=total_bases
    df['total run reads']=total_reads
    df['Sample name']=df['file']
    df.rename(columns={'sum_len':'total sample bases','num_seqs':'total sample reads'},inplace=True)

    cols=['Sample name','total sample bases','total sample reads','total run bases','total run reads']
    df=df[cols]
    print(df)

    return df

def makePathogenReport(df, output, pathogen_reduced):
    '''Pivot table by batch,Sample name, pathogen_reduced and pass'''
    #df2=df[df['isSpike']==0]
    df2=df[['batch','Sample name','pathogen_reduced','pass']]
    batchID=df2['batch'].unique()[0]
    for pathogen in pathogen_reduced.values():
        for sampleID in df2['Sample name'].unique():
            if pathogen not in df2[df2['Sample name']==sampleID]['pathogen_reduced'].values:
                df2=df2._append({'batch':batchID,'Sample name':sampleID,'pathogen_reduced':pathogen,'pass':False},ignore_index=True)
   
    df3=df2.pivot_table(index=['batch','Sample name'], columns='pathogen_reduced', values='pass', aggfunc='first')
    df3.reset_index(inplace=True)
    df3.fillna(False,inplace=True)
    
    #df3.dropna(subset=['Sample name'],inplace=True)
    df3.to_csv(output,index=False)

def run(args):
    df=getSamples(args.input)
    sample_dfs=getSeqkits(args.tsv_seqkits)
    df=pd.merge(df,sample_dfs,on='Sample name',how='left')

    ref_map, pathogen_reduced = getMeta(args)
    df['pathogen']=df['chrom'].map(ref_map)
    df['pathogen_reduced']=df['pathogen'].map(pathogen_reduced)

    thresholds=calculateThresholds(df, args.percentage_type_reads, args.percentage_run_reads)
    df=pd.merge(df,thresholds,on='chrom',how='left')

    df.rename(columns={'counts':'sample num reads'},inplace=True)

    # get spike results
    df=getSpikes(df)

    df=applyThresholds(df)
    
    df.sort_values(by=['Sample name','chrom'], ascending=[True,True], inplace=True)
    
    df['batch']=args.batch
    rename_dict={'position cov1':'Cov1','position cov3':'Cov3','position cov5':'Cov5','position cov10':'Cov10',
                 'covBreadth1x': 'Cov1_perc', 'covBreadth3x':'Cov3_perc','covBreadth5x':'Cov5_perc','covBreadth10x':'Cov10_perc'}
    df.rename(columns=rename_dict,inplace=True)
    cols=['batch','Sample name','chrom', 'length','pathogen', 'pathogen_reduced', 'mapQ',
          'meanDepth', 'meanDepth_trunc5', 'meanDepth_trunc10','AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc', 
          'Cov1', 'Cov3', 'Cov5', 'Cov10', 
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc',
          'sample num reads','total run reads mapped', 'total run reads inc unmapped', 
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'median_read_length', 'median_aligned_length', 'mean_read_length', 'mean_aligned_length',
          'Sample_num_reads_200', 'Sample_num_reads_300', 'Sample_num_reads_400',
          'orthoreovirus', 'zika', 'MS2', 'murine_respirovirus',
          'orthoreovirus passed', 'zika passed', 'MS2 passed', 'murine_respirovirus passed',
          'total sample bases','total sample reads','total run bases','total run reads',
          'isSpike',
          'PC_PASSES','PCs_passed','pass']
    #df.to_csv('test.csv')
    df=df[cols]
    df.to_csv(args.output,index=False)

    makePathogenReport(df, f'{args.batch}_pathogen_report.csv', pathogen_reduced)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='combine reports from all samples for overall run report and stats')

    argparser.add_argument('-i', '--input', help='Input files', required=True, nargs='+')
    argparser.add_argument('-t', '--tsv_seqkits', help='Input seqkit tsv files', required=True, nargs='+')
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
    