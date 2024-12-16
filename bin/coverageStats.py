#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def getData(f,f2, multiSegmentDict):
    names=['original_chrom','position','depth']
    df=pd.read_csv(f,sep='\t',names=names)
    
    df['chrom']=df['original_chrom'].map(multiSegmentDict)
    df['chrom']=df['chrom'].fillna(df['original_chrom'])

    df3=pd.read_csv(f2) 
    
    df3['chrom2']=df3['chrom'].map(multiSegmentDict)
    df3['chrom2']=df3['chrom2'].fillna(df3['chrom'])
    df3['chrom']=df3['chrom2']
    df3.drop(columns=['chrom2'],inplace=True)

    df3['non-unique counts']=df3.groupby(['Sample name','chrom'])['non-unique counts'].transform('sum')
    df3['counts']=df3.groupby(['Sample name','chrom'])['counts'].transform('sum')
    df3.drop_duplicates(subset=['Sample name','chrom'],keep='first', inplace=True)
    #df3.drop(columns=['Sample name'],inplace=True)
    #df=df.merge(df3,on=['chrom'], how='outer')
    
    return df, df3

def coverageStats(df, df3, df4):
    cov1=df[df.depth >= 1].groupby(['chrom']).count()
    cov3=df[df.depth >= 3].groupby(['chrom']).count()
    cov5=df[df.depth >= 5].groupby(['chrom']).count()
    cov10=df[df.depth >= 10].groupby(['chrom']).count()

    bases=df.groupby(['chrom'])['depth'].sum()
    chromLens=df.groupby(['chrom'])['position'].count()
    meanDepth=bases/cov1['depth']

    # Additional columns
    df['trunc5']=np.where(df['depth']>5,5,df['depth'])
    df['trunc10']=np.where(df['depth']>10,10,df['depth'])
    AuG_trunc5=df.groupby(['chrom'])['trunc5'].sum()
    AuG_trunc10=df.groupby(['chrom'])['trunc10'].sum()
    meanDepth_trunc5=AuG_trunc5/cov1['depth']
    meanDepth_trunc10=AuG_trunc10/cov1['depth']

    cov1.reset_index(inplace=True)
    cov3.reset_index(inplace=True)
    cov3.rename(columns={'depth':'depth cov3','position':'position cov3'},inplace=True)
    cov5.reset_index(inplace=True)
    cov5.rename(columns={'depth':'depth cov5','position': 'position cov5' },inplace=True)
    cov10.reset_index(inplace=True)
    print(cov3)

    df2=cov1.merge(cov10,on='chrom',suffixes=[' cov1',' cov10'],how='outer')
    df2=df2.merge(cov3,on='chrom',how='outer')
    df2=df2.merge(cov5,on='chrom',how='outer')
    print(df2.columns)
    df2['length']=df2.chrom.map(chromLens)
    df2['covBreadth1x']=(df2['depth cov1']/df2['length'])*100
    df2['covBreadth3x']=(df2['depth cov3']/df2['length'])*100
    df2['covBreadth5x']=(df2['depth cov5']/df2['length'])*100
    df2['covBreadth10x']=(df2['depth cov10']/df2['length'])*100
    df2['meanDepth']=df2.chrom.map(meanDepth)
    df2['AuG_trunc5_sum']=df2.chrom.map(AuG_trunc5)
    df2['AuG_trunc10_sum']=df2.chrom.map(AuG_trunc10)
    df2['AuG_trunc5']=df2['AuG_trunc5_sum']/df2['length']
    df2['AuG_trunc10']=df2['AuG_trunc10_sum']/df2['length']
    df2['meanDepth_trunc5']=df2.chrom.map(meanDepth_trunc5)
    df2['meanDepth_trunc10']=df2.chrom.map(meanDepth_trunc10)
    df2['Sample name']=sys.argv[2]
    df2['bases']=df2.chrom.map(bases)
    df2['AuG']=df2['bases']/df2['length']
    df2['Bases_perc']=df2['AuG']*100

    df2=df2.merge(df3,on=['Sample name', 'chrom'],how='outer')

    df2=df2[['Sample name','chrom','length',
             'counts', 'non-unique counts', 'bases','meanDepth',
             'position cov1', 'position cov3','position cov5', 'position cov10',
             'covBreadth1x', 'covBreadth3x','covBreadth5x','covBreadth10x',
             'meanDepth_trunc5','meanDepth_trunc10','AuG','AuG_trunc5','AuG_trunc10']]
    

    if len(df4)==0:
        print('no alignment info')
        df2['median_read_length']=None
        df2['median_aligned_length']=None
        df2['mean_read_length']=None
        df2['mean_aligned_length']=None
        for c in [200, 300, 400]:
            df2[f'Sample_num_reads_{c}']=None
        df2.to_csv('coverage_stats.csv',index=False)
        return

    median_read_length=df4.groupby('ref')['read_len'].median()
    mean_read_length=df4.groupby('ref')['read_len'].mean()
    median_aligned_length=df4.groupby('ref')['align_len'].median()
    mean_aligned_length=df4.groupby('ref')['align_len'].mean()
    print(mean_aligned_length)
    df2['median_read_length']=df2['chrom'].map(median_read_length)
    df2['median_aligned_length']=df2['chrom'].map(median_aligned_length)
    df2['mean_read_length']=df2['chrom'].map(mean_read_length)
    df2['mean_aligned_length']=df2['chrom'].map(mean_aligned_length)

    for c in [200, 300, 400]:
        df5=df4[df4['align_len']>c]
        sample_num_reads=df5.groupby('ref')['query'].count()
        df2[f'Sample_num_reads_{c}']=df2['chrom'].map(sample_num_reads)


    df2.to_csv('coverage_stats.csv',index=False)
    #print(df2)
    return
    


multiSegment=pd.read_csv(sys.argv[4])
path_dict={}
for i,r in multiSegment.iterrows():
    path_dict.setdefault(r['pathogen'],[]).append(r['reference'])
# reverse path_dict

multiSegmentDict={v:k for k,values in path_dict.items() for v in values}
df, df3=getData(sys.argv[1], sys.argv[3], multiSegmentDict)

try:
    print('reading aligntment info')
    df4=pd.read_csv(sys.argv[5])
    df4['chrom2']=df4['ref'].map(multiSegmentDict)
    df4['chrom2']=df4['chrom2'].fillna(df4['ref'])
    df4['ref']=df4['chrom2']
except:
    print('no alignment info')
    df4=pd.DataFrame()

coverageStats(df, df3, df4)
