#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

positive_controls=['MS2',
'murine_respirovirus',
'orthoreovirus_1',
'orthoreovirus_2',
'orthoreovirus_3',
'orthoreovirus_4',
'orthoreovirus_5',
'orthoreovirus_6',
'orthoreovirus_7',
'orthoreovirus_8',
'orthoreovirus_9',
'orthoreovirus_10',
'zika']

def getDataFrame(input):
    dfs=[]
    for folder in input:
        run=folder.split('/')[-1]
        df=pd.read_csv(f'{folder}/depth_summary/{run}_depth_summary.csv')
        dfs.append(df)
    df=pd.concat(dfs)
    # remove SQK-RBK114-96_barcode from Sample name
    df['Sample name']=df['Sample name'].str.replace('SQK-RBK114-96_barcode','')
    # remove unclassified samples
    df=df[df['Sample name']!='unclassified']
    # convert Sample name to int
    df['barcode']=df['Sample name'].map(int)
    # remove sup from batch
    df['run']=df['batch'].str.replace('_sup','')
    return df

def checkThresholds(df, percentage_type_reads, percentage_run_reads ):
    '''Calculate the thresholds for the run and filter'''
    cols=['run', 'barcode', 'chrom','sample num reads','batch']
    df=df[cols]
    total_reads=df['sample num reads'].sum()
        
    g=df.groupby(['chrom'])[['sample num reads']].sum()
    g['type numreads threshold']=(g['sample num reads']/100)*float(percentage_type_reads)
    g['run numreads threshold']=(total_reads/100)*float(percentage_run_reads)
    g.rename(columns={'sample num reads':'run num reads'},inplace=True)

    df=pd.merge(df,g, on='chrom',how='left')

    ## Apply thresholds
    df['type pass']=np.where(df['sample num reads']>=df['type numreads threshold'],True,False)
    df['run pass']=np.where(df['sample num reads']>=df['run numreads threshold'],True,False)
    df['pass']=np.where(df['type pass'] & df['run pass'],True,False)
    df=df[df['pass']==True]
    df['percentage_type_reads']=percentage_type_reads
    df['percentage_run_reads']=percentage_run_reads
    return df

def getMeta(meta, pathogens):
    df=pd.read_csv(meta)
    # strip whitespace from values
    df['pathogen 1'] = df['pathogen 1'].str.strip()

    df2=pd.read_csv(pathogens)
    
    path_dict={}
    for i,r in df2.iterrows():
        path_dict.setdefault(r['pathogen'],[]).append(r['reference'])
    
    metaDict={}
    for i,r in df.iterrows():
        runBar=r['Run']+'_'+str(r['barcode'])
        if r['pathogen 1'] in path_dict:
            metaDict.setdefault(runBar,[]).extend(path_dict[r['pathogen 1']])
        if r['pathogen 2'] in path_dict:
            metaDict[runBar].extend(path_dict[r['pathogen 2']])
        if r['pathogen 3'] in path_dict:
            metaDict[runBar].extend(path_dict[r['pathogen 3']])
    return metaDict, df 


def checkSensitivity(df, metaDict):
    barcodes=df['barcode'].unique()
    dfs=[]
    for barcode in barcodes:
        dfB=df[df['barcode']==barcode]
        #dfB.copy()
        dfB['runBar']=dfB['run']+'_'+dfB['barcode'].astype(str)
        runBar=dfB['runBar'].unique()[0]
        if runBar not in metaDict:
            continue
        species=metaDict[runBar]
        species_poscontrol=species + positive_controls
        print(species)
        dfB['TP']=np.where(dfB['chrom'].isin(species),True,False)
        dfB['pcTP']=np.where(dfB['chrom'].isin(positive_controls),True,False)
        dfB['FP']=np.where(~dfB['chrom'].isin(species_poscontrol),True,False)
        print(dfB)
        dfs.append(dfB)
    if len(dfs)>0:
        df=pd.concat(dfs)
    return df


def checkSensSpec(data, metaDict, metaDF):
    runs=metaDF['Run'].unique()
    dfs=[]
    for run in runs:
        df=data[data['run']==run]
        df2=checkThresholds(df, 1.0, 0.5)
        
        # check sensitivity
        df2=checkSensitivity(df2, metaDict)
        dfs.append(df2)
    
    data=pd.concat(dfs)
    return data


def main(args):
    df=getDataFrame(args.input)
    metaDict, metaDF=getMeta(args.meta, args.pathogens)
    print(metaDict)
    df=checkSensSpec(df, metaDict, metaDF)
    df.to_csv(args.output, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='Input folders')
    parser.add_argument('-m', '--meta', required=True,
                        help='Metadata file')
    parser.add_argument('-p', '--pathogens', required=True,
                        help='Pathogens file')
    parser.add_argument('-o', '--output', required=True,
                         help='Output file')
    args = parser.parse_args()

    main(args)

    