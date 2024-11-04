#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics

import seaborn as sns
import matplotlib.pyplot as plt

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

def getMeta(meta, pathogens, pathogens_reduced):
    df=pd.read_csv(meta)
    # strip whitespace from values
    df['pathogen 1'] = df['pathogen 1'].str.strip()

    df2=pd.read_csv(pathogens)
    
    path_dict={}
    for i,r in df2.iterrows():
        path_dict.setdefault(r['pathogen'],[]).append(r['reference'])

    # reverse path_dict
    path_dict_rev={v:k for k,values in path_dict.items() for v in values}
    
    metaDict={}
    for i,r in df.iterrows():
        runBar=r['Run']+'_'+str(r['barcode'])
        if r['pathogen 1'] in path_dict:
            metaDict.setdefault(runBar,[]).extend(path_dict[r['pathogen 1']])
        if r['pathogen 2'] in path_dict:
            metaDict[runBar].extend(path_dict[r['pathogen 2']])
        if r['pathogen 3'] in path_dict:
            metaDict[runBar].extend(path_dict[r['pathogen 3']])

    # read pathogenes reduced and create a dictionary with pathogen as key and pathogen_reduced as value
    df3=pd.read_csv(pathogens_reduced)
    df3.set_index('pathogen',inplace=True)
    path_dict_reduced=df3.to_dict()['pathogen_reduced']

    ## count number of pathogens in the test to caluculate TN
    num_pathogens=len(df3['pathogen_reduced'].unique())
    
    return metaDict, df, path_dict_rev, path_dict_reduced, num_pathogens

def checkFN(chroms, species, path_dict, path_dict_reduced):
    '''Check for false negatives. Convert chroms and species to same species and check for differences'''
    chroms=set(path_dict_reduced[path_dict[c.rstrip()].rstrip()] for c in chroms)
    species=set(path_dict_reduced[path_dict[s.rstrip()].rstrip()] for s in species)
    FN=0

    new_rows=[]
    for s in species:
        if s not in chroms:
            FN+=1
            d={'chrom':s,'sample num reads':0, 'pcTP':0, 'TP':0, 'FP':0, 'FN':1}
            new_rows.append(d)
    return FN, new_rows

def checkSensitivity(df, metaDict, path_dict, path_dict_reduced):
    '''Iterate through barcodes from the run and filter out any barcodes not in the metaDict,
    then calculate the true positives, false positives'''
    barcodes=df['barcode'].unique()
    dfs=[]
    for barcode in barcodes:
        dfB=df[df['barcode']==barcode].copy()
        dfB['runBar']=dfB['run']+'_'+dfB['barcode'].astype(str)
        runBar=dfB['runBar'].unique()[0]
        
        # only use barcodes in the metaDict
        if runBar not in metaDict:
            continue
        
        # Get a list of species and positive controls
        species=metaDict[runBar]
        species_poscontrol=species + positive_controls
        
        dfB['TP']=np.where(dfB['chrom'].isin(species),1,0)
        dfB['pcTP']=np.where(dfB['chrom'].isin(species_poscontrol),1,0)
        dfB['FP']=np.where(~dfB['chrom'].isin(species_poscontrol),1,0)
       
        # Look for false negatives
        chroms=dfB['chrom'].unique()
        dfB['FN'],new_rows =checkFN(chroms, species, path_dict, path_dict_reduced)
        if len(new_rows)>0:
            dfNR=pd.DataFrame(new_rows)
            dfNR['run']=dfB['run'].unique()[0]
            dfNR['barcode']=barcode
            dfNR['runBar']=runBar
            dfB=pd.concat([dfB,dfNR])
        dfs.append(dfB)

        # add in positive controls if not seen
        pc_missed=[]
        for pc in positive_controls:
            if pc not in chroms:
                d={'run':dfB['run'].unique()[0],
                   'barcode':barcode,
                   'runBar':runBar,
                   'chrom':pc,
                   'sample num reads':0, 
                   'pcTP':1, 'TP':0, 'FP':0, 'FN':0}
                pc_missed.append(d)
        if len(pc_missed)>0:
            dfPCmissed=pd.DataFrame(pc_missed)
            dfB=pd.concat([dfB,dfPCmissed])


    if len(dfs)>0:
        df=pd.concat(dfs)
    return df


def checkSensSpec(data, metaDict, metaDF, path_dict, path_dict_reduced, num_pathogens, percentage_type_reads=1.0, percentage_run_reads=0.5):
    runs=metaDF['Run'].unique()
    samples=metaDF.groupby(['Run'])['barcode'].nunique().reset_index()
    samples.rename(columns={'barcode':'samples','Run':'run'},inplace=True)
    dfs=[]
    for run in runs:
        df=data[data['run']==run]
        df2=checkThresholds(df, percentage_type_reads, percentage_run_reads)
        
        # check sensitivity
        df2=checkSensitivity(df2, metaDict, path_dict, path_dict_reduced)
        if len(df2)>0:
            dfs.append(df2)
    
    if len(dfs)>0:
        data=pd.concat(dfs)

        # Create a second dataframe of one sample per row with TP, FP, FN, TN
        #print(data)
        TP=data.groupby(['run','barcode'])[['TP']].sum().reset_index()
        pcTP=data.groupby(['run','barcode'])[['pcTP']].sum().reset_index()
        FP=data.groupby(['run','barcode'])[['FP']].sum().reset_index()
        FN=data.groupby(['run','barcode'])[['FN']].max().reset_index()
        
        #print(samples)

        df=pd.merge(TP,pcTP, on=['run','barcode'],how='left')
        df=pd.merge(df,FP, on=['run','barcode'],how='left')
        df=pd.merge(df,FN, on=['run','barcode'],how='left')
        
        

        # calculate the TPR and FPR rates
        df2=df.groupby('run')[['TP','pcTP','FP','FN']].sum().reset_index()
        df2['TPR']=df2['TP']/(df2['pcTP']+df2['FN'])
        df2['pcTPR']=df2['TP']/(df2['TP']+df2['FN'])
        df2=pd.merge(df2,samples, on=['run'],how='left')

        # TN either number of pathogens to look for minus FP,
        # Or number of samples minus FP
        df2['TN']=(num_pathogens*df2['samples'])-df2['FP']
        #df2['TN']=df2['samples']-df2['FP']
        df2['TN']=np.where(df2['TN']<0,0,df2['TN'])

        df2['FPR']=df2['FP']/(df2['FP']+df2['TN'])
        df2['TNR']=df2['TN']/(df2['FP']+df2['TN'])
        df2['1-specificity']=1-df2['TNR']
        df2['percentage_type_reads']=percentage_type_reads
        df2['percentage_run_reads']=percentage_run_reads
    
        return data, df, df2
    else:
        return None, None, None

def plotROC(df2):
    '''Plot the ROC curve'''
    #print(df)
    #chrom='MS2'
    #df2=df[df['chrom']==chrom]
    X=df2[['sample num reads']]
    y=df2['pcTP']
    #split the dataset into training (70%) and testing (30%) sets
    X_train,X_test,y_train,y_test = train_test_split(X,y,test_size=0.3,random_state=0) 

    #print(X_train)
    #print(y_train)
    #instantiate the model
    log_regression = LogisticRegression()

    #fit the model using the training data
    log_regression.fit(X_train,y_train)

    #define metrics
    y_pred_proba = log_regression.predict_proba(X)[::,1]
    fpr, tpr, _ = metrics.roc_curve(y,  y_pred_proba)
    auc = metrics.roc_auc_score(y, y_pred_proba)

    #create ROC curve
    plt.plot(fpr,tpr, label="AUC="+str(auc))
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.legend(loc=4)
    plt.savefig('roc_curve.png')
    
def main(args):
    df=getDataFrame(args.input)
    metaDict, metaDF, path_dict, pathogens_reduced, num_pathogens=getMeta(args.meta, args.pathogens, args.pathogens_reduced)
    #print(metaDict)
    dfs=[]
    thresholds=[0.0, 0.1, 0.3, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2.4, 3.3, 4.2, 5.1, 6.0, 7.0, 8.0, 9.0, 10.0, 20, 30, 40, 50, 60, 70, 80, 90, 99]
    thresholds=[0.0]
    for p in thresholds:
        for r in thresholds:
            print(p,r)
            data, df2, sensSpec=checkSensSpec(df, metaDict, metaDF, path_dict, pathogens_reduced, num_pathogens, percentage_type_reads=p, percentage_run_reads=r)
            data.to_csv(f'{args.output}_{p}_{r}.csv')
            if sensSpec is not None:
                dfs.append(sensSpec)
    sensSpec=pd.concat(dfs)
    sensSpec.to_csv(f'{args.output}_sensSpec.csv')
    plotROC(data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='Input folders')
    parser.add_argument('-m', '--meta', required=True,
                        help='Metadata file')
    parser.add_argument('-p', '--pathogens', required=True,
                        help='Pathogens file')
    parser.add_argument('-pr', '--pathogens_reduced', required=True,
                        help='Pathogens reduced file')
    parser.add_argument('-o', '--output', required=True,
                         help='Output file')
    args = parser.parse_args()

    main(args)
    #df=pd.read_csv(args.output)
    #plotROC(df)

    