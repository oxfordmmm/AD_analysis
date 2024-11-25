#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
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
'orthoreovirus',
'zika']

def reduce_fluA(df):
    df2=df.copy()
    df2=df[df['pathogen_reduced'].isin(['Influenza A/H1-2009','Influenza A/H3','Influenza A/H1'])]
    df2.sort_values(by=['sample num reads','Cov1_perc'],ascending=False, inplace=True)
    df2.drop_duplicates(subset=['run','barcode'],inplace=True, keep='first')
    df2['pathogen']='Influenza A'
    df2['pathogen_reduced']='Influenza A'
    #df=df[df['pathogen_reduced']!='Influenza A']
    # assign df2 to df by 'Run','barcode'
    df=pd.concat([df,df2])
    return df

def getDataFrame(input):
    dfs=[]
    for folder in input:
        run=folder.split('/')[-1]
        df=pd.read_csv(f'{folder}/depth_summary/{run}_depth_summary.csv')
        dfs.append(df)
    df=pd.concat(dfs)
    # remove SQK-RBK114-96_barcode from Sample name
    df['Sample name']=df['Sample name'].str.replace('SQK-RBK114-96_barcode','')
    # remove unclassified samples and nan
    df=df[df['Sample name']!='unclassified']
    df=df[~df['Sample name'].isnull()]
    # convert Sample name to int
    print(df['Sample name'].unique())
    df['barcode']=df['Sample name'].map(int)
    # remove sup from batch
    df['run']=df['batch'].str.replace('_sup','')
    # remove SARS_coronavirus_Tor2
    df=df[df['chrom']!='SARS_coronavirus_Tor2']

    # reduce fluA to one species 
    df=reduce_fluA(df)
    return df

def checkThresholds(df, percentage_type_reads, percentage_run_reads ):
    '''Calculate the thresholds for the run and filter'''
    #cols=['run', 'barcode', 'chrom','sample num reads','batch','length','bases', 'position cov1', 'position cov10','avDepth', 'covBreadth1x']
    #df=df[cols]
    total_reads=df['sample num reads'].sum()
        
    g=df.groupby(['chrom'])[['sample num reads']].sum()
    g['type numreads threshold']=(g['sample num reads']/100)*float(percentage_type_reads)
    g['run numreads threshold']=(total_reads/100)*float(percentage_run_reads)
    g.rename(columns={'sample num reads':'type num reads'},inplace=True)
    g['run num reads']=total_reads

    df=pd.merge(df,g, on='chrom',how='left')

    ## Apply thresholds
    df['sample reads percent of type']=(df['sample num reads']/df['type num reads'])*100
    df['sample reads percent of run']=(df['sample num reads']/df['run num reads'])*100
    df['type pass']=np.where(df['sample num reads']>=df['type numreads threshold'],True,False)
    df['run pass']=np.where(df['sample num reads']>=df['run numreads threshold'],True,False)
    df['pass']=np.where(df['type pass'] & df['run pass'],True,False)
    df=df[df['pass']==True]
    df['percentage_type_reads']=percentage_type_reads
    df['percentage_run_reads']=percentage_run_reads
    return df

def getMeta(meta, pathogens, pathogens_reduced,biofire):
    df=pd.read_csv(meta)
    # strip whitespace from values
    df['pathogen 1'] = df['pathogen 1'].str.strip()

    biofire_pathogens=open(biofire).read().splitlines()

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

    ## Create meta data biofire output
    metaDF=df.copy()
    # add Negavtive control where pathogen 1 is empty
    metaDF['pathogen 1']=np.where(metaDF['pathogen 1'].isnull(),'Negative control',metaDF['pathogen 1'])
    metaDF=metaDF.melt(id_vars=['Run','barcode','sample_name','seq_name'],
                       value_vars=['pathogen 1','pathogen 2','pathogen 3'],
                       var_name='pathogen number',value_name='pathogen')
    metaDF['pathogen reduced']=metaDF['pathogen'].map(path_dict_reduced)
    metaDF.dropna(subset=['pathogen'],inplace=True)
    metaDF['Biofire positive']=np.where(metaDF['pathogen']=='Negative control',0,1)
  
    # add all species not found as biofire negative
    #biofire_pathogens=df3['pathogen_reduced'].unique()
    #biofire_pathogens=list(set(biofire_pathogens)-set(positive_controls))
    runs=metaDF['Run'].unique()
    additonal_rows=[]
    for run in runs:
        barcodes=metaDF[metaDF['Run']==run]['barcode'].unique()
        for barcode in barcodes:
            seq_name=metaDF[(metaDF['Run']==run) & (metaDF['barcode']==barcode)]['seq_name'].unique()[0]
            pathogens=metaDF[(metaDF['Run']==run) & (metaDF['barcode']==barcode)]['pathogen reduced'].unique()
            missed_pathogens=list(set(biofire_pathogens) - set(pathogens))
            for p in missed_pathogens:
                d={'Run':run,'barcode':barcode,'seq_name':seq_name,'pathogen number':None,'pathogen':p,'pathogen reduced':p,'Biofire positive':0}
                additonal_rows.append(d)
    if len(additonal_rows)>0:
        ar_df=pd.DataFrame(additonal_rows)
        metaDF=pd.concat([metaDF,ar_df])

    metaDF.sort_values(by=['Run','barcode','pathogen'],inplace=True)

    metaDF=metaDF[metaDF['pathogen']!='Negative control']
    metaDF=metaDF[metaDF['pathogen']!='SARS_coronavirus_Tor2']
    metaDF=metaDF[metaDF['pathogen']!='orthoreovirus']
    #keep_runs=['expt10_03072024', 'expt11_150824','expt10A_17072024']
    #metaDF=metaDF[metaDF['Run'].isin(keep_runs)]
    
    cols=['Run','barcode','seq_name','pathogen reduced','Biofire positive']
    metaDF=metaDF[cols]
    metaDF.rename(columns={'pathogen reduced':'pathogen'},inplace=True)

    print(metaDF)
    metaDF.to_csv('metaDF.csv', index=False)

    metaDFbiofire_only=metaDF[metaDF['pathogen'].isin(biofire_pathogens)]
    metaDFbiofire_only.to_csv('metaDF_biofire_only.csv', index=False)
    
    return metaDFbiofire_only

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
        dfB=dfB[dfB['chrom']!='unmapped']
        if len(dfB)==0:
            continue
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
                   'pcTP':0, 'TP':0, 'FP':0, 'FN':1}
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
        data['pathogen']=data['chrom'].map(path_dict)
        data['pathogen_reduced']=data['pathogen'].map(path_dict_reduced)
        data['pathogen_reduced']=np.where(data['pathogen_reduced'].isnull(),data['chrom'],data['pathogen_reduced'])

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

def getAdditionalMetrics(fpr, tpr, thresholds, pathogen, metric):
    df=pd.DataFrame({'fpr':fpr, 'tpr':tpr, 'thresholds':thresholds})
    df['sensitivity']=df['tpr']
    df['specificity']=1-df['fpr']
    df['1-specificity']=df['fpr']
    df['pathogen']=pathogen
    df['metric']=metric
    #df['accuracy']=(df['tpr']+df['specificity'])/2
    df['Youden J statistic']=df['sensitivity']+df['specificity']-1
    df['F1 score']=(2*df['sensitivity']*df['specificity'])/(df['sensitivity']+df['specificity'])
    df.to_csv(f'additional_stats/{pathogen}_{metric}_thresholds.csv')


def plotROC(df, pathogen, metric):
    '''Plot the ROC curve'''
    print(pathogen, metric)
    if pathogen!='All_pathogens':
        df2=df[df['pathogen']==pathogen]
    else:
        df2=df
    df2.fillna(0,inplace=True)
    #X=df2[['sample num reads','sample reads percent of type', 'sample reads percent of run']]
    X=df2[metric]
    y=df2['Biofire positive']

    # check if y contail only one class
    if len(y.unique())==1 or len(X)<2:
        print(y)
        return
    
    # show distribution of the data
    #g=df.groupby(['pathogen','Biofire positive'])[metric].count()
    #g.to_csv(f'sample_counts/{pathogen}_{metric}_counts.csv')
    #print(g)

    # ROC curve from sklearn
    fpr, tpr, thresholds = roc_curve(y, X)
    getAdditionalMetrics(fpr, tpr, thresholds, pathogen, metric)
    # Compute the AUC (area under the curve)
    roc_auc = auc(fpr, tpr)
    #print(thresholds)

    #create ROC curve
    plt.plot(fpr,tpr, label="AUC="+str(roc_auc))
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.legend(loc=4)
    # add diagonal line
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    metric_=metric.replace(' ','_')
    plt.savefig(f'rocs/{pathogen}_roc_{metric_}.pdf')
    plt.clf()

def box_plots(df):
    '''Melt wide to long for metrics'''
    metrics=['meanDepth', 'meanDepth_trunc5', 'meanDepth_trunc10','AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc', 
          'Cov1', 'Cov3', 'Cov5', 'Cov10', 
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc',
          'sample num reads','total run reads mapped', 'total run reads inc unmapped', 
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample']
    
    df2=df.melt(id_vars=['Run','seq_name','barcode','pathogen','Biofire positive'],
                value_vars=metrics,
                var_name='metric',value_name='value')
    
    df2=df2[~df2['metric'].isin(['total run reads mapped', 'total run reads inc unmapped'])]
    # facetted box plot with each metric per facet and biofire positive or negative on x axis
    g=sns.FacetGrid(df2, col='metric', col_wrap=4, sharey=False)
    g.map(sns.boxplot, 'Biofire positive', 'value', order=[0,1])
    g.set_titles("{col_name}")
    # set all y-axis to start at 0 
    for ax in g.axes.flat:
        ax.set_ylim(0)
        #ax.set_yscale('log')

    #g.set_xticklabels(rotation=90)
    df2.to_csv('box_plots.csv')
    plt.tight_layout()
    plt.savefig('box_plots.pdf')

    # scatter plots
    plot_metrics=['sample num reads','meanDepth', 'AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc',  
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc']
    df3=df[plot_metrics]
    g2 = sns.PairGrid(df3)
    g2.map(sns.scatterplot)
    plt.tight_layout()
    plt.savefig('scatter_plots.pdf')
    plt.clf()



def main(args):
    df=getDataFrame(args.input)
    metaDFbiofire_only=getMeta(args.meta, args.pathogens, args.pathogens_reduced, args.biofire)
    df.to_csv('all_results.csv')
    # remove duplicates for pathogen_reduced and keep row with highest sample num reads
    df=df.sort_values(by=['batch','Sample name','pathogen_reduced','sample num reads'],ascending=False)
    df.drop_duplicates(subset=['batch','Sample name','pathogen_reduced'],inplace=True,keep='first')
    df.to_csv('all_results_no_duplicates.csv')

    # remove _sup from batch
    df['Run']=df['batch'].str.replace('_sup','')

    # merge with metaDFbiofire_only
    df.drop(columns=['pathogen'],inplace=True)
    df2=metaDFbiofire_only.merge(df,left_on=['Run','barcode','pathogen'],
                                right_on=['Run','barcode','pathogen_reduced'],
                                how='left')
    df2.drop(columns=['pathogen_reduced','run','Sample name','batch','chrom'],inplace=True)

    # make box plots
    box_plots(df2)

    df2.fillna(0,inplace=True)
    df2.to_csv('biofire_results_merged.csv', index=False)

    # plot ROC curves
    metrics=['sample num reads','meanDepth', 'AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc',  
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc']
    for metric in metrics:
        plotROC(df2, 'All_pathogens', metric)

    #print(metaDict)
    #dfs=[]
    #thresholds=[0.0, 0.1, 0.3, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 2.4, 3.3, 4.2, 5.1, 6.0, 7.0, 8.0, 9.0, 10.0, 20, 30, 40, 50, 60, 70, 80, 90, 99]
    #thresholds=[0.0]
    #for p in thresholds:
    #    for r in thresholds:
    #        print(p,r)
    #        data, df2, sensSpec=checkSensSpec(df, metaDict, metaDF, path_dict, pathogens_reduced, num_pathogens, percentage_type_reads=p, percentage_run_reads=r)
    #        data.to_csv(f'{args.output}_{p}_{r}.csv')
    #        if sensSpec is not None:
    #            dfs.append(sensSpec)
    #sensSpec=pd.concat(dfs)
    #sensSpec.to_csv(f'{args.output}_sensSpec.csv')
    
    

    #pathogens=data['pathogen_reduced'].unique()
    #for pathogen in pathogens:
    #    plotROC(data, pathogen, metric)


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
    parser.add_argument('-bf', '--biofire', required=True,
                        help='list of biofire pathogens')
    parser.add_argument('-o', '--output', required=True,
                         help='Output file')
    args = parser.parse_args()

    main(args)
    #df=pd.read_csv(args.output)
    #plotROC(df)

    