#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import log_loss 
from sklearn.metrics import accuracy_score
from sklearn import metrics
from scipy.stats import spearmanr
import os
import joblib
from patsy import dmatrices
import statsmodels.api as sm 

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

def getMeta(meta, pathogens, pathogens_reduced,biofire,keep_runs):
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
    metaDF=metaDF[metaDF['Run'].isin(keep_runs)]
    
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
    '''Calculate sensitivity, specificity, Youden J statistic, F1 score and accuracy and return the top threshold'''
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
    df.sort_values(by='Youden J statistic',ascending=False,inplace=True)
    top_threshold=df.head(1)
    return top_threshold

def getPPV_NPV_MLR(thresh_df, df, pathogen, metric):
    '''Calculate the positive predictive value and negative predictive value'''
    df['TP']=np.where((df['gold_standard']==1) & (df['logit_pred']==1),1,0)
    df['FP']=np.where((df['gold_standard']==0) & (df['logit_pred']==1),1,0)
    df['TN']=np.where((df['gold_standard']==0) & (df['logit_pred']==0),1,0)
    df['FN']=np.where((df['gold_standard']==1) & (df['logit_pred']==0),1,0)
    TP=df['TP'].sum()
    FP=df['FP'].sum()
    TN=df['TN'].sum()
    FN=df['FN'].sum()
    PPV=TP/(TP+FP)
    NPV=TN/(TN+FN)
    thresh_df['PPV']=PPV
    thresh_df['NPV']=NPV
    print(f'PPV: {PPV}, NPV: {NPV}')
    thresh_df.to_csv(f'additional_stats/{pathogen}_{metric}_predictions.csv')
    return df, thresh_df

def getPPV_NPV(thresh_df, df, pathogen, metric):
    '''Calculate the positive predictive value and negative predictive value'''
    print(thresh_df)
    threshold=thresh_df['thresholds'].values[0]
    df['prediction']=np.where(df[metric]>threshold,1,0)
    df['TP']=np.where((df['gold_standard']==1) & (df['prediction']==1),1,0)
    df['FP']=np.where((df['gold_standard']==0) & (df['prediction']==1),1,0)
    df['TN']=np.where((df['gold_standard']==0) & (df['prediction']==0),1,0)
    df['FN']=np.where((df['gold_standard']==1) & (df['prediction']==0),1,0)
    TP=df['TP'].sum()
    FP=df['FP'].sum()
    TN=df['TN'].sum()
    FN=df['FN'].sum()
    PPV=TP/(TP+FP)
    NPV=TN/(TN+FN)
    thresh_df['PPV']=PPV
    thresh_df['NPV']=NPV
    print(f'PPV: {PPV}, NPV: {NPV}')
    thresh_df.to_csv(f'additional_stats/{pathogen}_{metric}_predictions.csv')
    return thresh_df

def plotROC(df, pathogen, metric, remove_no_data=True):
    '''Plot the ROC curve'''
    print(pathogen, metric)
    # remove rows where sample num reads is 0 and biofire is positive for that sample
    if remove_no_data:
        df=df[~((df['sample num reads']==0) & (df['Biofire positive']==1))]
    if pathogen!='All_pathogens':
        df2=df[df['pathogen']==pathogen]
    else:
        df2=df
    df2.fillna(0,inplace=True)
    #X=df2[['sample num reads','sample reads percent of type', 'sample reads percent of run']]
    X=df2[metric]
    y=df2['gold_standard']

    # check if y contail only one class
    if len(y.unique())==1 or len(X)<2:
        print(y)
        return
    
    # ROC curve from sklearn
    fpr, tpr, thresholds = roc_curve(y, X)
    threshold=getAdditionalMetrics(fpr, tpr, thresholds, pathogen, metric)
    thresh_df=getPPV_NPV(threshold, df2, pathogen, metric)
    thresh_df['remove_no_data']=remove_no_data
    # Compute the AUC (area under the curve)
    roc_auc = auc(fpr, tpr)
    #print(thresholds)
    thresh_df['AUC']=roc_auc

    #create ROC curve
    plt.plot(fpr,tpr, label="AUC="+str(roc_auc))
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    metric_=metric.replace(' ','_')

    if remove_no_data:
        plt.savefig(f'rocs/{pathogen}_roc_{metric_}_no_zero.pdf')
    else:
        plt.savefig(f'rocs/{pathogen}_roc_{metric_}.pdf')
    plt.clf()
    plt.close()
    return thresh_df

def plotROCs(df2):
    metrics=['sample num reads','meanDepth', 'AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc',  
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc',
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'median_read_length', 'median_aligned_length', 'mean_read_length', 'mean_aligned_length',
          'Sample_num_reads_200', 'Sample_num_reads_300', 'Sample_num_reads_400']
    thresh_dfs=[]
    for metric in metrics:
        thresh_df=plotROC(df2, 'All_pathogens', metric, remove_no_data=False)
        if thresh_df is not None:
            thresh_dfs.append(thresh_df)

        # again but remove samples with 0 reads and biofire positive
        thresh_df=plotROC(df2, 'All_pathogens', metric, remove_no_data=True)
        if thresh_df is not None:
            thresh_dfs.append(thresh_df)
    thresh_df=pd.concat(thresh_dfs)
    thresh_df.to_csv('all_pathogens_thresholds.csv')

def reg_coef(x,y,label=None,color=None,hue=None,**kwargs):
    ax = plt.gca()
    # remove nan values
    x2=[]
    y2=[]
    #print('X:  ',x.name)
    #print('Y:  ',y.name)
    #import sys
    #sys.exit()
    for ix, iy in zip(x, y):
        if np.isnan(ix) or np.isnan(iy):
            continue
        else:
            x2.append(ix)
            y2.append(iy)
    p = spearmanr(x2,y2)
    #print(list(x2),list(y2), p.correlation, p.pvalue)
    ax.annotate('x: {}\ny: {}\n \u03C1={:.2f}'.format(x.name, y.name,p.correlation ), xy=(0.5,0.5), xycoords='axes fraction', ha='center')
    #ax.set_axis_off()

    # cirlce plots
    
    #print(c_size)
    x_max = max(x2)
    y_max = max(y2)
    c_size = (len(x2)/len(x))*x_max
    if p.correlation>0:
        c_colour='blue'
    elif p.correlation==0:
        c_colour='black'
    else:
        c_colour='red'
    #Drawing_uncolored_circle = plt.Circle((x_max/2, y_max/2), radius=c_size, fill=False, color=c_colour, lw=2, alpha=0.5)
    #c_size2 = len(x2)/len(x)
    #Drawing_uncolored_circle = plt.Circle((0.5, 0.5), radius=c_size2, fill=False, color='black', lw=2, alpha=0.5)
    #ax.add_patch(Drawing_uncolored_circle)
    #ax.set(xlim=(0,1),ylim=(0, 1))
    
def reg_coef_inside(x,y,label=None,color=None,hue=None,**kwargs):
    ax = plt.gca()
    # remove nan values
    x2=[]
    y2=[]
    #print('X:  ',x.name)
    #print('Y:  ',y.name)
    #import sys
    #sys.exit()
    for ix, iy in zip(x, y):
        if np.isnan(ix) or np.isnan(iy):
            continue
        else:
            x2.append(ix)
            y2.append(iy)
    p = spearmanr(x2,y2)
    #print(list(x2),list(y2), p.correlation, p.pvalue)
    ax.annotate('P={:.2f}'.format(p.correlation ), xy=(0.7,0.1), xycoords='axes fraction', ha='center')
    ax.set_axis_off()

def box_plots(df):
    '''Melt wide to long for metrics'''
    metrics=['meanDepth', 'meanDepth_trunc5', 'meanDepth_trunc10','AuG','AuG_trunc5','AuG_trunc10',
          'bases', 'bases_perc', 
          'Cov1', 'Cov3', 'Cov5', 'Cov10', 
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc',
          'sample num reads','total run reads mapped', 'total run reads inc unmapped', 
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'median_read_length', 'median_aligned_length', 'mean_read_length', 'mean_aligned_length',
          'Sample_num_reads_200', 'Sample_num_reads_300', 'Sample_num_reads_400']
    
    df2=df.melt(id_vars=['Run','seq_name','barcode','pathogen','Biofire positive','gold_standard'],
                value_vars=metrics,
                var_name='metric',value_name='value')
    
    df2=df2[~df2['metric'].isin(['total run reads mapped', 'total run reads inc unmapped'])]
    # facetted box plot with each metric per facet and biofire positive or negative on x axis
    g=sns.FacetGrid(df2, col='metric', col_wrap=6, sharey=False)
    g.map(sns.boxplot, 'gold_standard', 'value', order=[0,1])
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
          'Cov1_perc','Cov3_perc','Cov5_perc','Cov10_perc',
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'median_read_length', 'median_aligned_length', 'mean_read_length', 'mean_aligned_length',
          'Sample_num_reads_200', 'Sample_num_reads_300', 'Sample_num_reads_400']
    df3=df[plot_metrics]
    g2 = sns.PairGrid(df3)
    g2.map_upper(reg_coef, hue=None)
    #g2.map_diag(sns.kdeplot)
    g2.map_lower(sns.scatterplot)
    #g2.map_lower(sns.regplot,  scatter_kws = {"color": "black", "alpha": 0.5},
    #        line_kws = {"color": "red"})
    #g2.map_lower(reg_coef_inside, hue=None)
    
    # make both axis log scale
    #for ax in g2.axes.flat:
    #if ax.get_xlabel() in log_columns:
    #    ax.set(xscale="log")

    plt.tight_layout()
    plt.savefig('scatter_plots.pdf')
    plt.clf()

    ## create matrix
    df4=df3.corr(method='spearman')
    df4.to_csv('spearman_correlation_matrix.csv')

    # melt matrix
    df5=df4.reset_index().melt(id_vars='index')
    df5.to_csv('spearman_correlation_matrix_melted.csv')

def adjust_gold_standard(df):
    '''Create a new column with the gold standard
    change the Influnza subtypes gold standard to 1 if Influnza A positive and sequencing maps to Influenza A'''
    df['gold_standard']=df['Biofire positive']
    
    # FLU_A_POS is 1 if biofire positive for Influenza A for the run and barcode
    df['condition'] = (df['Biofire positive'] == 1) & (df['pathogen'] == 'Influenza A')
    df['FLU_A_POS'] = df.groupby(['Run', 'barcode'])['condition'].transform('any')

    flu_A_options=['Influnza A','Influenza A/H1-2009','Influenza A/H3','Influenza A/H1']

    filtered_df = df[df['pathogen'].isin(flu_A_options)]

    # Step 2: Group by 'Run' and 'barcode' and find the index of the row with the highest 'num_reads'
    filtered_df['TOP_FLU_A'] = filtered_df.groupby(['Run', 'barcode'])['sample num reads'].transform('idxmax')

    # Step 3: Set the 'top_flu' column to 1 for the row with the highest 'num_reads', and 0 for others
    df['TOP_FLU_A'] = df.index.isin(filtered_df['TOP_FLU_A'])
    
    df['gold_standard']=np.where( ( df['pathogen_reduced'].isin(flu_A_options) ) & ( df['FLU_A_POS']==1) & (df['TOP_FLU_A']==1), 1, df['gold_standard'] )
    return df

def plot_read_length_hist(input, pathogens, pathogens_reduced):
    dfs=[]
    for folder in input:
        run=folder.split('/')[-1]
        barcodes=os.listdir(f'{folder}/alignment_info/')
        barcodes=[b for b in barcodes if b.endswith('_all_alignments.csv')]
        for barcode in barcodes:
            df=pd.read_csv(f'{folder}/alignment_info/{barcode}')
            df['Run']=run
            dfs.append(df)
            b=barcode.split('_')[0]
            df['barcode']=b
            dfs.append(df)

    df=pd.concat(dfs)

    # reduce ref to pathogens
    df2=pd.read_csv(pathogens)
    
    path_dict={}
    for i,r in df2.iterrows():
        path_dict.setdefault(r['pathogen'],[]).append(r['reference'])

    # reverse path_dict
    path_dict_rev={v:k for k,values in path_dict.items() for v in values}
    df['Pathogen']=df['ref'].map(path_dict_rev)

    df3=pd.read_csv(pathogens_reduced)
    df3.set_index('pathogen',inplace=True)
    path_dict_reduced=df3.to_dict()['pathogen_reduced']

    df['Pathogen_reduced']=df['Pathogen'].map(path_dict_reduced)

    # plot histogram of read length by species
    g=sns.FacetGrid(df, col='Run', row='Pathogen_reduced', sharey=False)
    g.map(sns.histplot, 'read_len')
    plt.savefig('plots/read_length_hist.pdf')
    plt.clf()

    # plot histogram of align_len length by species
    g=sns.FacetGrid(df, col='Run', row='Pathogen_reduced', sharey=False)
    g.map(sns.histplot, 'align_len')
    g.set_titles("{row_name}|{col_name}")
    plt.savefig('plots/align_length_hist.pdf')
    plt.clf()

    # box plot of read length by species
    fig, ax = plt.subplots(figsize=(16, 9))
    sns.boxplot(data=df, x='Pathogen_reduced', y='read_len')
    plt.xticks(rotation=90)
    # set y axis range from 0 to 1000
    plt.ylim(0,1000)
    plt.tight_layout()

    # set dimensions of plot
    plt.savefig('plots/read_length_boxplot.pdf')
    plt.clf()

    # box plot of align length by species
    fig, ax = plt.subplots(figsize=(16, 9))
    sns.boxplot(data=df, x='Pathogen_reduced', y='align_len')
    plt.xticks(rotation=90)
    plt.ylim(0,1000)
    plt.tight_layout()
    
    plt.savefig('plots/align_length_boxplot.pdf')
    plt.clf()

def multivariate_logistic_regression(df, remove_no_data=False):
    metrics=['AuG_trunc10', 'sample num reads', 'Cov1_perc',
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 
          'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'mean_read_length','mean_aligned_length']
    
    if remove_no_data:
        suffix='_no_zero'
        df=df[~((df['sample num reads']==0) & (df['gold_standard']==1))]
    else:
        suffix=''

    X=df[metrics]
    y=df['gold_standard']
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)
    logreg = LogisticRegression(max_iter=1000)
    logreg.fit(X, y)
    y_pred = logreg.predict(X)
    y_probs = logreg.predict_proba(X)[:,1]
    df['logit_probs']=y_probs
    df['logit_pred']=y_pred

    # save model
    joblib.dump(logreg, f'logit_model{suffix}.pkl')

    #print('Accuracy of logistic regression classifier on test set: {:.2f}'.format(logreg.score(X_test, y_test)))

    # confusion matrix
    #confusion_matrix = confusion_matrix(y, y_pred)
    #print(confusion_matrix)

    # calculate AIC
    n = len(y)
    k = len(metrics)+1
    ll = log_loss(y, y_probs)
    aic = 2*k - 2*ll
    print(f'AIC: {aic}')


    # classification report
    print(classification_report(y, y_pred))

    # ROC curve
    logit_roc_auc = roc_auc_score(y, logreg.predict(X))
    fpr, tpr, thresholds = roc_curve(y, logreg.predict_proba(X)[:,1])
    plt.figure()
    plt.plot(fpr, tpr, label='Logistic Regression (area = %0.2f)' % logit_roc_auc)
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver operating characteristic {suffix}')
    plt.legend(loc="lower right")
    plt.savefig(f'rocs/Log_ROC{suffix}.pdf')
    plt.clf()
    plt.close()

    # make data frame for fpr, tpr, thresholds
    df2=pd.DataFrame({'fpr':fpr, 'tpr':tpr})
    #for i, threshold in enumerate(thresholds):
    df2[f'thresholds'] = thresholds
    df2['sensitivity']=df2['tpr']
    df2['specificity']=1-df2['fpr']
    df2['1-specificity']=df2['fpr']
    df2['Youden J statistic']=df2['sensitivity']+df2['specificity']-1
    df2['F1 score']=(2*df2['sensitivity']*df2['specificity'])/(df2['sensitivity']+df2['specificity'])
    df2.to_csv(f'logit_thresholds{suffix}.csv')
    #print(thresholds)

    # print feature importance
    # Get feature importance
    feature_importance = np.abs(logreg.coef_[0])
    feature_names = X.columns #if hasattr(X, 'columns') else [f"Feature {i}" for i in range(X.shape[1])]
    # Display sorted feature importance
    sorted_indices = np.argsort(feature_importance)[::-1]
    fi_list=[]
    for idx in sorted_indices:
        fi_list.append({'feature': feature_names[idx], 'importance':feature_importance[idx]})
        print(f"{feature_names[idx]}: {feature_importance[idx]}")
    fi_df=pd.DataFrame(fi_list)
    fi_df.to_csv(f'logit_feature_importance{suffix}.csv')
    # bar plot of feature importance
    plt.figure(figsize=(10, 6))
    sns.barplot(x=fi_df['importance'], y=fi_df['feature'])
    plt.title('Feature Importance')
    plt.xlabel('Importance')
    plt.ylabel('Feature')
    #plt.xscale('log')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'features/logit_feature_importance{suffix}.pdf')    
    plt.clf()
    plt.close()

    X=X.to_numpy()
    for feature_index in range(X.shape[1]):
        #feature_index = 0  # Feature of interest (index of the feature)
        feature_values = X[:, feature_index]

        name=metrics[feature_index]
        # Create a scatter plot to show the relationship
        plt.figure(figsize=(8, 6))
        plt.scatter(feature_values, y_probs, alpha=0.5, s=10, c=y, cmap='coolwarm')
        plt.title(f"Effect of Feature {name} on Predicted Probability")
        plt.xlabel(f"Feature {name} Value")
        plt.ylabel("Predicted Probability (P(y=1))")
        plt.colorbar(label='True Labels')  # Color represents true labels (y_test)
        plt.grid(True)
        plt.savefig(f'features/logit_feature_{name}{suffix}.pdf')
        plt.clf()
        plt.close()

        # Sort feature values and predicted probabilities
        sorted_indices = np.argsort(feature_values)
        sorted_feature_values = feature_values[sorted_indices]
        sorted_probs = y_probs[sorted_indices]

        # Plot a smoother curve
        plt.figure(figsize=(8, 6))
        plt.plot(sorted_feature_values, sorted_probs, color='blue', linewidth=2)
        plt.title(f"Effect of Feature {name} on Predicted Probability (Smoothed)")
        plt.xlabel(f"Feature {name} Value")
        plt.ylabel("Predicted Probability (P(y=1))")
        plt.grid(True)
        plt.savefig(f'features/logit_feature_{name}_smoothed{suffix}.pdf')
        plt.clf()
        plt.close()

    ## PPV and NPV
    pathogen=f'All_pathogens{suffix}'
    metric='logit_probs'
    threshold=getAdditionalMetrics(fpr, tpr, thresholds, pathogen, metric)
    df,thresh_df=getPPV_NPV_MLR(threshold, df, pathogen, metric)
    #df.to_csv('logit_predictions.csv')
    df.to_csv(f'logit_predictions{suffix}.csv')
    thresh_df['remove_no_data']='remove_no_data'  
    thresh_df['AUC']=logit_roc_auc
    thresh_df.to_csv(f'logit_thresholds_PPV_NPV{suffix}.csv')

    # plot preds and probs by TP, FP, TN, FN
    cols=['logit_probs','logit_pred','gold_standard', 'pathogen', 'TP', 'FP', 'TN', 'FN']
    df3=df[cols]
    # melt to long
    df3=df3.melt(id_vars=['logit_probs','logit_pred','gold_standard', 'pathogen'],
                               value_vars=['TP','FP','TN','FN'],
                               var_name='Test result type',value_name='test result')
    df3=df3[df3['test result']==1]
    g2=sns.boxplot(data=df3, x='Test result type', y='logit_probs', color='grey', fill=False)
    # add dots to boxplot
    sns.stripplot(data=df3, x='Test result type', y='logit_probs', hue='pathogen', alpha=0.5)
    # move legend outside of plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(f'plots/logit_predictions_boxplot{suffix}.pdf')
    plt.savefig(f'plots/logit_predictions_boxplot{suffix}.svg')
    plt.clf()
    plt.close()

def statsmodels_logistic_regression(df, remove_no_data=False):
    metrics=['AuG_trunc10', 'bases', 'sample_num_reads', 'Cov1_perc',
          'Sample_reads_percent_of_run', 'Sample_reads_percent_of_refs', 
          'Sample_reads_percent_of_type_run', 'Sample_reads_percent_of_type_sample',
          'mean_read_length','mean_aligned_length']
    
    if remove_no_data:
        suffix='_no_zero'
        df=df[~((df['sample num reads']==0) & (df['gold_standard']==1))]
    else:
        suffix=''
    df.rename(columns={'sample num reads':'sample_num_reads'},inplace=True)
    #X=df[metrics]
    #y=df['gold_standard']

    metrics_str=' + '.join(metrics)
    print(metrics_str)
    # AuG_trunc10 + bases + sample_num_reads + Cov1_perc + Sample_reads_percent_of_run + Sample_reads_percent_of_refs + Sample_reads_percent_of_type_run + Sample_reads_percent_of_type_sample + mean_read_length + mean_aligned_length
    y, X = dmatrices( f'gold_standard ~ AuG_trunc10 + sample_num_reads + Cov1_perc + Sample_reads_percent_of_run + Sample_reads_percent_of_refs + Sample_reads_percent_of_type_run + Sample_reads_percent_of_type_sample + mean_read_length + mean_aligned_length',
                      data=df, return_type='dataframe')
    mod = sm.Logit(y, X)
    res = mod.fit(maxiter=1000)
    print('AIC: ',res.aic)
    print(res.summary())

    # plot odds ratios
    fig, ax = plt.subplots(figsize=(10, 6))
    conf = res.conf_int()
    conf['OR'] = res.params
    conf.columns = ['2.5%', '97.5%', 'OR']
    conf = np.exp(conf)
    conf = conf.sort_values('OR', ascending=True)
    conf = conf.dropna()
    conf['OR'].plot(kind='barh', color='blue', ax=ax)
    ax.set_title('Odds Ratios')
    plt.tight_layout()
    plt.savefig(f'SM/features/logit_odds_ratios{suffix}.pdf')
    plt.clf()
    plt.close()

    # plot ROC
    y_probs = res.predict(X)
    fpr, tpr, thresholds = roc_curve(y, y_probs)
    logit_roc_auc = roc_auc_score(y, y_probs)
    plt.figure()
    plt.plot(fpr, tpr, label='Logistic Regression (area = %0.2f)' % logit_roc_auc)
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver operating characteristic {suffix}')
    plt.legend(loc="lower right")
    plt.savefig(f'SM/rocs/Log_ROC{suffix}.pdf')
    plt.clf()
    plt.close()

    return res, mod


def main(args):
    df=getDataFrame(args.input)
    keep_runs=[run.split('/')[-1] for run in args.input]
    metaDFbiofire_only=getMeta(args.meta, args.pathogens, args.pathogens_reduced, args.biofire, keep_runs)
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
    df2=adjust_gold_standard(df2)
    df2.drop(columns=['pathogen_reduced','run','Sample name','batch','chrom'],inplace=True)

    # make box plots
    #box_plots(df2)

    df2.fillna(0,inplace=True)
    df2.to_csv('biofire_results_merged.csv', index=False)

    # plot ROC curves
    #plotROCs(df2)

    # Multivariate logistic regression
    multivariate_logistic_regression(df2)
    multivariate_logistic_regression(df2, remove_no_data=True)

    res,mod=statsmodels_logistic_regression(df2)
    #return res,mod, df2

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

    #res, mod2, df2 = main(args)
    main(args)
    #plot_read_length_hist(args.input, args.pathogens, args.pathogens_reduced)
    #df=pd.read_csv(args.output)
    #plotROC(df)

    