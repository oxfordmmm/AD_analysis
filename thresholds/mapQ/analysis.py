#!/usr/bin/env python3
import argparse
import pandas as pd
import pysam
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, roc_auc_score


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

def readBamFiles(bamFile):
    sam=pysam.AlignmentFile(bamFile, 'rb')
    contig_scores=[]
    for read in sam.fetch():
        # get contig name and mapQ and added to a list
        contig_scores.append({'chrom':read.reference_name, 'readID':read.query_name,
                              'edit_distance':read.get_tag('NM'), 
                              'aln_length':read.query_alignment_length,
                               'mapQ': read.mapping_quality})
        
    df=pd.DataFrame(contig_scores)
    return df


def getBams(runs: list):
    bams = []
    for run in runs:
        barcodes=os.listdir(f'{run}/bams/minimap2_unfiltered')
        bamFiles=[barcode for barcode in barcodes if barcode.endswith('minimap.bam')]
        for bam in bamFiles:
            df=readBamFiles(f'{run}/bams/minimap2_unfiltered/{bam}')
            df['Sample name']=bam.split('/')[-1]
            df['run']=run.split('/')[-1]
            bams.append(df)

        #bams.append(pysam.AlignmentFile(run, 'rb'))
    df=pd.concat(bams)
    df['Sample name']=df['Sample name'].str.replace('SQK-RBK114-96_barcode','')
    df['Sample name']=df['Sample name'].str.replace('_minimap.bam','')
    df['edit_percentage']=(df['edit_distance']/df['aln_length'])*100

    return df

def getPathogens(pathogens, pathogens_reduced, meta):
    df2=pd.read_csv(pathogens)
    
    path_dict={}
    for i,r in df2.iterrows():
        path_dict.setdefault(r['pathogen'],[]).append(r['reference'])

    # reverse path_dict
    path_dict_rev={v:k for k,values in path_dict.items() for v in values}

    df3=pd.read_csv(pathogens_reduced)
    df3.set_index('pathogen',inplace=True)
    path_dict_reduced=df3.to_dict()['pathogen_reduced']


    df=pd.read_csv(meta)
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

    #print(metaDF)
    metaDF.to_csv('metaDF.csv', index=False)

    metaDFbiofire_only=metaDF[metaDF['pathogen'].isin(biofire_pathogens)]
    
     
    
    return path_dict_rev, path_dict_reduced, metaDFbiofire_only

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
    
    return metaDFbiofire_only, path_dict_rev, path_dict_reduced

def adjust_gold_standard(df):
    '''Create a new column with the gold standard
    change the Influnza subtypes gold standard to 1 if Influnza A positive and sequencing maps to Influenza A'''
    df['gold_standard']=df['Biofire positive']
    
    # FLU_A_POS is 1 if biofire positive for Influenza A for the run and barcode
    df['condition'] = (df['Biofire positive'] == 1) & (df['pathogen_y'] == 'Influenza A')
    df['FLU_A_POS'] = df.groupby(['Run', 'barcode'])['condition'].transform('any')

    flu_A_options=['Influnza A','Influenza A/H1-2009','Influenza A/H3','Influenza A/H1']

    filtered_df = df[df['pathogen_y'].isin(flu_A_options)]

    # Step 2: Group by 'Run' and 'barcode' and find the index of the row with the highest 'num_reads'
    filtered_df['TOP_FLU_A'] = filtered_df.groupby(['Run', 'barcode'])['sample num reads'].transform('idxmax')

    # Step 3: Set the 'top_flu' column to 1 for the row with the highest 'num_reads', and 0 for others
    df['TOP_FLU_A'] = df.index.isin(filtered_df['TOP_FLU_A'])
    
    df['gold_standard']=np.where( ( df['pathogen_reduced'].isin(flu_A_options) ) & ( df['FLU_A_POS']==1) & (df['TOP_FLU_A']==1), 1, df['gold_standard'] )
    return df



def logistic_regression(df):
    # filter out Flu
    flu_A_options=['Influnza A','Influenza A/H1-2009','Influenza A/H3','Influenza A/H1']
    df=df.query('pathogen_reduced not in @flu_A_options')


    X=df[['mapQ']]
    y=df['True species hit']
    #clf = LogisticRegression(random_state=0).fit(X, y)
    #y_pred=clf.predict(X)
    #y_pred_prob=clf.predict_proba(X)[:,1]
    fpr, tpr, thresholds = roc_curve(y, X)
    roc_auc = auc(fpr, tpr)
    print(f'AUC: {roc_auc}')
    # plot ROC curve
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
                lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.savefig('MAPQ_ROC.pdf')

    # get Youden's J statistic
    df2=pd.DataFrame({'fpr':fpr, 'tpr':tpr, 'thresholds':thresholds})
    df2['sensitivity']=df2['tpr']
    df2['specificity']=1-df2['fpr']
    df2['J']=df2['sensitivity']+df2['specificity']-1
    print(df2)
    max_J=df2[df2['J']==df2['J'].max()]
    print(max_J)
    print(f'Youden J statistic: {max_J["J"].values[0]}')
    df2.to_csv('roc.csv', index=False)

    return roc_auc


def plot_mapQ(df):
    #g=sns.swarmplot(data=df, x='pathogen_reduced', y='mapQ', hue='True species hit')
    g1=sns.boxplot(data=df, x='pathogen_reduced', y='mapQ', hue='True species hit', showfliers=False)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('mapQ.pdf')
    plt.close()
    plt.clf()

    # plot just HRE

    df2=df[df['pathogen_reduced']=='Rhinovirus/enterovirus']
    g2=sns.boxplot(data=df2, x='chrom', y='mapQ', hue='True species hit')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('mapQ_HRE.pdf')
    plt.close()
    plt.clf()

    # plot edit ditance over mapQ for each species in a facet plot
    g=sns.FacetGrid(df, col='pathogen_reduced', col_wrap=4, hue='True species hit')
    g.map(sns.scatterplot, 'edit_percentage', 'mapQ')
    # log x axis
    #g.set(xscale='log')
    plt.tight_layout()
    plt.savefig('edit_percentage_mapQ.pdf')
    plt.close()
    plt.clf()


def main(args):
    bams=getBams(args.input)
    meta, pathogens, pathogens_reduced=getMeta(args.meta, args.pathogens, args.pathogens_reduced, args.biofire)
    bams['pathogen']=bams['chrom'].map(pathogens)
    bams['pathogen_reduced']=bams['pathogen'].map(pathogens_reduced)
    bams=bams[bams['Sample name']!='unclassified']
    bams['Sample name']=bams['Sample name'].map(int)

    df=bams.merge(meta, left_on=['run','Sample name', 'pathogen_reduced'], right_on=['Run', 'barcode', 'pathogen'], how='left')
    df['True species hit']=np.where(df['pathogen_reduced'].isin(positive_controls),1,0)
    df['True species hit']=np.where(df['Biofire positive']==1,1,df['True species hit'])
    # remove flu A subtypes
    #flu_A_options=['Influnza A','Influenza A/H1-2009','Influenza A/H3','Influenza A/H1']
    #df=df.query('pathogen_y not in @flu_A_options')
    df.to_csv('df.csv', index=False)
    print(df)

    # subsample to 1000 reads for each 'True species hit'
    bacteria=['Bordetella_parapertussis','Bordetella_pertussis']
    # only keep pathogen_reduced in bacteria and True species hit is 0, plus all True species hit is 1
    #df=df[(df['pathogen_reduced'].isin(bacteria) & (df['True species hit']==0)) | (df['True species hit']==1)]

    # add column for duplicates readID
    df['readID_dup']=df['readID'].duplicated()
    df.sort_values(by=['readID', 'mapQ'], ascending=False, inplace=True)
    df.drop_duplicates(subset=['readID'], keep='first', inplace=True)

    df2=df.groupby(['True species hit', 'pathogen_reduced']).apply(lambda x: x.sample(n=1000, replace=True)).reset_index(drop=True)
    df2.drop_duplicates(subset=['readID','chrom'], keep='first', inplace=True)
    df2.to_csv('df2.csv', index=False)
    print(df2)

    plot_mapQ(df2)
    roc_auc=logistic_regression(df)

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
    args = parser.parse_args()

    #res, mod2, df2 = main(args)
    main(args)